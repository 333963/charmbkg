#!/usr/bin/env python3
"""DNET etap 6 — sieć D0 z torami otoczenia BEZ toru zgodnego z wolnym pionem D* (maska dn6_mask: Δm < DM_MAX, ładunek zgodny; bez prawdy MC).
Cała reszta = D0.tr.H+e04.mc.d5 (ta sama pula, foldy, wagi, DisCo masy λ = 5, 20 epok, ziarna 1–3; kod treningu etapu 5a: dn5_train fold/merge).
Maska nakładana RAZ przy wczytaniu tabel (a_E i p_E wszystkich próbek) → trening, score OOF tabel etapu 1 (w tym pary p_ tagu D*) i tło widzą
te same wejścia. Tagi: D0.tr.H+e04np.mc.d5.s{1,2,3} (bloki np1…np3).

  dn6_train.py fold  --data DIR1/data --env DIR2/parts --out DIR --block npS --fold k [--device cpu|cuda|auto]
  dn6_train.py merge --data … --env … --models DIR/models --folds-dir DIR/folds --out DIR --block npS
  dn6_train.py list | selftest --dir KATALOG   (selftest: dane syntetyczne, 2 foldy × 1 epoka, maska, zapis/odczyt, dn6_apply)"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
DN5 = os.environ.get("DN5") or os.path.join(_Q, "dnet5")
sys.path.insert(0, DN5)
sys.path.insert(0, HERE)
import dn5_lib as L5  # noqa: E402
import dn5_train as T5  # noqa: E402
import dn6_mask as MK  # noqa: E402

FEATS = "H+e04np"


def _np(seed):
    c = L5._v("D0", seed=seed)
    c.update({"feats": FEATS, "tag": f"D0.tr.{FEATS}.mc.{c['dec']}.s{seed}", "mask": {"kind": "dstar_dm", "dm_max": MK.DM_MAX, "sign": MK.SIGN, "mpi": MK.MPI}})
    return c


BLOCKS = {f"np{s}": _np(s) for s in (1, 2, 3)}
L5.BLOCKS.update(BLOCKS)
_load0 = L5.load


def _load(ddir, edir, sp, samples=None):
    S, full = _load0(ddir, edir, sp, samples)
    MK.mask_tables(S, log=lambda s: print(s, flush=True))
    return S, full


L5.load = _load                                   # dn5_train.prepare woła L5.load → maska w treningu i w score OOF


def selftest(d):
    """Mechanika na danych syntetycznych etapu 4 (dn4_train.synth): maska zmienia tokeny, fold × 2 + merge, score w [0, 1], meta z maską,
    dn6_apply ≡ score OOF (te same wiersze), dn4_lib.score_rows na NIEzamaskowanych torach ≠ (kontrola, że maska jest używana)."""
    import json
    import subprocess
    import dn4_train as T4
    import dn4_lib as L4
    ddir, edirs = os.path.join(d, "data"), {"D0": os.path.join(d, "parts"), "Dplus": os.path.join(d, "parts_dp")}
    T4.synth(ddir, edirs, n=3000)
    S0, _ = _load0(ddir, edirs["D0"], "D0")
    S1, _ = L5.load(ddir, edirs["D0"], "D0")
    T0, T1 = S0["fd_nd"], S1["fd_nd"]
    ch = (T0["a_E"] != T1["a_E"]).any((1, 2))
    dm = MK.token_dm(T0["a_E"], T0["a_pt"], T0["a_eta"], T0["a_phi"], T0["a_m"])
    assert ch.sum() == (T1["a_nrm"] > 0).sum(), "licznik maski"
    assert not MK.veto(T1["a_E"], T1["a_Dd"], T1["a_pt"], T1["a_eta"], T1["a_phi"], T1["a_m"]).any(), "po masce nie zostaje żaden tor w oknie"
    print(f"selftest maska: wierszy zmienionych {int(ch.sum())}/{len(ch)}, Δm syntetyczne min {np.nanmin(dm):.3f} mediana {np.nanmedian(dm):.3f}", flush=True)
    # czteropęd: tor współliniowy z D0 o p → Δm z definicji
    E = np.zeros((1, 2, 6), np.float16); E[0, 0] = [0.5, 0.0, 0.0, 1, 1, 10]; E[0, 1] = [0.5, 0.02, -0.01, -1, 1, 10]
    Dd = np.zeros((1, 2, 8), np.float32); Dd[0, 0, 7] = -1
    pt, eta, phi, m = np.array([6.0]), np.array([0.3]), np.array([1.0]), np.array([1.865])
    v = MK.token_dm(E, pt, eta, phi, m)[0]
    def ref(tpt, te, tp):
        pD = np.array([6 * np.cos(1.0), 6 * np.sin(1.0), 6 * np.sinh(0.3)]); ED = np.sqrt(1.865 ** 2 + pD @ pD)
        pt_ = np.array([tpt * np.cos(tp), tpt * np.sin(tp), tpt * np.sinh(te)]); Et = np.sqrt(MK.MPI ** 2 + pt_ @ pt_)
        return np.sqrt((ED + Et) ** 2 - (pD + pt_) @ (pD + pt_)) - 1.865
    r = [ref(0.5, 0.3, 1.0), ref(0.5, 0.3 + float(E[0, 1, 1]), 1.0 + float(E[0, 1, 2]))]
    assert np.allclose(v, r, atol=2e-5), (v, r)
    vt = MK.veto(E, Dd, pt, eta, phi, m, 0.155, 1)[0]
    assert list(vt) == [bool(r[0] < 0.155), False], vt                  # tor 0: ładunek +1 = −q_K → usuwany, tor 1: ładunek = q_K → zostaje
    print(f"selftest czteropęd: Δm {v} (odniesienie {r}), weto {vt}", flush=True)
    com = ["--data", ddir, "--env", edirs["D0"], "--env-dp", edirs["Dplus"], "--folds", "2", "--epochs", "1", "--threads", "2", "--hsc", d]
    out = os.path.join(d, "out")
    env = {**os.environ, "DN_REDO": "1"}
    for k in range(2):
        subprocess.run([sys.executable, "-u", os.path.abspath(__file__), "fold", "--out", out, "--block", "np1", "--fold", str(k)] + com, check=True, env=env)
    subprocess.run([sys.executable, "-u", os.path.abspath(__file__), "merge", "--out", out, "--block", "np1", "--models", os.path.join(out, "models"),
                    "--folds-dir", os.path.join(out, "folds")] + com, check=True, env=env)
    tag = BLOCKS["np1"]["tag"]
    z = np.load(os.path.join(out, "sc", tag + ".npz"))
    meta = json.load(open(os.path.join(out, "meta", tag + ".json")))
    assert meta["mask"]["dm_max"] == MK.DM_MAX and meta["feats"] == FEATS and meta["tag"] == tag, meta.get("mask")
    assert {"fd_nd__a", "fd_nd__p", "mb_nd__a"} <= set(z.files), z.files
    assert all(np.nanmin(z[k]) >= 0 and np.nanmax(z[k]) <= 1 and np.isfinite(z[k]).any() for k in z.files)
    _, ms = L4.load_models(os.path.join(out, "models"), tag)                        # format etapu 4
    T = S1["data_428770"]; U = S0["data_428770"]
    s_m = L4.score_rows("tr", ms, T["a_uid"], T["a_X"], T["a_Dd"], T["a_E"], T["a_eta"], T["a_phi"])
    s_u = L4.score_rows("tr", ms, U["a_uid"], U["a_X"], U["a_Dd"], U["a_E"], U["a_eta"], U["a_phi"])
    d_m = float(np.max(np.abs(s_m - z["data_428770__a"][T["a_idx"]])))
    chd = (T["a_nrm"] > 0)
    d_u = float(np.max(np.abs(s_u - s_m)[chd])) if chd.any() else 0.0
    assert d_m < 1e-6, d_m
    assert float(np.max(np.abs(s_u - s_m)[~chd])) < 1e-6
    print(f"selftest score: OOF ≡ dn4_lib na zamaskowanych torach (max|Δ| {d_m:.2g}); bez maski różni się na {int(chd.sum())} zmienionych wierszach (max|Δ| {d_u:.3g})", flush=True)
    import dn6_apply as A6
    md = os.path.join(d, "models_app"); os.makedirs(md, exist_ok=True)
    if not os.path.exists(os.path.join(md, tag)):
        os.symlink(os.path.join(out, "models", tag), os.path.join(md, tag))
    A6.A4.main(["--sample", "data_428770", "--data", ddir, "--env", edirs["D0"], "--models", md, "--species", "D0", "--out", os.path.join(d, "app.npz"), "--threads", "2"])
    y = np.load(os.path.join(d, "app.npz"))
    lg = y["D0__lg"][:, 0].astype(np.float64)
    assert len(lg) == len(s_m) and float(np.max(np.abs(lg - L4.to_lg16(s_m).astype(np.float64)))) == 0.0, "dn6_apply ≠ score z maską"
    print("selftest dn6_apply: logit ≡ score z maską — OK", flush=True)
    print("SELFTEST DN6 OK", flush=True)
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        d = sys.argv[sys.argv.index("--dir") + 1]
        sys.exit(selftest(d))
    sys.exit(T5.main())
