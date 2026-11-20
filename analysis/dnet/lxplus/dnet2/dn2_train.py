#!/usr/bin/env python3
"""DNET etap 2 — trening modeli na zbiorze torów (dn2_models: ds = DeepSets, tr = transformer) dla D0, na wierszach etapu 1 z wyciągniętym
otoczeniem (dn2_extract). Ta sama pula, foldy (hash event_uid), wagi i strata (BCE + λ·DisCo) co w etapie 1 (dn_train.pool / dn_models.fit).

  dn2_train.py train --data DIR1/data --env DIR2/parts --out DIR --block NAZWA [--device auto]

Tokeny: córki (K, π; log pT, Δη, Δφ względem D, ładunek względem K, dE/dx, npts, d0, z0sinθ) + do NT torów otoczenia (te same pola bez d0/z0 —
drzewo Tracks ich nie ma) w ΔR < 0,4 (e04) lub 1,0 (e10); „dau” = same córki. Cechy globalne: 18 cech H (wierzchołek — nie do odtworzenia z torów).
Score zapisywany dla WSZYSTKICH wierszy tabel etapu 1 (NaN poza podzbiorem z otoczeniem) → analyze.py porównuje na wierszach wspólnych."""
import argparse
import glob
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.environ.get("DN1") or os.path.join(os.path.dirname(HERE), "dnet1"))
import dn_feat as DF  # noqa: E402
import dn_train as DT  # noqa: E402

LAM = DF.LAM[0]


def _c(arch, env, bkg="mc", dec=f"d{LAM}", seed=1, full=False):
    return {"sp": "D0", "model": arch, "feats": f"H+{env}", "env": env, "bkg": bkg, "dec": dec, "seed": seed, "full": full,
            "tag": f"D0.{arch}.H+{env}.{bkg}.{dec}.s{seed}"}


BLOCKS = {
    "dau": [_c("ds", "dau"), _c("tr", "dau", full=True)],                      # tylko córki (pełne parametry) + H
    "dauS": [_c("tr", "dau", seed=2), _c("tr", "dau", seed=3)],
    "ds04": [_c("ds", "e04")],                                                 # córki + otoczenie ΔR < 0,4 (wyuczona izolacja)
    "tr04": [_c("tr", "e04", full=True)],
    "tr04S2": [_c("tr", "e04", seed=2)],
    "tr04S3": [_c("tr", "e04", seed=3)],
    "tr10": [_c("tr", "e10")],                                                 # otoczenie ΔR < 1,0
    "tr04x": [_c("tr", "e04", dec="n")],                                       # bez dekorelacji (kontrola)
    "tr04sb": [_c("tr", "e04", bkg="sbmc")],                                   # tło z SB danych, dekorelacja na tle MC
}
KEYS = ["X", "Dd", "m", "pt", "eta", "phi", "uid", "truth", "ok", "swap", "m_swap", "w"]


def load(ddir, edir, sp="D0"):
    """Tabele etapu 1 ograniczone do wierszy z otoczeniem + E, nE; zwraca (S2, długości pełnych tabel)."""
    S, full = {}, {}
    hc = [DF.ALLF[sp].index(f) for f in DF.H]
    for s in DT.SAMPLES:
        p = os.path.join(ddir, f"dn_{s}.npz")
        parts = sorted(glob.glob(os.path.join(edir, f"dn2_{s}_*.npz")))
        if not os.path.isfile(p) or not parts:
            continue
        z = np.load(p, allow_pickle=False)
        T = {}
        for pre in ("a_", "p_"):
            k = f"{sp}__{pre}"
            if k + "m" not in z.files:
                continue
            ix, E, nE = [], [], []
            for f in parts:
                e = np.load(f, allow_pickle=False)
                if k + "idx" in e.files:
                    ix.append(e[k + "idx"]); E.append(e[k + "E"]); nE.append(e[k + "nE"])
            if not ix:
                continue
            ix, E, nE = np.concatenate(ix), np.concatenate(E), np.concatenate(nE)
            o = np.argsort(ix)
            ix, E, nE = ix[o], E[o], nE[o]
            full[(s, pre)] = len(z[k + "m"])
            for c in KEYS:
                if k + c in z.files:
                    v = z[k + c][ix]
                    T[pre + c] = v[:, hc] if c == "X" else v
            T[pre + "idx"], T[pre + "E"], T[pre + "nE"] = ix, E, nE
        S[s] = T
    return S, full


def tokens(Dd, E, eta, phi, env):
    """→ tok [n, nd + NT', 11] float32, mask [n, nd + NT'] float32. Dd: (pt, eta, phi, d0, z0, dedx, npts, charge); E: (pt, dη, dφ, q, dedx, npts)."""
    n, nd = Dd.shape[:2]
    Dd = np.nan_to_num(Dd.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    qK = np.where(Dd[:, 0, 7] != 0, Dd[:, 0, 7], 1.0)
    d = np.zeros((n, nd, 11), np.float32)
    d[:, :, 0] = np.log(np.clip(Dd[:, :, 0], 0.05, None))
    d[:, :, 1] = Dd[:, :, 1] - eta[:, None]
    d[:, :, 2] = (Dd[:, :, 2] - phi[:, None] + np.pi) % (2 * np.pi) - np.pi
    d[:, :, 3] = Dd[:, :, 7] * qK[:, None]
    d[:, :, 4] = np.clip(Dd[:, :, 5], 0, 10); d[:, :, 5] = Dd[:, :, 6]; d[:, :, 6] = Dd[:, :, 3]; d[:, :, 7] = Dd[:, :, 4]
    d[:, 0, 8] = 1.0; d[:, 1:, 9] = 1.0
    md = np.ones((n, nd), np.float32)
    if env == "dau":
        return d, md
    E = np.nan_to_num(E.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    rmax = {"e04": 0.4, "e10": 1.0}[env]
    me = ((E[:, :, 0] > 0) & (np.hypot(E[:, :, 1], E[:, :, 2]) < rmax)).astype(np.float32)
    e = np.zeros((n, E.shape[1], 11), np.float32)
    e[:, :, 0] = np.log(np.clip(E[:, :, 0], 0.05, None)); e[:, :, 1] = E[:, :, 1]; e[:, :, 2] = E[:, :, 2]
    e[:, :, 3] = E[:, :, 3] * qK[:, None]; e[:, :, 4] = np.clip(E[:, :, 4], 0, 10); e[:, :, 5] = E[:, :, 5]; e[:, :, 10] = 1.0
    e *= me[:, :, None]
    if env == "e04":                                   # w ΔR < 0,4 rzadko > 12 torów — krótsza sekwencja = szybszy trening
        e, me = e[:, :12], me[:, :12]
    return np.concatenate([d, e], 1), np.concatenate([md, me], 1)


def run_config(c, S, full, a, log):
    import torch
    import dn_models as DM
    import dn2_models as D2
    K = a.folds
    P = DT.pool(S, "D0", c["bkg"], a.cap_sig, a.cap_bkg)
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in P["src"]])  # noqa: E731
    G = g("X")
    tok, mask = tokens(g("Dd"), g("E"), g("eta"), g("phi"), c["env"])
    fold = DF.fold_of(P["uid"], K)
    val = DF.unit_hash(P["uid"], salt=23) >= 0.9
    dec, lam = c["dec"][0], float(c["dec"][1:] or 0)
    t0 = time.time()
    log(f"[{c['tag']}] pula {len(P['y'])} (S {int(P['y'].sum())}, B {int((P['role'] == 'B').sum())}, D {int((P['role'] == 'D').sum())}), tokenów {tok.shape[1]}, "
        f"średnio ważnych {mask.sum(1).mean():.1f}, urządzenie {a.device}")
    models, info = [], []
    for k in list(range(K)) + (["full"] if c["full"] else []):
        trn = (fold != k) if k != "full" else np.ones(len(fold), bool)
        u = trn & (P["w"] > 0)
        med, sc = DM.prep_stats(G[u])
        mv = mask[u] > 0
        tv = tok[u][:, :, :D2.NCONT][mv]
        tmean, tstd = tv.mean(0), tv.std(0)
        net = D2.SetNet(G.shape[1], c["model"], tmean, np.where(tstd > 0, tstd, 1.0), med, sc)
        net, h = DM.fit(net, (tok[trn], mask[trn], G[trn]), P["y"][trn], P["w"][trn], P["m"][trn], P["dmask"][trn], P["dw"][trn], val[trn], dec, lam,
                        1000 * c["seed"] + (K if k == "full" else k), log, device=a.device)
        models.append(net); info.append(h)
    out = {}
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "uid" not in T:
                continue
            f = DF.fold_of(T[pre + "uid"], K)
            sc_ = np.full(full[(s, pre)], np.nan, np.float32)
            sub = np.full(len(f), np.nan, np.float32)
            for a0 in range(0, len(f), 400_000):
                sl = slice(a0, a0 + 400_000)
                tk, mk = tokens(T[pre + "Dd"][sl], T[pre + "E"][sl], T[pre + "eta"][sl], T[pre + "phi"][sl], c["env"])
                for k in range(K):
                    w = f[sl] == k
                    if w.any():
                        sub[sl][w] = DM.predict(models[k], (tk[w], mk[w], T[pre + "X"][sl][w]), device=a.device)
            sc_[T[pre + "idx"]] = sub
            out[f"{s}__{pre[0]}"] = sc_
    os.makedirs(os.path.join(a.out, "sc"), exist_ok=True); os.makedirs(os.path.join(a.out, "meta"), exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "sc", c["tag"] + ".npz"), **out)
    meta = {**c, "features": DF.H, "token": D2.TOK, "n_feat": f"{len(DF.H)} + {tok.shape[1]}×{len(D2.TOK)}", "n_pool": int(len(P["y"])), "n_sig": int(P["y"].sum()),
            "n_bkg": int((P["role"] == "B").sum()), "nparams": DM.nparams(models[0]), "npar_unit": "wagi", "folds": K, "fit": info,
            "sec": round(time.time() - t0, 1), "threads": a.threads, "device": a.device, "hp": DM.HP, "arch": c["model"]}
    if c["full"]:
        os.makedirs(os.path.join(a.out, "onnx"), exist_ok=True)
        meta["onnx"] = DM.export_onnx(models[K], (tok[:200], mask[:200], G[:200]), os.path.join(a.out, "onnx", c["tag"] + ".onnx"), ["tokens", "mask", "globals"])
        meta["onnx"]["inputs"] = {"tokens": D2.TOK, "globals": DF.H}
        log(f"  ONNX: {meta['onnx']}")
    json.dump(meta, open(os.path.join(a.out, "meta", c["tag"] + ".json"), "w"), indent=1)
    log(f"[{c['tag']}] gotowe, {time.time() - t0:.0f}s, parametrów {meta['nparams']}")
    del torch


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train", "list"])
    ap.add_argument("--data"); ap.add_argument("--env"); ap.add_argument("--out"); ap.add_argument("--block"); ap.add_argument("--only", default="")
    ap.add_argument("--folds", type=int, default=int(os.environ.get("DN_FOLDS", "5")))
    ap.add_argument("--cap-sig", type=int, default=int(os.environ.get("DN2_CAP_SIG", "300000")))
    ap.add_argument("--cap-bkg", type=int, default=int(os.environ.get("DN2_CAP_BKG", "500000")))
    ap.add_argument("--epochs", type=int, default=int(os.environ.get("DN2_EPOCHS", "20")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("DN_CPUS", "8")))
    ap.add_argument("--device", default="auto")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for b, cs in BLOCKS.items():
            print(b, len(cs), " ".join(c["tag"] for c in cs))
        return 0
    import torch
    import dn_models as DM
    torch.set_num_threads(max(1, a.threads))
    DM.HP["epochs"] = a.epochs
    if a.device == "auto":
        a.device = "cuda" if torch.cuda.is_available() else "cpu"
    if a.device == "cuda":  # 18.11: test dymny — karta bez jądra CUDA dla tej wersji torch (np. V100) → CPU zamiast padać
        try:
            _t = torch.randn(64, 8, device="cuda")
            float((torch.cat([_t, _t], -1) @ torch.randn(16, 4, device="cuda")).sum())
        except Exception as e:  # noqa: BLE001
            print(f"UWAGA: CUDA nie działa na {torch.cuda.get_device_name(0)} ({type(e).__name__}) → CPU", flush=True)
            a.device = "cpu"
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    log(f"torch {torch.__version__}, urządzenie {a.device}" + (f" ({torch.cuda.get_device_name(0)})" if a.device == "cuda" else ""))
    S, full = load(a.data, a.env)
    log("wczytano: " + ", ".join(f"{s} a={len(T.get('a_m', []))}/{full.get((s, 'a_'))} p={len(T.get('p_m', []))}/{full.get((s, 'p_'))}" for s, T in S.items()))
    for c in BLOCKS[a.block]:
        if a.only and c["tag"] not in a.only.split(","):
            continue
        if os.path.isfile(os.path.join(os.environ.get("DN_DONE_DIR") or a.out, "meta", c["tag"] + ".json")) and os.environ.get("DN_REDO") != "1":
            log(f"[{c['tag']}] już jest — pomijam")
            continue
        run_config(c, S, full, a, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
