#!/usr/bin/env python3
"""DNET etap 4 — trening z ZAPISEM modeli foldów (potrzebne do zastosowania na pełnej statystyce danych: metoda R, bin 2–3 GeV).
Ta sama pula, foldy (hash event_uid), wagi, model (dn2_models.SetNet „tr”: 18 cech H + córki + tory otoczenia ΔR < 0,4) i strata
(BCE + λ·DisCo masy, λ = DN_LAM[0]) co w etapach 2–3; kod etapów 1–3 bez zmian.

  dn4_train.py train --data DIR1/data --env DIR2/parts --env-dp DIR4/parts_dp --out DIR --block NAZWA
  dn4_train.py list
  dn4_train.py selftest [--need-torch] [--dir KATALOG]    dane syntetyczne (D0 i D+): pule, tokeny, drzewo z zapisem; z torch: wszystkie bloki
                                                           1 epoka + odczyt modeli ≡ score OOF + dn4_apply na tabeli syntetycznej

Bloki (tag = <gatunek>.<model>.<cechy>.<tło>.<dek>.s<ziarno>):
  xgbD0, xgbDp   odniesienie: drzewo XGB H, tło z dalekich pasm danych 428770 (= D0/Dplus.xgb.H.sb.n.s1 etapu 1, ta sama pula i ziarno), z zapisem
  d0s1..d0s3     D0.tr.H+e04.mc.d5.s{1,2,3} (powtórka etapu 2 z zapisem modeli)
  dps1..dps3     Dplus.tr.H+e04.mc.d5.s{1,2,3} (trzy córki + otoczenie; tory otoczenia D+ z węzłów extdp)
  d0sb           D0.tr.H+e04.sb2mc.d5.s1: BCE na dalekich pasmach danych OBU runów (po połowie wagi), DisCo na tle MC w całym zakresie mas
Wyjście: sc/<tag>.npz (score OOF wierszy tabel etapu 1; NaN poza podzbiorem z torami), meta/<tag>.json, models/<tag>/{fold<k>.pt|json, meta.json}."""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn4_lib as L4  # noqa: E402  (ustawia sys.path na ../dnet1, ../dnet2)
import dn_feat as DF  # noqa: E402
import dn_train as DT  # noqa: E402
import dn2_train as D2T  # noqa: E402

LAM = DF.LAM[0]
SB_RUNS = ("data_428770", "data_427929")


def _tr(sp, bkg="mc", seed=1):
    return {"sp": sp, "model": "tr", "feats": f"H+{L4.ENV}", "env": L4.ENV, "bkg": bkg, "dec": f"d{LAM}", "seed": seed, "full": False,
            "lam": float(LAM), "tag": f"{sp}.tr.H+{L4.ENV}.{bkg}.d{LAM}.s{seed}"}


def _xgb(sp):
    return {"sp": sp, "model": "xgb", "feats": "H", "bkg": "sb", "dec": "n", "seed": 1, "full": False, "tag": f"{sp}.xgb.H.sb.n.s1"}


BLOCKS = {
    "xgbD0": [_xgb("D0")], "xgbDp": [_xgb("Dplus")],
    "d0s1": [_tr("D0", seed=1)], "d0s2": [_tr("D0", seed=2)], "d0s3": [_tr("D0", seed=3)],
    "dps1": [_tr("Dplus", seed=1)], "dps2": [_tr("Dplus", seed=2)], "dps3": [_tr("Dplus", seed=3)],
    "d0sb": [_tr("D0", bkg="sb2mc")],
}


# ----------------------------------------------------------------------------------------------------------- pula
def sb_mask(S, sp, s):
    """Dalekie pasma boczne danych próbki s (jak dn_train.masks dla 428770)."""
    m0 = DF.cfg["species"][sp]["m0"]
    da = S[s]
    b = da["a_ok"] & (np.abs(da["a_m"] - m0) > DF.FAR[sp])
    if sp == "D0":
        b &= np.abs(da["a_m_swap"] - m0) > DF.REFL_VETO
    return b


def pool4(S, sp, bkg, cap_sig, cap_bkg):
    """dn_train.pool + tło „sb2mc”: BCE na dalekich pasmach danych obu runów (każdy run przeważony do widma (pT, |η|) sygnału i znormalizowany
    do połowy sumy wag sygnału), człon dekorelacji na tle MC w całym zakresie mas (rola D, waga BCE = 0). Brak przecieku do oceny:
    score wiersza danych = model foldu, który jego zdarzenia nie widział (dn4_lib.score_rows)."""
    if bkg != "sb2mc":
        return DT.pool(S, sp, bkg, cap_sig, cap_bkg)
    for s in SB_RUNS:
        if s not in S:
            raise SystemExit(f"brak {s} dla tła sb2mc")
    sig, bmc, _ = DT.masks(S, sp)

    def take(s, mk, cap, salt):
        i = np.flatnonzero(mk)
        if cap and len(i) > cap:
            i = i[DF.unit_hash(S[s]["a_uid"][i], salt) < cap / len(i)]
        return i
    src = [("fd_nd", take("fd_nd", sig, cap_sig, 11), "S")]
    for j, s in enumerate(SB_RUNS):
        src.append((s, take(s, sb_mask(S, sp, s), cap_bkg // 2, 13 + j), "B"))
    src.append(("mb_nd", take("mb_nd", bmc, cap_bkg // 2, 12), "D"))
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in src])  # noqa: E731
    role = np.concatenate([np.full(len(i), r) for _, i, r in src])
    part = np.concatenate([np.full(len(i), j) for j, (_, i, _) in enumerate(src)])
    P = {"src": src, "role": role, "y": role == "S", "m": g("m").astype(np.float32), "uid": g("uid"), "pt": g("pt"), "eta": g("eta")}
    w0 = np.where(role == "S", g("w"), 1.0).astype(np.float64)
    w = np.zeros(len(role))
    s_ = role == "S"
    w[s_] = w0[s_]
    nB = sum(1 for _, _, r in src if r == "B")
    for j, (_, _, r) in enumerate(src):
        if r == "S":
            continue
        b_ = part == j
        if b_.any():
            w[b_] = DT.reweight(P["pt"][b_], P["eta"][b_], w0[b_], P["pt"][s_], P["eta"][s_], w0[s_]) / (nB if r == "B" else 1)
    w /= w[s_ | (role == "B")].mean()
    P["dmask"] = role == "D"
    P["dw"] = np.where(P["dmask"], w, 0.0).astype(np.float32)
    P["w"] = np.where(role == "D", 0.0, w).astype(np.float32)
    return P


# ----------------------------------------------------------------------------------------------------------- trening
def _finish(c, a, out, meta, log, t0):
    os.makedirs(os.path.join(a.out, "sc"), exist_ok=True); os.makedirs(os.path.join(a.out, "meta"), exist_ok=True)
    np.savez_compressed(os.path.join(a.out, "sc", c["tag"] + ".npz"), **out)
    meta["sec"] = round(time.time() - t0, 1)
    json.dump(meta, open(os.path.join(a.out, "models", c["tag"], "meta.json"), "w"), indent=1)
    json.dump(meta, open(os.path.join(a.out, "meta", c["tag"] + ".json"), "w"), indent=1)        # meta/ = znacznik „gotowe” (kopiowany na końcu)
    log(f"[{c['tag']}] gotowe, {time.time() - t0:.0f}s, parametrów {meta['nparams']}")
    return meta


def run_xgb(c, a, log):
    """Drzewo odniesienia na PEŁNYCH tabelach etapu 1 (jak dn_train.run_config) + zapis modeli foldów."""
    sp, K = c["sp"], a.folds
    S = DT.load_species(a.data, sp)
    cols = [DF.ALLF[sp].index(f) for f in DF.H]
    P = DT.pool(S, sp, c["bkg"], a.cap_sig_xgb, a.cap_bkg_xgb)
    X = DT.pool_x(S, P, "X", cols)
    fold = DF.fold_of(P["uid"], K)
    val = DF.unit_hash(P["uid"], salt=23) >= 0.9
    t0 = time.time()
    log(f"[{c['tag']}] pula {len(P['y'])} (S {int(P['y'].sum())}, B {int((P['role'] == 'B').sum())}), {len(cols)} cech")
    mdir = os.path.join(a.out, "models", c["tag"])
    os.makedirs(mdir, exist_ok=True)
    models, info = [], []
    for k in range(K):
        ub = (fold != k) & (P["w"] > 0)
        tm = time.time()
        m = DT.train_xgb(X[ub], P["y"][ub], P["w"][ub], val[ub], 1000 * c["seed"] + k, a.threads)
        m.save_model(os.path.join(mdir, f"fold{k}.json"))
        models.append(m); info.append({"best_it": int(m.best_iteration), "sec": round(time.time() - tm, 1)})
        log(f"  fold {k}: best_it {m.best_iteration}, {time.time() - tm:.0f}s")
    out = {}
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "uid" in T:
                out[f"{s}__{pre[0]}"] = L4.score_rows("xgb", models, T[pre + "uid"], T[pre + "X"][:, cols])
    meta = {**c, "features": DF.H, "n_feat": len(DF.H), "n_pool": int(len(P["y"])), "n_sig": int(P["y"].sum()), "n_bkg": int((P["role"] == "B").sum()),
            "nparams": int(len(models[0].get_booster().get_dump())), "npar_unit": "drzewa (fold 0)", "folds": K, "fit": info, "threads": a.threads,
            "device": "cpu", "hp": {k: v for k, v in DF.cfg["training"]["xgb"].items()}}
    return _finish(c, a, out, meta, log, t0), models, S


def run_tr(c, S, full, a, log):
    import dn_models as DM
    import dn2_models as D2
    sp, K = c["sp"], a.folds
    P = pool4(S, sp, c["bkg"], a.cap_sig, a.cap_bkg)
    g = lambda k: np.concatenate([S[s]["a_" + k][i] for s, i, _ in P["src"]])  # noqa: E731
    G = g("X")
    tok, mask = D2T.tokens(g("Dd"), g("E"), g("eta"), g("phi"), c["env"])
    nd = int(DF.cfg["species"][sp]["n_daughters"])
    fold = DF.fold_of(P["uid"], K)
    val = DF.unit_hash(P["uid"], salt=23) >= 0.9
    dec, lam = c["dec"][0], float(c["lam"])
    t0 = time.time()
    log(f"[{c['tag']}] pula {len(P['y'])} (S {int(P['y'].sum())}, B {int((P['role'] == 'B').sum())}, D {int((P['role'] == 'D').sum())}; "
        + ", ".join(f"{s}:{r}={len(i)}" for s, i, r in P["src"]) + f"), tokenów {tok.shape[1]} (córek {nd}), średnio ważnych {mask.sum(1).mean():.2f}, "
        f"urządzenie {a.device}")
    mdir = os.path.join(a.out, "models", c["tag"])
    models, info = [], []
    for k in range(K):
        trn = fold != k
        u = trn & (P["w"] > 0)
        med, sc = DM.prep_stats(G[u])
        mv = mask[u] > 0
        tv = tok[u][:, :, :D2.NCONT][mv]
        tmean, tstd = tv.mean(0), tv.std(0)
        net = D2.SetNet(G.shape[1], c["model"], tmean, np.where(tstd > 0, tstd, 1.0), med, sc)
        net, h = DM.fit(net, (tok[trn], mask[trn], G[trn]), P["y"][trn], P["w"][trn], P["m"][trn], P["dmask"][trn], P["dw"][trn], val[trn], dec, lam,
                        1000 * c["seed"] + k, log, device=a.device)
        L4.save_tr(mdir, k, net, G.shape[1], c["model"])
        models.append(net); info.append(h)
    del tok, mask
    out = {}
    for s, T in S.items():
        for pre in ("a_", "p_"):
            if pre + "uid" not in T:
                continue
            sc_ = np.full(full[(s, pre)], np.nan, np.float32)
            sc_[T[pre + "idx"]] = L4.score_rows("tr", models, T[pre + "uid"], T[pre + "X"], T[pre + "Dd"], T[pre + "E"], T[pre + "eta"], T[pre + "phi"],
                                                device=a.device)
            out[f"{s}__{pre[0]}"] = sc_
    meta = {**c, "features": DF.H, "token": D2.TOK, "n_feat": f"{len(DF.H)} + {nd + 12}×{len(D2.TOK)}", "n_pool": int(len(P["y"])), "n_sig": int(P["y"].sum()),
            "n_bkg": int((P["role"] == "B").sum()), "n_src": {f"{s}:{r}": int(len(i)) for s, i, r in P["src"]}, "nparams": DM.nparams(models[0]),
            "npar_unit": "wagi", "folds": K, "fit": info, "threads": a.threads, "device": a.device, "hp": dict(DM.HP), "arch": c["model"]}
    return _finish(c, a, out, meta, log, t0), models


# ----------------------------------------------------------------------------------------------------------- test mechaniki
NT = 20


def _synth_table(rng, n, kind, sp):
    scfg = DF.cfg["species"][sp]
    nd, m0, (lo, hi) = int(scfg["n_daughters"]), scfg["m0"], scfg["mass_range"]
    nf = len(DF.ALLF[sp])
    truth = {"fd": np.where(rng.random(n) < 0.85, 0, 8), "mb": np.where(rng.random(n) < 0.1, 0, 8), "data": np.full(n, -1)}[kind].astype(np.int16)
    sig = (truth >= 0) & (truth < 8)
    if kind == "data":
        sig = rng.random(n) < 0.1
    X = rng.normal(size=(n, nf)).astype(np.float32)
    X[:, :3] += sig[:, None] * 1.0
    pt = np.clip(3.0 + rng.exponential(3.0, n), 3.0, 29.0)
    eta, phi = rng.uniform(-2.4, 2.4, n), rng.uniform(-np.pi, np.pi, n)
    m = np.where(sig, m0 + rng.normal(0, 0.02, n), rng.uniform(lo, hi, n)).clip(lo, hi - 1e-4)
    Dd = np.zeros((n, nd, 8), np.float32)
    fr = rng.dirichlet(np.ones(nd) * 3, n)
    Dd[:, :, 0] = pt[:, None] * fr
    Dd[:, :, 1] = eta[:, None] + rng.normal(0, 0.2, (n, nd)); Dd[:, :, 2] = phi[:, None] + rng.normal(0, 0.2, (n, nd))
    Dd[:, :, 3] = rng.normal(0, 0.1, (n, nd)); Dd[:, :, 4] = rng.normal(0, 0.2, (n, nd)); Dd[:, :, 5] = rng.uniform(0.5, 3, (n, nd))
    Dd[:, :, 6] = rng.integers(5, 20, (n, nd)); q = rng.choice([-1.0, 1.0], n); Dd[:, 0, 7] = q; Dd[:, 1:, 7] = -q[:, None]
    T = {"X": X, "Dd": Dd, "m": m.astype(np.float32), "pt": pt.astype(np.float32), "eta": eta.astype(np.float32), "phi": phi.astype(np.float32),
         "uid": rng.integers(1, 2 ** 62, n).astype(np.uint64), "truth": truth, "ok": np.ones(n, bool), "swap": np.zeros(n, np.int8),
         "m_swap": (m + 0.1).astype(np.float32), "w": np.ones(n, np.float32), "e_grl": np.ones(n, np.int32),
         "e_tp_HLT_mb_mbts_L1MBTS_1": (rng.random(n) < 0.5).astype(np.int8), "e_tp_HLT_noalg_L1J12": (rng.random(n) < 0.5).astype(np.int8)}
    ntrk = np.minimum(rng.poisson(np.where(sig, 4.0, 10.0)), NT)
    pe = (0.1 + rng.exponential(0.8, (n, NT))) * (np.arange(NT)[None, :] < ntrk[:, None])
    pe = -np.sort(-pe, axis=1)
    r, an = np.sqrt(rng.random((n, NT))), rng.uniform(0, 2 * np.pi, (n, NT))
    E = np.stack([pe, r * np.cos(an), r * np.sin(an), rng.choice([-1.0, 1.0], (n, NT)), rng.uniform(0.5, 3, (n, NT)),
                  rng.integers(3, 20, (n, NT)).astype(float)], 2)
    E *= (pe > 0)[:, :, None]
    in04 = (pe > 0) & (r < 0.4)
    T["E"], T["nE"] = E.astype(np.float16), np.stack([in04.sum(1), ntrk], 1).astype(np.int16)
    T["pc40"] = (pe * (in04 & (pe > 0.5))).sum(1).astype(np.float32)
    return T


def synth(ddir, edirs, seed=5, n=5000):
    """Tabele dn_<próbka>.npz (D0 i D+) i części dn2_<próbka>_<k>.npz (D0 → edirs['D0'], D+ → edirs['Dplus']) z liczb losowych."""
    rng = np.random.default_rng(seed)
    os.makedirs(ddir, exist_ok=True)
    for d in edirs.values():
        os.makedirs(d, exist_ok=True)
    keys = D2T.KEYS + ["e_grl", "e_tp_HLT_mb_mbts_L1MBTS_1", "e_tp_HLT_noalg_L1J12"]
    for s, kind in (("fd_nd", "fd"), ("mb_nd", "mb"), ("data_428770", "data"), ("data_427929", "data")):
        z = {}
        parts = {sp: [{}, {}] for sp in edirs}
        for sp in edirs:
            for pre, nn_, frac in (("a_", n, 0.8), ("p_", n // 4, 1.0)):
                if pre == "p_" and (kind == "mb" or sp != "D0"):
                    continue
                T = _synth_table(rng, nn_, kind, sp)
                k = f"{sp}__{pre}"
                for c in keys:
                    z[k + c] = T[c]
                idx = rng.permutation(np.flatnonzero(rng.random(nn_) < frac))
                for j, ii in enumerate((idx[: len(idx) // 2], idx[len(idx) // 2:])):
                    parts[sp][j].update({k + "idx": ii, k + "E": T["E"][ii], k + "nE": T["nE"][ii], k + "pc40": T["pc40"][ii]})
        np.savez(os.path.join(ddir, f"dn_{s}.npz"), **z)
        for sp, d in edirs.items():
            for j, p in enumerate(parts[sp]):
                np.savez(os.path.join(d, f"dn2_{s}_{j}.npz"), **p)


def selftest(a):
    import tempfile
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    ctx = None
    if a.dir:
        d = a.dir
        os.makedirs(d, exist_ok=True)
    else:
        ctx = tempfile.TemporaryDirectory(prefix="dn4_selftest_")
        d = ctx.name
    try:
        ddir, edirs, odir = os.path.join(d, "data"), {"D0": os.path.join(d, "parts"), "Dplus": os.path.join(d, "parts_dp")}, os.path.join(d, "out")
        synth(ddir, edirs)
        DF.cfg["training"]["xgb"].update(n_estimators=30, learning_rate=0.3)
        ns = argparse.Namespace(folds=2, cap_sig=0, cap_bkg=0, cap_sig_xgb=0, cap_bkg_xgb=0, device="cpu", threads=2, out=odir, data=ddir)
        tags = [c["tag"] for cs in BLOCKS.values() for c in cs]
        assert len(set(tags)) == len(tags)
        SS = {}
        for sp in ("D0", "Dplus"):
            S, full = D2T.load(ddir, edirs[sp], sp)
            SS[sp] = (S, full)
            log(f"{sp}: " + ", ".join(f"{s} a={len(T.get('a_m', []))}/{full.get((s, 'a_'))}" for s, T in S.items()))
            # pule: mc i sb2mc
            P = pool4(S, sp, "sb2mc", 0, 0)
            B, Dm, Sg = P["role"] == "B", P["role"] == "D", P["y"]
            assert B.any() and Dm.any() and (P["w"][Dm] == 0).all() and (P["dw"][~Dm] == 0).all() and (P["dw"][Dm] > 0).any(), "sb2mc: role/wagi"
            srcB = [(s, i) for s, i, r in P["src"] if r == "B"]
            assert [s for s, _ in srcB] == list(SB_RUNS)
            n0 = len(P["src"][0][1]); n1 = len(srcB[0][1])
            w1, w2 = P["w"][n0:n0 + n1].sum(), P["w"][n0 + n1:n0 + n1 + len(srcB[1][1])].sum()
            assert abs(w1 / w2 - 1) < 1e-3 and abs((w1 + w2) / P["w"][Sg].sum() - 1) < 1e-3, f"sb2mc: wagi runów {w1} {w2} {P['w'][Sg].sum()}"
            for s, i in srcB:                                        # tło = tylko dalekie pasma (poza oknem i pasmem pośrednim)
                assert (np.abs(S[s]["a_m"][i] - DF.cfg["species"][sp]["m0"]) > DF.FAR[sp]).all()
            nd = int(DF.cfg["species"][sp]["n_daughters"])
            T = S["mb_nd"]
            tok, mask = D2T.tokens(T["a_Dd"], T["a_E"], T["a_eta"], T["a_phi"], L4.ENV)
            assert tok.shape[1:] == (nd + 12, 11) and mask[:, :nd].all() and (tok[:, 0, 8] == 1).all() and (tok[:, 1:nd, 9] == 1).all() \
                and (tok[:, nd:, 10][mask[:, nd:] > 0] == 1).all(), f"{sp}: tokeny"
        log("numpy OK: pule (mc, sb2mc: wagi runów po połowie, tylko dalekie pasma), tokeny D0 (2 córki) i D+ (3 córki)")
        # drzewa odniesienia z zapisem → odczyt ≡ score OOF
        for b in ("xgbD0", "xgbDp"):
            c = BLOCKS[b][0]
            meta, models, S1 = run_xgb(c, ns, log)
            _, ms = L4.load_models(os.path.join(odir, "models"), c["tag"])
            z = np.load(os.path.join(odir, "sc", c["tag"] + ".npz"))
            cols = [DF.ALLF[c["sp"]].index(f) for f in DF.H]
            for s, T in S1.items():
                v = L4.score_rows("xgb", ms, T["a_uid"], T["a_X"][:, cols])
                assert np.max(np.abs(v - z[f"{s}__a"])) < 1e-6 and np.isfinite(v).all(), f"{c['tag']} {s}: model z pliku ≠ score OOF"
        log("xgb OK: modele foldów zapisane, odczyt ≡ score OOF")
        try:
            import torch
            import dn_models as DM
        except ImportError:
            if a.need_torch:
                raise
            log("torch niedostępny — sieci (trening, zapis/odczyt, dn4_apply dla sieci) NIE przetestowane")
            _apply_test(d, ddir, edirs, odir, log, with_tr=False)
            return 0
        torch.set_num_threads(2)
        hp0 = dict(DM.HP)
        DM.HP.update(epochs=1, batch=256)
        try:
            for b, cs in BLOCKS.items():
                for c in cs:
                    if c["model"] == "xgb":
                        continue
                    S, full = SS[c["sp"]]
                    meta, models = run_tr(c, S, full, ns, log)
                    _, ms = L4.load_models(os.path.join(odir, "models"), c["tag"])
                    z = np.load(os.path.join(odir, "sc", c["tag"] + ".npz"))
                    for s, T in S.items():
                        for pre in ("a_", "p_"):
                            if pre + "idx" not in T:
                                continue
                            v = z[f"{s}__{pre[0]}"]
                            fin = np.isfinite(v)
                            assert len(v) == full[(s, pre)] and fin.sum() == len(T[pre + "idx"]) and fin[T[pre + "idx"]].all(), f"{c['tag']} {s} {pre}: wzór NaN"
                            r = L4.score_rows("tr", ms, T[pre + "uid"], T[pre + "X"], T[pre + "Dd"], T[pre + "E"], T[pre + "eta"], T[pre + "phi"])
                            assert np.max(np.abs(r - v[T[pre + "idx"]])) < 1e-5, f"{c['tag']} {s} {pre}: model z pliku ≠ score OOF"
                    assert len(meta["fit"]) == 2
            # szybkość wnioskowania (do planowania węzłów app)
            S, full = SS["Dplus"]
            T = S["mb_nd"]
            _, ms = L4.load_models(os.path.join(odir, "models"), BLOCKS["dps1"][0]["tag"])
            rep = 20
            big = {k: np.concatenate([T["a_" + k]] * rep) for k in ("uid", "X", "Dd", "E", "eta", "phi")}
            tt = time.time()
            L4.score_rows("tr", ms, big["uid"], big["X"], big["Dd"], big["E"], big["eta"], big["phi"])
            log(f"wnioskowanie tr (D+, 15 tokenów, 2 wątki): {len(big['uid']) / (time.time() - tt):.0f} wierszy/s")
        finally:
            DM.HP.clear(); DM.HP.update(hp0)
        log("torch OK: wszystkie bloki (1 epoka, 2 foldy, CPU), zapis/odczyt modeli ≡ score OOF")
        _apply_test(d, ddir, edirs, odir, log, with_tr=True)
    finally:
        if ctx is not None:
            ctx.cleanup()
    return 0


def _apply_test(d, ddir, edirs, odir, log, with_tr):
    """dn4_apply na tabeli syntetycznej jednej próbki (D0 + D+ w jednym katalogu części) ≡ score z treningu."""
    import shutil
    import dn4_apply as AP
    ed = os.path.join(d, "parts_both")
    os.makedirs(ed, exist_ok=True)
    s = "data_428770"
    for j in range(2):                                             # dn2_extract --species D0,Dplus zapisuje oba gatunki w jednym pliku
        z = {}
        for sp, e in edirs.items():
            z.update(dict(np.load(os.path.join(e, f"dn2_{s}_{j}.npz"))))
        np.savez(os.path.join(ed, f"dn2_{s}_{j}.npz"), **z)
    one = os.path.join(d, "data_one")
    os.makedirs(one, exist_ok=True)
    shutil.copy(os.path.join(ddir, f"dn_{s}.npz"), one)
    out = os.path.join(d, "app.npz")
    AP.main(["--sample", s, "--data", one, "--env", ed, "--models", os.path.join(odir, "models"), "--out", out])
    z = np.load(out, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    n = 0
    for sp in ("D0", "Dplus"):
        tags = meta["tags"][sp]
        assert tags and (with_tr is False or any(".tr." in t for t in tags)), f"{sp}: brak modeli w app"
        full = np.load(os.path.join(one, f"dn_{s}.npz"))
        uid = z[f"{sp}__uid"]
        o = np.argsort(full[f"{sp}__a_uid"]); pos = o[np.searchsorted(full[f"{sp}__a_uid"][o], uid)]
        assert np.array_equal(full[f"{sp}__a_uid"][pos], uid) and np.allclose(full[f"{sp}__a_m"][pos], z[f"{sp}__m"])
        for j, t in enumerate(tags):
            ref = np.load(os.path.join(odir, "sc", t + ".npz"))[f"{s}__a"][pos]
            assert np.array_equal(L4.to_lg16(ref), z[f"{sp}__lg"][:, j]), f"{t}: app ≠ score treningu"
            n += 1
        assert z[f"{sp}__cls"].dtype == np.uint8 and (z[f"{sp}__cls"] & 1).all()
    log(f"apply OK: {n} modeli, logit float16 ≡ score treningu, klasy trygera zapisane")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train", "list", "selftest"])
    ap.add_argument("--data"); ap.add_argument("--env"); ap.add_argument("--env-dp"); ap.add_argument("--out"); ap.add_argument("--block")
    ap.add_argument("--folds", type=int, default=int(os.environ.get("DN_FOLDS", "5")))
    ap.add_argument("--cap-sig", type=int, default=int(os.environ.get("DN4_CAP_SIG", "300000")))
    ap.add_argument("--cap-bkg", type=int, default=int(os.environ.get("DN4_CAP_BKG", "500000")))
    ap.add_argument("--cap-sig-xgb", type=int, default=int(os.environ.get("DN_CAP_SIG", "400000")))
    ap.add_argument("--cap-bkg-xgb", type=int, default=int(os.environ.get("DN_CAP_BKG", "800000")))
    ap.add_argument("--epochs", type=int, default=int(os.environ.get("DN4_EPOCHS", "20")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("DN_CPUS", "8")))
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--need-torch", action="store_true"); ap.add_argument("--dir", default="")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for b, cs in BLOCKS.items():
            print(b, len(cs), " ".join(c["tag"] for c in cs))
        return 0
    if a.cmd == "selftest":
        return selftest(a)
    if a.block not in BLOCKS:
        raise SystemExit(f"nieznany blok {a.block}; są: {' '.join(BLOCKS)}")
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    cache = {}
    for c in BLOCKS[a.block]:
        if os.path.isfile(os.path.join(os.environ.get("DN_DONE_DIR") or a.out, "meta", c["tag"] + ".json")) and os.environ.get("DN_REDO") != "1":
            log(f"[{c['tag']}] już jest — pomijam")
            continue
        if c["model"] == "xgb":
            run_xgb(c, a, log)
            continue
        import torch
        import dn_models as DM
        torch.set_num_threads(max(1, a.threads))
        DM.HP["epochs"] = a.epochs
        if c["sp"] not in cache:
            cache.clear()
            edir = a.env if c["sp"] == "D0" else a.env_dp
            cache[c["sp"]] = D2T.load(a.data, edir, c["sp"])
            S, full = cache[c["sp"]]
            log(f"torch {torch.__version__}; wczytano {c['sp']} ({edir}): "
                + ", ".join(f"{s} a={len(T.get('a_m', []))}/{full.get((s, 'a_'))} p={len(T.get('p_m', []))}/{full.get((s, 'p_'))}" for s, T in S.items()))
            if "fd_nd" not in S or "mb_nd" not in S:
                raise SystemExit(f"brak torów otoczenia {c['sp']} dla fd_nd/mb_nd w {edir}")
        run_tr(c, *cache[c["sp"]], a, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
