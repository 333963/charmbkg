#!/usr/bin/env python3
"""FMM — udział MisMatch* (f_MM) w tle pod pikiem D0: prawda MC vs estymatory na danych z par bliźniaczych K↔π (projekt B).
Otwarty punkt METHODS_NOTE §7: „f_MM in data can be measured directly from hypo_swap … not yet done”.

Wejście: fmm_<próbka>.npz z etapu ORCH `fmm` (Q_lxplus/fmm/fmm_extract.py; próbki mb_nd, fd_nd, data_428770, data_427929).
Wyjście: RAPORT.md + results.json w --out.  Użycie:
  ~/Documents/charm/.venv_h5/bin/python charm_work/FMM_data/analyze.py [--in FMM_data/runs/lxplus] [--out FMM_data] [--tag ""]

Definicje (jak B1-D0): m0 = 1,86484; SR |m−m0| < 55 MeV; SB |m−m0| > 75 MeV w oknie 1,70–2,05; weto refleksji |m_swap−m0| > 25 MeV
wszędzie; tło = truth ≥ 8; MisMatch* = truth 8–11; f_MM = N(8–11)/N(≥ 8). WP: progi score H z meta modelu (per bin pT).
Estymatory f_MM na danych (żaden nie korzysta z prawdy w danych; `hypo_swap` istnieje TYLKO w MC):
  E1  „twin-rate”: r = odsetek kandydatów, których bliźniak K↔π też leży w oknie ntupla; f = (r − r_O)/(r_MM − r_O) (kalibracja MC).
  E2  dopasowanie szablonów rozkładu dyskryminanty D (XGB MisMatch* vs reszta tła, uczony na ¼ MC) — 3 składowe: sygnał (fd),
      MisMatch*, reszta tła; warianty cech: twin (tylko obserwable bliźniacze), pid (dE/dx + pędy torów), all (wszystko + cechy H).
  Wielkość „lokalna” (ta, do której odnoszą się czułości B1 cat_local): Δf_loc = (f̂_SR − f̂_SB)_dane − (f_SR − f_SB)_MC.
Domknięcie: kalibracja (trening ćw. q0, szablony q1) vs test (q2+q3 jako pseudo-dane): nominalne, wstrzyknięcie lokalne MisMatch*×1,5·g(m),
  oraz test „stress”: pseudo-dane przeważone wagami dane/MC (klasyfikator w SB, jak SALAD) — prawdziwy skład znany → obciążenie estymatora
  wywołane samym niedomodelowaniem cech = systematyka modelu estymatora.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import xgboost as xgb

M0, SR, SBMIN, VETO = 1.86484, 0.055, 0.075, 0.025
PT = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0), "3-30": (3.0, 30.0)}
WPS = ["0.80", "0.70", "0.50"]
SENS = {"D0": {"0.80": -1.0, "0.70": -1.7, "0.50": -4.6}, "D+": {"0.80": -3.3, "0.70": -4.1, "0.50": -6.4}}   # % ε_B na +1 pp (B1 §2, §6e)
ASSUMED = {"D0": {"0.80": 5.0, "0.70": 10.0, "0.50": 15.0}, "D+": {"0.80": 5.0, "0.70": 10.0, "0.50": 25.0}}   # przyjęta syst. modelu tła [%]
MBTS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
NBOOT = 200
RNG = np.random.default_rng(20261117)


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------------------------------------------------------ wczytanie
def load(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    d = {k: z[k] for k in z.files if k != "meta"}
    n = len(d["m"])
    d["dm"] = np.abs(d["m"].astype(np.float64) - M0)
    d["veto"] = np.abs(d["m_swap"].astype(np.float64) - M0) > VETO
    d["sr"] = (d["dm"] < SR) & d["veto"]
    d["sb"] = (d["dm"] > SBMIN) & d["veto"]
    d["tw_in"] = (d["m_swap"] >= meta["window"][0]) & (d["m_swap"] <= meta["window"][1])
    thr = np.full(n, np.inf)
    d["wp"] = {}
    for wp in WPS:
        t = np.full(n, np.inf)
        for b, (lo, hi) in PT.items():
            if b in meta["thresholds"]:
                s = (d["pt"] >= lo) & (d["pt"] < hi)
                t[s] = meta["thresholds"][b][wp]
        d["wp"][wp] = d["score"] >= t
    u = d["uid"].astype(np.uint64)
    with np.errstate(over="ignore"):
        x = (u ^ (u >> np.uint64(33))) * np.uint64(0xFF51AFD7ED558CCD)
        x = x ^ (x >> np.uint64(33))
    d["q"] = (x % np.uint64(4)).astype(np.int8)          # ćwiartka po zdarzeniu (bliźniaki razem)
    d["w1"] = np.ones(n)
    return d, meta


def cats(d):
    t, hs = d["truth"], d["hypo_swap"]
    return {"sig": (t <= 1) & (t >= 0) & (hs != 1), "swap": (t <= 1) & (t >= 0) & (hs == 1), "c27": (t >= 2) & (t <= 7),
            "MM": (t >= 8) & (t <= 11), "O": t >= 12}


def feats(d, variant):
    f = lambda k: d[k].astype(np.float32)   # noqa: E731
    k0 = d["k_idx"] == 0
    pK, pPi = np.where(k0, f("dpt0"), f("dpt1")), np.where(k0, f("dpt1"), f("dpt0"))
    dK, dPi = np.where(k0, f("dedx0"), f("dedx1")), np.where(k0, f("dedx1"), f("dedx0"))
    eK, ePi = np.where(k0, f("deta0"), f("deta1")), np.where(k0, f("deta1"), f("deta0"))
    twin = [f("m_swap"), f("m_swap") - f("m"), (pK - pPi) / (pK + pPi), f("f_abs_cos_theta_star")]
    pid = [dK, dPi, np.where(k0, d["npts0"], d["npts1"]).astype(np.float32), np.where(k0, d["npts1"], d["npts0"]).astype(np.float32),
           pK * np.cosh(eK), pPi * np.cosh(ePi), np.abs(eK), np.abs(ePi)]
    if variant == "twin":
        cols = twin
    elif variant == "pid":
        cols = pid
    else:
        hk = sorted(k for k in d if k.startswith("f_"))
        cols = twin + pid + [f(k) for k in hk] + [f("score"), f("dd00"), f("dd01"), f("dz00"), f("dz01")]
    return np.column_stack(cols)


def train(X, y, w=None, n=200):
    cl = xgb.XGBClassifier(n_estimators=n, max_depth=4, learning_rate=0.1, subsample=0.8, colsample_bytree=0.8, tree_method="hist",
                           n_jobs=os.cpu_count() or 4, min_child_weight=20, verbosity=0)
    cl.fit(X, y, sample_weight=w)
    return cl


def auc(y, s, w=None):
    o = np.argsort(s, kind="mergesort")
    y, w = np.asarray(y)[o], (np.ones(len(o)) if w is None else np.asarray(w)[o])
    wn = np.cumsum(w * (1 - y))
    den = (w * y).sum() * (w * (1 - y)).sum()
    return float((w * y * wn).sum() / den) if den > 0 else float("nan")


# ------------------------------------------------------------------------------------------------------------------ dopasowanie
def hist(x, edges, w):
    i = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
    nb = len(edges) - 1
    return np.bincount(i, w, nb), np.bincount(i, w * w, nb)


def _solve(d, vd, T, VT, prior=None):
    """WLS: d ≈ Σ a_k T_k (T_k znormalizowane do 1); wariancja = dane + szablony (iteracyjnie)."""
    a = np.full(T.shape[1], d.sum() / T.shape[1])
    for _ in range(4):
        v = np.maximum(vd, 1.0) + (VT * a ** 2).sum(1)
        s = 1.0 / np.sqrt(v)
        A, y = T * s[:, None], d * s
        if prior is not None:                                         # wiersz kary: (a_S − N_S0)/σ_S (sygnał = kolumna 0)
            r = np.zeros(T.shape[1]); r[0] = 1.0 / prior[1]
            A, y = np.vstack([A, r]), np.r_[y, prior[0] / prior[1]]
        a = np.linalg.lstsq(A, y, rcond=None)[0]     # bez ograniczenia a ≥ 0: f̂ nieobciążone przy granicy (może wyjść < 0)
    chi2 = float((((d - T @ a) ** 2) / v).sum())
    return a, chi2


def tfit(d, vd, tm, with_sig=True, nboot=NBOOT, prior=None):
    """tm: {'S': (h, v), 'MM': (h, v), 'R': (h, v), 'c': udział truth<8 w R}. Zwraca f̂ (MM/(MM+O)), błąd bootstrap (dane+szablony), N_S."""
    names = (["S"] if with_sig else []) + ["MM", "R"]

    def norm(k, fluct):
        h, v = tm[k]
        h = h + (RNG.normal(0, np.sqrt(np.maximum(v, 1e-12))) if fluct else 0.0)
        h = np.maximum(h, 0)
        s = max(h.sum(), 1e-9)
        return h / s, v / s ** 2

    def one(fluct):
        cols = [norm(k, fluct) for k in names]
        T, VT = np.column_stack([c[0] for c in cols]), np.column_stack([c[1] for c in cols])
        dd = d + (RNG.normal(0, np.sqrt(np.maximum(vd, 1.0))) if fluct else 0.0)
        pr = None
        if with_sig and prior is not None:
            pr = (prior[0] + (RNG.normal(0, prior[1]) if fluct else 0.0), prior[1])
        a, chi2 = _solve(dd, vd, T, VT, pr)
        aMM, aR = a[-2], a[-1]
        den = aMM + aR * (1 - tm["c"])
        return (aMM / den if den > 0 else np.nan), (a[0] if with_sig else 0.0), chi2

    f, ns, chi2 = one(False)
    B = np.array([one(True)[:2] for _ in range(nboot)])
    return {"f": float(f), "e": float(np.nanstd(B[:, 0])), "ns": float(ns), "ens": float(np.nanstd(B[:, 1])), "chi2": chi2,
            "ndf": int(len(d) - len(names)), "n": float(d.sum())}


def sig_prior(n_sr, n_sb, T):
    """Łagodne ograniczenie sygnału w SR: N_S0 = N_SR − N_SB·(SR/SB)_MC tła; σ = 0,5·N_S0 + 1 % N_SR + stat (stosunek SR/SB tła w danych
    nie jest znany lepiej niż ~1 % → sygnał nie jest ustalony z liczenia; ograniczenie tylko stabilizuje dopasowanie)."""
    rho = (T["sr"][4] + T["sr"][5]) / max(T["sb"][4] + T["sb"][5], 1)
    ns0 = max(n_sr - n_sb * rho, 0.0)
    return ns0, 0.5 * ns0 + 0.01 * n_sr + float(np.sqrt(n_sr + rho * rho * n_sb))


def edges_for(x, w, n):
    nb = 12 if n >= 50_000 else 8 if n >= 10_000 else 5 if n >= 2_000 else 3
    q = np.quantile(x, np.linspace(0, 1, nb + 1)) if len(x) else np.linspace(0, 1, nb + 1)
    q = np.unique(q)
    if len(q) < 3:
        q = np.array([0.0, 0.5, 1.0])
    q[0], q[-1] = -np.inf, np.inf
    return q


def templates(mb, cm, fd, cf, D_mb, D_fd, sel_mb, sel_fd, w_mb=None):
    """Szablony D w komórce (sel_*: maski); zwraca (edges, tm, f_true, n_bkg)."""
    w = mb["w1"] if w_mb is None else w_mb
    R = sel_mb & (cm["O"] | cm["c27"] | cm["swap"])
    MMm = sel_mb & cm["MM"]
    Om = sel_mb & cm["O"]
    e = edges_for(D_mb[R | MMm], None, int((R | MMm).sum()))
    tm = {"MM": hist(D_mb[MMm], e, w[MMm]), "R": hist(D_mb[R], e, w[R]), "S": hist(D_fd[sel_fd & cf["sig"]], e, fd["w1"][sel_fd & cf["sig"]])}
    nR, nO, nMM = w[R].sum(), w[Om].sum(), w[MMm].sum()
    tm["c"] = float(1 - nO / nR) if nR > 0 else 0.0
    ftrue = float(nMM / (nMM + nO)) if nMM + nO > 0 else float("nan")
    eft = float(np.sqrt(ftrue * (1 - ftrue) / max(nMM + nO, 1)))
    return e, tm, ftrue, eft, int(MMm.sum()), int(R.sum())


def cell_masks(d, b, sel):
    lo, hi = PT[b]
    m = (d["pt"] >= lo) & (d["pt"] < hi)
    if sel != "all":
        m &= d["wp"][sel]
    return m


def data_classes(d, run):
    n = len(d["m"])
    tp = lambda c: (d[f"tp_{c}"] == 1) if f"tp_{c}" in d else np.zeros(n, bool)   # noqa: E731
    mb_any = np.zeros(n, bool)
    for c in MBTS:
        mb_any |= tp(c)
    j12 = tp("HLT_noalg_L1J12")
    grl = d["grl"] == 1
    if run == 427929:
        return {"427929 MB (mbts_L1MBTS_1∧GRL)": tp("HLT_mb_mbts_L1MBTS_1") & grl, "427929 GRL (wszystkie)": grl}
    return {"428770 all": np.ones(n, bool), "428770 J12only": j12 & ~mb_any, "428770 MBTSbased": mb_any}


# ------------------------------------------------------------------------------------------------------------------ główne
def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=os.path.join(here, "runs", "lxplus"))
    ap.add_argument("--out", default=here)
    ap.add_argument("--tag", default="", help="sufiks plików wyjściowych (np. _TEST)")
    ap.add_argument("--variants", default="twin,pid,all")
    a = ap.parse_args(argv)
    P = lambda s: os.path.join(a.inp, f"fmm_{s}.npz")   # noqa: E731
    for s in ("mb_nd", "fd_nd"):
        if not os.path.exists(P(s)):
            sys.exit(f"BŁĄD: brak {P(s)} (etap ORCH fmm nie pobrany?)")
    mb, meta = load(P("mb_nd"))
    fd, _ = load(P("fd_nd"))
    cm, cf = cats(mb), cats(fd)
    DATA = {}
    for s, run in (("data_428770", 428770), ("data_427929", 427929)):
        if os.path.exists(P(s)):
            d, mt = load(P(s))
            DATA[s] = (d, run, mt)
    log(f"mb {len(mb['m'])}, fd {len(fd['m'])}, dane {[(s, len(v[0]['m'])) for s, v in DATA.items()]}")
    res = {"meta": {"in": a.inp, "n_mb": int(len(mb["m"])), "n_fd": int(len(fd["m"])), "window": meta["window"], "lxy_min": meta["lxy_min"],
                    "n_data": {s: int(len(v[0]["m"])) for s, v in DATA.items()}, "frac": {s: v[2]["frac"] for s, v in DATA.items()},
                    "n_files": {"mb_nd": len(meta["files"]), **{s: len(v[2]["files"]) for s, v in DATA.items()}}}}
    bkg = cm["MM"] | cm["O"]

    # ---- 0. czym jest hypo_swap / MisMatch* (tabele prawdy MC, SR)
    sr = mb["sr"]
    tr = {}
    for k in ("sig", "swap", "c27", "MM", "O"):
        s = sr & cm[k]
        hs = mb["hypo_swap"][s]
        tr[k] = {"n": int(s.sum()), "hs-1": float((hs == -1).mean()) if s.any() else 0, "hs0": float((hs == 0).mean()) if s.any() else 0,
                 "hs1": float((hs == 1).mean()) if s.any() else 0, "tw_in": float(mb["tw_in"][s].mean()) if s.any() else 0,
                 "eps80": float(mb["wp"]["0.80"][s].mean()) if s.any() else 0}
    s = sr & cm["MM"]
    charm = lambda p: np.isin(np.abs(p), [411, 421, 431, 413, 423, 433, 4122])   # noqa: E731
    c0, c1 = charm(mb["tpar0"][s]), charm(mb["tpar1"][s])
    tr["MM_parent"] = {"D_tpdg==0": float((mb["D_tpdg"][s] == 0).mean()) if s.any() else 0, "one_charm_daughter": float((c0 ^ c1).mean()) if s.any() else 0,
                       "two": float((c0 & c1).mean()) if s.any() else 0, "none": float((~c0 & ~c1).mean()) if s.any() else 0}
    s = sr & cm["O"]
    c0, c1 = charm(mb["tpar0"][s]), charm(mb["tpar1"][s])
    tr["O_parent"] = {"one_charm_daughter": float((c0 ^ c1).mean()) if s.any() else 0, "none": float((~c0 & ~c1).mean()) if s.any() else 0}
    res["truth_sr"] = tr

    # ---- 1. f_MM prawda MC per bin pT × selekcja (SR, SB)
    res["f_mc"] = {}
    for b in PT:
        for sel in ["all"] + WPS:
            c = cell_masks(mb, b, sel)
            o = {}
            for rn in ("sr", "sb"):
                nMM, nB = int((c & mb[rn] & cm["MM"]).sum()), int((c & mb[rn] & bkg).sum())
                f = nMM / nB if nB else float("nan")
                o[rn] = {"f": f, "e": float(np.sqrt(f * (1 - f) / nB)) if nB else float("nan"), "n_bkg": nB, "n_MM": nMM}
            res["f_mc"][f"{b}|{sel}"] = o

    # ---- 2. E1 twin-rate
    res["E1"] = {}
    for b in PT:
        c = cell_masks(mb, b, "all")
        r = {k: float(mb["tw_in"][c & sr & cm[k]].mean()) if (c & sr & cm[k]).any() else float("nan") for k in ("MM", "O", "sig")}
        nMM, nO = int((c & sr & cm["MM"]).sum()), int((c & sr & cm["O"]).sum())
        slope = r["MM"] - r["O"]
        eslope = float(np.sqrt(r["MM"] * (1 - r["MM"]) / max(nMM, 1) + r["O"] * (1 - r["O"]) / max(nO, 1)))
        rsb = {rn: float(mb["tw_in"][c & mb[rn] & bkg].mean()) if (c & mb[rn] & bkg).any() else float("nan") for rn in ("sr", "sb")}
        ratio = (c & sr & bkg).sum() / max((c & mb["sb"] & bkg).sum(), 1)
        o = {"r_MM": r["MM"], "r_O": r["O"], "r_sig": r["sig"], "slope_pp_per_pp": slope, "e_slope": eslope, "r_mc_sr": rsb["sr"], "r_mc_sb": rsb["sb"], "data": {}}
        for sname, (d, run, _m) in DATA.items():
            for cn, cmk in data_classes(d, run).items():
                cc = cell_masks(d, b, "all") & cmk
                nsr, nsb = int((cc & d["sr"]).sum()), int((cc & d["sb"]).sum())
                if nsr < 200 or nsb < 200:
                    continue
                ns = max(nsr - nsb * ratio, 0.0)                         # sygnał: N_SR − N_SB·(SR/SB)_MC tła
                rsr_raw, rsbd = float(d["tw_in"][cc & d["sr"]].mean()), float(d["tw_in"][cc & d["sb"]].mean())
                rb = (rsr_raw * nsr - ns * r["sig"]) / max(nsr - ns, 1)
                e = float(np.sqrt(rsr_raw * (1 - rsr_raw) / nsr + rsbd * (1 - rsbd) / nsb))
                dloc = (rb - rsbd) - (rsb["sr"] - rsb["sb"])
                o["data"][cn] = {"r_sr_bkg": rb, "r_sb": rsbd, "ns_frac": ns / nsr, "d_r_glob": rb - rsb["sr"], "d_r_loc": dloc, "e_r": e,
                                 "df_glob_pp": 100 * (rb - rsb["sr"]) / slope if slope else float("nan"),
                                 "df_loc_pp": 100 * dloc / slope if slope else float("nan"), "e_df_pp": 100 * e / abs(slope) if slope else float("nan")}
        res["E1"][b] = o

    # ---- 3. E2 szablony D (warianty cech) + domknięcie
    res["E2"], res["closure"], res["auc"] = {}, {}, {}
    q = mb["q"]
    for v in a.variants.split(","):
        log(f"wariant {v}: trening")
        X = feats(mb, v)
        tr_m = (q == 0) & bkg & mb["veto"]
        if cm["MM"][tr_m].sum() < 50:
            log("  za mało MisMatch* do treningu — pomijam")
            continue
        y = cm["MM"][tr_m].astype(int)
        wt = np.where(y == 1, (1 - y.mean()) / max(y.mean(), 1e-9), 1.0)
        cl = train(X[tr_m], y, wt)
        D_mb = cl.predict_proba(X)[:, 1]
        D_fd = cl.predict_proba(feats(fd, v))[:, 1]
        te = (q != 0) & bkg & sr
        res["auc"][v] = {"auc_MM_vs_O_SR_test": auc(cm["MM"][te], D_mb[te]), "n_train": int(tr_m.sum()), "n_MM_train": int(y.sum())}
        log(f"  AUC(MM vs O, SR, test) = {res['auc'][v]['auc_MM_vs_O_SR_test']:.4f}")
        D_data = {s: cl.predict_proba(feats(d, v))[:, 1] for s, (d, run, _m) in DATA.items()}
        out = {}
        for b in PT:
            for sel in ["all"] + WPS:
                key = f"{b}|{sel}"
                cmb, cfd = cell_masks(mb, b, sel) & (q != 0), cell_masks(fd, b, sel)
                T = {}
                for rn in ("sr", "sb"):
                    T[rn] = templates(mb, cm, fd, cf, D_mb, D_fd, cmb & mb[rn], cfd & fd[rn])
                o = {"f_mc_sr": T["sr"][2], "e_f_mc_sr": T["sr"][3], "f_mc_sb": T["sb"][2], "e_f_mc_sb": T["sb"][3], "n_MM_tmpl": T["sr"][4], "data": {}}
                for sname, (d, run, _m) in DATA.items():
                    for cn, cmk in data_classes(d, run).items():
                        cc = cell_masks(d, b, sel) & cmk
                        if (cc & d["sr"]).sum() < 300 or (cc & d["sb"]).sum() < 300 or T["sr"][4] < 30 or T["sb"][4] < 30:
                            continue
                        r = {}
                        pr = sig_prior(float((cc & d["sr"]).sum()), float((cc & d["sb"]).sum()), T)
                        for rn in ("sr", "sb"):
                            e, tm = T[rn][0], T[rn][1]
                            h, hv = hist(D_data[sname][cc & d[rn]], e, d["w1"][cc & d[rn]])
                            r[rn] = tfit(h, hv, tm, with_sig=(rn == "sr"), prior=pr)
                        dg = r["sr"]["f"] - T["sr"][2]
                        dl = (r["sr"]["f"] - r["sb"]["f"]) - (T["sr"][2] - T["sb"][2])
                        o["data"][cn] = {"f_sr": r["sr"]["f"], "e_sr": r["sr"]["e"], "f_sb": r["sb"]["f"], "e_sb": r["sb"]["e"],
                                         "ns_frac": r["sr"]["ns"] / max(r["sr"]["n"], 1), "chi2_sr": r["sr"]["chi2"], "ndf_sr": r["sr"]["ndf"],
                                         "chi2_sb": r["sb"]["chi2"], "ndf_sb": r["sb"]["ndf"], "n_sr": r["sr"]["n"],
                                         "df_glob_pp": 100 * dg, "e_glob_pp": 100 * float(np.hypot(r["sr"]["e"], T["sr"][3])),
                                         "df_loc_pp": 100 * dl, "e_loc_pp": 100 * float(np.sqrt(r["sr"]["e"] ** 2 + r["sb"]["e"] ** 2 + T["sr"][3] ** 2 + T["sb"][3] ** 2))}
                out[key] = o
        res["E2"][v] = out

        # domknięcie: szablony q1, pseudo-dane q2+q3 (wszystkie kategorie mb), komórki 'all'
        log(f"  domknięcie {v}")
        cal, tst = q == 1, q >= 2
        g = np.exp(-0.5 * ((mb["m"].astype(float) - M0) / 0.035) ** 2)
        W = {"nominal": mb["w1"], "inject_local(MM×(1+0,5g))": np.where(cm["MM"], 1 + 0.5 * g, 1.0),
             "inject_global(MM×1,5)": np.where(cm["MM"], 1.5, 1.0)}
        if DATA:                                                      # stress: wagi dane/MC z klasyfikatora w SB (cechy 'all')
            sname = "data_428770" if "data_428770" in DATA else next(iter(DATA))
            dd = DATA[sname][0]
            Xa = feats(mb, "all")
            ms, ds = mb["sb"] & (q != 0) & (cm["MM"] | cm["O"] | cm["c27"]), dd["sb"]
            if ms.sum() > 2000 and ds.sum() > 2000:
                Xd = feats(dd, "all")[ds]
                if len(Xd) > 600_000:
                    Xd = Xd[RNG.choice(len(Xd), 600_000, replace=False)]
                Xm = Xa[ms & (q == 1)] if (ms & (q == 1)).sum() > 2000 else Xa[ms]     # uczone na q1 (kalibracja), stosowane do testu
                ww = np.r_[np.ones(len(Xd)), np.full(len(Xm), len(Xd) / len(Xm))]
                rw = train(np.r_[Xd, Xm], np.r_[np.ones(len(Xd)), np.zeros(len(Xm))].astype(int), ww, n=150)
                p = np.clip(rw.predict_proba(Xa)[:, 1], 1e-3, 1 - 1e-3)
                wst = np.clip(p / (1 - p), 0.05, 20.0)
                W[f"stress(wagi {sname}/MC z SB)"] = wst
                pt_ = rw.predict_proba(Xa[ms & tst])[:, 1]
                res["closure"].setdefault("_stress", {})["auc_data_vs_mc_SB"] = auc(np.r_[np.ones(min(len(Xd), 200_000)), np.zeros(len(pt_))],
                                                                                  np.r_[rw.predict_proba(Xd[:200_000])[:, 1], pt_])
                res["closure"]["_stress"]["sample"] = sname
        cl_out = {}
        for b in PT:
            cmb, cfd = cell_masks(mb, b, "all"), cell_masks(fd, b, "all")
            T = {rn: templates(mb, cm, fd, cf, D_mb, D_fd, cmb & cal & mb[rn], cfd & fd[rn]) for rn in ("sr", "sb")}
            if T["sr"][4] < 30 or T["sb"][4] < 30:
                continue
            for wn, w in W.items():
                r, ft = {}, {}
                pr = sig_prior(float(w[cmb & tst & mb["sr"]].sum()), float(w[cmb & tst & mb["sb"]].sum()), T)
                for rn in ("sr", "sb"):
                    s = cmb & tst & mb[rn]
                    h, hv = hist(D_mb[s], T[rn][0], w[s])
                    r[rn] = tfit(h, hv, T[rn][1], with_sig=(rn == "sr"), prior=pr)
                    ft[rn] = float(w[s & cm["MM"]].sum() / max(w[s & bkg].sum(), 1e-9))
                dl_hat = (r["sr"]["f"] - r["sb"]["f"]) - (T["sr"][2] - T["sb"][2])
                dl_true = (ft["sr"] - ft["sb"]) - (T["sr"][2] - T["sb"][2])
                cl_out[f"{b}|{wn}"] = {"f_true_sr": ft["sr"], "f_hat_sr": r["sr"]["f"], "e_sr": r["sr"]["e"], "f_true_sb": ft["sb"], "f_hat_sb": r["sb"]["f"],
                                       "e_sb": r["sb"]["e"], "bias_glob_pp": 100 * (r["sr"]["f"] - ft["sr"]), "e_glob_pp": 100 * r["sr"]["e"],
                                       "dloc_true_pp": 100 * dl_true, "dloc_hat_pp": 100 * dl_hat, "bias_loc_pp": 100 * (dl_hat - dl_true),
                                       "e_loc_pp": 100 * float(np.hypot(r["sr"]["e"], r["sb"]["e"]))}
        res["closure"][v] = cl_out

    # ---- 4. systematyka ε_B: Δf_loc (wariant 'all', 3–30 i per pT) × czułość B1; błąd = stat ⊕ |obciążenie stress|
    res["syst"] = {}
    vb = "all" if "all" in res["E2"] else (next(iter(res["E2"])) if res["E2"] else None)
    if vb:
        for b in PT:
            cell = res["E2"][vb].get(f"{b}|all", {}).get("data", {})
            st = [abs(x["bias_loc_pp"]) for k, x in res["closure"].get(vb, {}).items() if k.startswith(f"{b}|stress")]
            stg = [abs(x["bias_glob_pp"]) for k, x in res["closure"].get(vb, {}).items() if k.startswith(f"{b}|stress")]
            est = [x["e_loc_pp"] for k, x in res["closure"].get(vb, {}).items() if k.startswith(f"{b}|stress")]
            smod = max(st[0], est[0]) if st else float("nan")           # obciążenie niewidoczne poniżej błędu stat. testu → bierzemy większe
            for cn, x in cell.items():
                tot = float(np.sqrt(x["e_loc_pp"] ** 2 + (smod if np.isfinite(smod) else 0) ** 2))
                bound = float(np.hypot(x["df_loc_pp"], tot))
                res["syst"][f"{b}|{cn}"] = {"df_loc_pp": x["df_loc_pp"], "e_stat_pp": x["e_loc_pp"], "e_model_pp": smod, "e_tot_pp": tot,
                                            "df_glob_pp": x["df_glob_pp"], "e_glob_stat_pp": x["e_glob_pp"], "e_glob_model_pp": stg[0] if stg else float("nan"),
                                            "bound_pp": bound, "variant": vb,
                                            "eps": {sp: {wp: {"shift": SENS[sp][wp] * x["df_loc_pp"], "unc": abs(SENS[sp][wp]) * tot,
                                                              "bound": abs(SENS[sp][wp]) * bound, "assumed": ASSUMED[sp][wp]} for wp in WPS} for sp in SENS}}
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(os.path.join(a.out, f"results{a.tag}.json"), "w"), indent=1)
    report(res, os.path.join(a.out, f"RAPORT{a.tag}.md"), a)
    log("zapisano", os.path.join(a.out, f"RAPORT{a.tag}.md"))


# ------------------------------------------------------------------------------------------------------------------ raport
def F(x, n=2, s=False):
    if x is None or not np.isfinite(x):
        return "—"
    return f"{x:+.{n}f}" if s else f"{x:.{n}f}"


def report(R, path, a):
    L = []
    me = R["meta"]
    test = min(me["n_files"].values()) < 20
    L += ["# FMM — udział MisMatch* pod pikiem D0 z par bliźniaczych K↔π: prawda MC vs dane (projekt B)", "",
          f"Wygenerowane {time.strftime('%Y-%m-%d %H:%M')} przez `FMM_data/analyze.py` z `{me['in']}`. "
          f"MC mb_nd {me['n_mb']} wierszy ({me['n_files']['mb_nd']} plików), fd_nd {me['n_fd']}; dane: " +
          ", ".join(f"{s} {n} ({me['n_files'][s]} plików, frac {me['frac'][s]})" for s, n in me["n_data"].items()) +
          f". Okno {me['window'][0]}–{me['window'][1]} GeV, L_xy ≥ {me['lxy_min']} mm, pT 3–30, weto refleksji 25 MeV, SR 55 MeV, SB > 75 MeV.", ""]
    if test:
        L += ["> **BIEG TESTOWY na kilku plikach** — liczby służą wyłącznie sprawdzeniu kodu; błędy statystyczne są duże, wnioski ilościowe "
              "dopiero po pełnym biegu etapu ORCH `fmm`.", ""]
    t = R["truth_sr"]
    L += ["## 0. Ustalenia: czym jest `hypo_swap` i czy identyfikuje MisMatch*", "",
          "- `Dmeson_hypo_swap` to **flaga prawdy MC** (konwerter `charmpp_h5_v20.py`): 1, gdy przypisanie K/π kandydata jest odwrotne do "
          "`truthPdgId` córek, 0 gdy zgodne, −1 gdy para nie jest prawdziwym K+π. **W danych zawsze −1** — nie da się z niej nic „zmierzyć wprost”.",
          "- Ntupla zapisuje obie hipotezy mas tej samej pary torów jako **osobne wiersze** (bliźniaki), o ile obie masy mieszczą się w oknie; "
          "masa bliźniaka = `m_swap` (liczona analitycznie dla każdego kandydata, zgodna z zapisanym bliźniakiem do < 1 MeV).",
          "- MisMatch* = `truth_category` 8–11 (kod B1). To **nie** są zamiany K↔π: to pary, w których (zwykle) jeden tor pochodzi z rozpadu "
          "hadronu powabnego, a drugi jest przypadkowy. `hypo_swap` jest w nich określone równie rzadko/często jak w OtherNoB.", "",
          "| kategoria (MC mb, SR) | N | hypo_swap −1 / 0 / 1 | bliźniak w oknie | ε(WP80) |", "|---|---|---|---|---|"]
    nm = {"sig": "sygnał (0/1, hipoteza dobra)", "swap": "refleksja (0/1, hypo_swap = 1; po wecie)", "c27": "truth 2–7", "MM": "**MisMatch* (8–11)**", "O": "Other (12–15)"}
    for k in ("sig", "swap", "c27", "MM", "O"):
        x = t[k]
        L.append(f"| {nm[k]} | {x['n']} | {F(x['hs-1'])} / {F(x['hs0'])} / {F(x['hs1'])} | {F(x['tw_in'], 3)} | {F(x['eps80'], 3)} |")
    mp, op = t["MM_parent"], t["O_parent"]
    L += ["", f"Rodzice torów (prawda): MisMatch* — dokładnie jeden tor z hadronu powabnego {F(100 * mp['one_charm_daughter'], 1)} %, oba "
          f"{F(100 * mp['two'], 1)} %, żaden {F(100 * mp['none'], 1)} %; Other — jeden {F(100 * op['one_charm_daughter'], 1)} %, żaden {F(100 * op['none'], 1)} %.", ""]

    L += ["## 1. f_MM w MC (prawda) per bin pT i WP", "", "f_MM = N(truth 8–11)/N(truth ≥ 8), po wecie refleksji; „all” = przed cięciem na H.", "",
          "| pT | selekcja | f_MM SR [%] | f_MM SB [%] | N tła SR | N MM SR |", "|---|---|---|---|---|---|"]
    for k, x in R["f_mc"].items():
        b, sel = k.split("|")
        L.append(f"| {b} | {'all' if sel == 'all' else 'WP' + sel[2:]} | {F(100 * x['sr']['f'])} ± {F(100 * x['sr']['e'])} | {F(100 * x['sb']['f'])} ± {F(100 * x['sb']['e'])} | {x['sr']['n_bkg']} | {x['sr']['n_MM']} |")

    L += ["", "## 2. Estymator E1 — sama obserwabla bliźniacza (odsetek kandydatów z bliźniakiem w oknie)", "",
          "Kalibracja MC: r = r_O + f_MM·(r_MM − r_O). „Nachylenie” = zmiana r [pp] na 1 pp f_MM. Δf = Δr/nachylenie; błąd tylko statystyczny danych.", "",
          "| pT | r_MM | r_O | r_sygnał | nachylenie ± stat MC | próbka danych | r tła SR | r SB | Δf_glob [pp] | Δf_loc [pp] | ± stat [pp] |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for b, x in R["E1"].items():
        first = True
        for cn, y in (x["data"] or {"—": None}).items():
            head = f"| {b} | {F(x['r_MM'], 3)} | {F(x['r_O'], 3)} | {F(x['r_sig'], 3)} | {F(x['slope_pp_per_pp'], 4)} ± {F(x['e_slope'], 4)} |" if first else "| | | | | |"
            first = False
            L.append(head + (f" {cn} | {F(y['r_sr_bkg'], 4)} | {F(y['r_sb'], 4)} | {F(y['df_glob_pp'], 1, True)} | {F(y['df_loc_pp'], 1, True)} | {F(y['e_df_pp'], 1)} |" if y else " — | | | | | |"))
    L += ["", "Odczyt: r_MM i r_O różnią się o ułamek punktu procentowego na punkt f_MM, a sygnał ma r kilkakrotnie większe — obserwabla bliźniacza "
          "jest kinematyczna (asymetria pędów torów) i **nie rozróżnia MisMatch* od tła kombinatorycznego**; każde niedomodelowanie kinematyki "
          "lub resztka sygnału rzędu 0,1 pp w r przekłada się na wiele pp f_MM. E1 nie jest pomiarem.", ""]

    L += ["## 3. Estymator E2 — dopasowanie szablonów dyskryminanty D (MisMatch* vs reszta tła)", "",
          "| wariant cech | AUC(MM vs Other, SR, niezależna część MC) | N treningu (MM) |", "|---|---|---|"]
    for v, x in R["auc"].items():
        L.append(f"| {v} | {F(x['auc_MM_vs_O_SR_test'], 3)} | {x['n_train']} ({x['n_MM_train']}) |")
    L += ["", "f̂ = udział MisMatch* w tle z dopasowania 3 szablonów (sygnał z fd, MisMatch*, reszta) do rozkładu D w danych; błąd = bootstrap "
          "(dane + statystyka szablonów); sygnał z łagodnym ograniczeniem z liczenia SR/SB (σ = 50 % ⊕ 1 % N_SR); f̂ bez ograniczenia do [0, 1]. Δf_glob = f̂_SR(dane) − f_SR(MC); Δf_loc = (f̂_SR − f̂_SB)_dane − (f_SR − f_SB)_MC.", ""]
    for v, out in R["E2"].items():
        L += [f"### wariant `{v}`", "", "| pT | sel. | f_MM MC SR [%] | próbka danych | f̂ SR [%] | f̂ SB [%] | Δf_glob [pp] | Δf_loc [pp] | sygnał/N_SR | χ²/ndf SR |", "|---|---|---|---|---|---|---|---|---|---|"]
        for k, x in out.items():
            b, sel = k.split("|")
            for cn, y in x["data"].items():
                L.append(f"| {b} | {'all' if sel == 'all' else 'WP' + sel[2:]} | {F(100 * x['f_mc_sr'])} ± {F(100 * x['e_f_mc_sr'])} | {cn} | {F(100 * y['f_sr'])} ± {F(100 * y['e_sr'])} | "
                         f"{F(100 * y['f_sb'])} ± {F(100 * y['e_sb'])} | {F(y['df_glob_pp'], 2, True)} ± {F(y['e_glob_pp'])} | {F(y['df_loc_pp'], 2, True)} ± {F(y['e_loc_pp'])} | "
                         f"{F(y['ns_frac'], 4)} | {F(y['chi2_sr'], 1)}/{y['ndf_sr']} |")
        L.append("")

    L += ["## 4. Domknięcie estymatora na MC (kalibracja: trening q0 + szablony q1; test: q2+q3 jako pseudo-dane)", "",
          "„stress” = pseudo-dane przeważone wagami dane/MC (klasyfikator uczony w SB na cechach, bez m) — skład prawdziwy znany z prawdy; "
          "obciążenie = f̂ − f_prawda pokazuje, ile fałszywego Δf_MM produkuje samo niedomodelowanie cech."
          + (f" AUC dane-vs-MC w SB: {F(R['closure'].get('_stress', {}).get('auc_data_vs_mc_SB'), 3)}." if "_stress" in R["closure"] else ""), "",
          "| wariant | pT | test | f prawda SR [%] | f̂ SR [%] | obciążenie glob. [pp] | Δf_loc prawda [pp] | Δf_loc est. [pp] | obciążenie lok. [pp] |", "|---|---|---|---|---|---|---|---|---|"]
    for v, out in R["closure"].items():
        if v.startswith("_"):
            continue
        for k, x in out.items():
            b, wn = k.split("|")
            L.append(f"| {v} | {b} | {wn} | {F(100 * x['f_true_sr'])} | {F(100 * x['f_hat_sr'])} ± {F(100 * x['e_sr'])} | {F(x['bias_glob_pp'], 2, True)} ± {F(x['e_glob_pp'])} | "
                     f"{F(x['dloc_true_pp'], 2, True)} | {F(x['dloc_hat_pp'], 2, True)} | {F(x['bias_loc_pp'], 2, True)} ± {F(x['e_loc_pp'])} |")
    L.append("")

    L += ["## 5. Δf_MM (dane − MC) → systematyka ε_B przez czułości B1", "",
          "Wielkość wchodząca do czułości B1 (`cat_local`): lokalna zmiana składu pod pikiem, Δf_loc, selekcja „all”. Błąd całkowity = stat ⊕ model "
          "(model = max(|obciążenie lokalne w teście stress|, błąd stat. tego testu)). Ograniczenie |Δf| = √(Δf² + σ_tot²). "
          "Czułości [% ε_B na +1 pp]: D0 −1,0/−1,7/−4,6, D⁺ −3,3/−4,1/−6,4 (WP80/70/50). „Przyjęte” = dotychczasowa systematyka modelu tła "
          "(D0 ±5/±10/±15 %, D⁺ ±5/±10/±25 %). Δf_glob podane informacyjnie — SALAD koryguje skład globalny tylko częściowo (B1 cat_global).", "",
          "| pT | próbka | Δf_loc [pp] | ± stat | ± model | Δf_glob [pp] (± stat, ± model) | ogr. |Δf_loc| [pp] | D0: przesunięcie ε_B [%] @80/70/50 | D0: ogr. [%] | D0 przyjęte [%] | D⁺-czułości: ogr. [%] | D⁺ przyjęte [%] |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, x in R["syst"].items():
        b, cn = k.split("|")
        e0, e1 = x["eps"]["D0"], x["eps"]["D+"]
        L.append(f"| {b} | {cn} | {F(x['df_loc_pp'], 2, True)} | {F(x['e_stat_pp'])} | {F(x['e_model_pp'])} | {F(x['df_glob_pp'], 2, True)} (± {F(x['e_glob_stat_pp'])}, ± {F(x['e_glob_model_pp'])}) | {F(x['bound_pp'])} | "
                 + " / ".join(F(e0[w]["shift"], 1, True) for w in WPS) + " | " + " / ".join(F(e0[w]["bound"], 1) for w in WPS) + " | "
                 + " / ".join(F(e0[w]["assumed"], 0) for w in WPS) + " | " + " / ".join(F(e1[w]["bound"], 1) for w in WPS) + " | " + " / ".join(F(e1[w]["assumed"], 0) for w in WPS) + " |")
    L += ["", "## 6. Wnioski i zastrzeżenia", "",
          "1. Zdanie z notatki („f_MM in data can be measured directly from `hypo_swap`”) jest **nieścisłe**: `hypo_swap` to prawda MC, a MisMatch* "
          "nie jest kategorią zamian K↔π. Obserwabla bliźniacza (m_swap / bliźniak w oknie) ma AUC ≈ 0,5 dla MisMatch* vs Other (tabela §3) — "
          "**nie identyfikuje** MisMatch*.",
          "2. To, co da się zrobić: **ograniczenie** (przedział) na lokalną zmianę składu pod pikiem z dopasowania szablonów dyskryminanty "
          "wielocechowej (wariant `all`), zakotwiczonego w SB (Δf_loc). Liczba jest warunkowa: zakłada, że niedomodelowanie kształtów "
          "cech w poszczególnych kategoriach jest gładkie w m (to samo założenie co SALAD) — test „stress” mierzy, ile to kosztuje.",
          "3. Δf_glob (bez kotwicy SB) jest zdominowane przez niedomodelowanie cech (obciążenie w teście stress) — nie interpretować jako pomiaru składu.",
          "4. Wariant `all` używa cech H — nie jest niezależny od kontroli ε_B(SR) z fitu danych (B1 §2 pkt 4); warianty `twin`/`pid` są "
          "niezależne, ale mają znikomą moc. Zgodność wariantów w granicach błędów jest testem spójności, nie dowodem.",
          "5. f_MM per WP w danych (tabele §3) wymaga odjęcia sygnału o udziale rosnącym z WP; błąd rośnie odpowiednio. Sygnał/N_SR z dopasowania "
          "jest swobodny — kontrola: powinien być zgodny z fitem masy.",
          "6. Czułości B1 pochodzą z jednego wstrzyknięcia (+3,5 pp D0; +4,7 pp D⁺) — liniowość założona (brak skanu, METHODS_NOTE F5).",
          "7. 428770 jest zdominowany triggerem L1_J12 (inna kinematyka niż MC minbias) — klasy J12only/MBTSbased i run 427929 podane osobno; "
          "szablony MC nie są przeważane do klasy triggera.", ""]
    open(path, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
