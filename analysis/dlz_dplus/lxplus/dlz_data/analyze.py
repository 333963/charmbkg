"""DLZ_dplus/DATA: fity masy D+ w danych po cięciu H vs H+DLZ (vs H+dalitz), per bin pT × WP; zysk σ(N)/N, ΔSF_rel, pik,
rzeźbienie, sPlot (ε warunkowe, cechy, C2ST sygnału), połowy danych, ziarna. KRYTERIA.md.
Użycie: python analyze.py RUN [--nboot 200]     (RUN: 427929 | 428770 | 428770loc)
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitlib as FL  # noqa: E402

BASE = os.environ.get("DLZ_BASE") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BINS = [(3, 5), (5, 8), (8, 30)]
WPS = [0.95, 0.90, 0.80, 0.70, 0.50]
TAGS = ["H_s12345", "HDLZ_s12345", "Hdal_s12345"]
SEED2 = ["H_s23456", "HDLZ_s23456"]
M0 = FL.M0
bn = lambda b: f"{b[0]}-{b[1]}"  # noqa: E731


def thresholds(s, w, pt):
    t = {}
    for lo, hi in BINS:
        a = (pt >= lo) & (pt < hi)
        ss, ww = s[a], w[a]
        o = np.argsort(-ss)
        c = np.cumsum(ww[o]) / ww.sum()
        for wp in WPS:
            t[(lo, hi, wp)] = float(ss[o][min(np.searchsorted(c, wp), len(o) - 1)])
    return t


def load_mc():
    b = np.load(os.path.join(BASE, "oof", "base.npz"))
    mc = {"m": b["S_m"].astype(float), "pt": b["S_pt"].astype(float), "w": b["S_w"].astype(float), "key": b["S_key"]}
    thr, sc = {}, {}
    for t in TAGS + SEED2:
        f = os.path.join(BASE, "oof", t + ".npz")
        if os.path.exists(f):
            sc[t] = np.load(f)["S"]
            thr[t] = thresholds(sc[t], mc["w"], mc["pt"])
    return mc, sc, thr, b


def wauc(s1, w1, s0, w0):
    """ważone AUC (dopuszcza wagi ujemne): P(s1 > s0)."""
    s = np.r_[s1, s0]
    lab = np.r_[np.ones(len(s1), bool), np.zeros(len(s0), bool)]
    w = np.r_[w1, w0]
    o = np.argsort(s, kind="stable")
    lab, w = lab[o], w[o]
    c0 = np.cumsum(np.where(lab, 0, w))
    return float(np.sum(np.where(lab, w * c0, 0)) / (w1.sum() * w0.sum()))


def sweights(m, res, sh):
    """sWeights sygnału dla kandydatów o masach m, z gęstości nominalnego fitu (kształty ustalone, plony wolne)."""
    eps = 1e-4
    fs = (FL.sig_cdf(m + eps, sh, res["dmu"], res["sc"]) - FL.sig_cdf(m - eps, sh, res["dmu"], res["sc"])) / (2 * eps)
    fb = (FL.bkg_cdf(m + eps, res["b"]) - FL.bkg_cdf(m - eps, res["b"])) / (2 * eps)
    ns, nb = res["N"], res["nb"]
    den = ns * fs + nb * fb
    A = np.array([[np.sum(fs * fs / den**2), np.sum(fs * fb / den**2)], [np.sum(fs * fb / den**2), np.sum(fb * fb / den**2)]])
    V = np.linalg.inv(A)
    return (V[0, 0] * fs + V[0, 1] * fb) / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--nboot", type=int, default=200)
    a = ap.parse_args()
    run = a.run
    mc, scS, thr, base = load_mc()
    d = np.load(os.path.join(BASE, "slim", f"data_{run}.npz"), allow_pickle=False)
    sel = np.ones(len(d["m"]), bool)
    H1 = d["hall_any"]
    if run == "427929":
        sel = (d["grl"] == 1) & (d["mbts"] == 1)
        H1 = d["hall_mb"]
    m, pt, ev = d["m"][sel].astype(float), d["pt"][sel].astype(float), d["event"][sel]
    tags = [t for t in TAGS + SEED2 if f"s_{t}" in d.files and t in thr]
    sc = {t: d[f"s_{t}"][sel] for t in tags}
    ib = np.digitize(pt, [3, 5, 8, 30]) - 1
    out = {"run": run, "n_loose": int(sel.sum()), "n_events_loose": int(len(np.unique(ev))), "n_cand_all": int(H1.sum()),
           "n_files": int(d["n_files"]), "fits": {}, "shapes": {}}
    print(f"run {run}: kandydatów po weto {H1.sum()}, luźnych {sel.sum()}, zdarzeń (luźne) {out['n_events_loose']}", flush=True)

    def tarr(t, wp):
        return np.array([thr[t][(lo, hi, wp)] for lo, hi in BINS])[ib]
    # --- kształty MC i fity nominalne
    shapes, fits, hists = {}, {}, {}
    for t in tags:
        for bi, b in enumerate(BINS):
            a_mc = (mc["pt"] >= b[0]) & (mc["pt"] < b[1])
            for wp in WPS:
                ps = a_mc & (scS[t] >= thr[t][(b[0], b[1], wp)])
                sh = FL.mc_shape(mc["m"][ps], mc["w"][ps])
                sh["eps_mc"] = float(mc["w"][ps].sum() / mc["w"][a_mc].sum())
                sh["n_mc"] = int(ps.sum())
                h, _ = np.histogram(m[(ib == bi) & (sc[t] >= thr[t][(b[0], b[1], wp)])], FL.E)
                r = FL.fit(h, sh)
                k = (t, bn(b), wp)
                shapes[k], fits[k], hists[k] = sh, r, h
                print(f"{t:12s} {bn(b):5s} WP{int(wp*100)}: N={r['N']:.0f}±{r['N_err']:.0f} ({100*r['N_err']/max(r['N'],1e-9):.1f} %) S/B={r['SB']:.3f} Z={r['Z']:.1f} "
                      f"dmu={1e3*r['dmu']:+.2f}±{1e3*r['dmu_err']:.2f} MeV sc={r['sc']:.3f}±{r['sc_err']:.3f} chi2/ndf={r['chi2_ndf']:.2f} valid={r['valid']}", flush=True)
    # warianty modelu tła / okna (tylko nominalne tagi, WP70 i WP90)
    alt = {}
    e2 = FL.E[(FL.E >= 1.72 - 1e-9) & (FL.E <= 2.10 + 1e-9)]
    for t in tags[:3]:
        for b in BINS:
            for wp in WPS:
                k = (t, bn(b), wp)
                h = hists[k]
                i0 = int(round((1.72 - FL.LO) / 0.005))
                alt[k] = {"pol1": FL.fit(h, shapes[k], order=1), "pol3": FL.fit(h, shapes[k], order=3),
                          "win": FL.fit(h[i0:i0 + len(e2) - 1], shapes[k], edges=e2), "fixshape": FL.fit(h, shapes[k], free_shape=False)}
    # --- bootstrap po zdarzeniach (zbiór luźny: przechodzi WP95 dowolnego nominalnego tagu)
    li = np.arange(len(m))
    ue, inv = np.unique(ev[li], return_inverse=True)
    masks = {(t, bn(b), wp): (ib[li] == bi) & (sc[t][li] >= thr[t][(b[0], b[1], wp)]) for t in tags[:3] for bi, b in enumerate(BINS) for wp in WPS}
    mbin = np.clip(((m[li] - FL.LO) / 0.005).astype(int), 0, FL.NB - 1)
    rng = np.random.default_rng(20261117)
    rep = {k: np.zeros((a.nboot, 5)) for k in masks}
    repx = {k: np.zeros((a.nboot, 5)) for k in masks}   # kształt sygnału ustalony na MC (dμ = 0, skala = 1)
    start = {k: [fits[k]["N"], fits[k]["dmu"], fits[k]["sc"], fits[k]["nb"]] + fits[k]["b"] for k in masks}
    fitsx = {k: alt[k]["fixshape"] for k in masks}
    startx = {k: [fitsx[k]["N"], 0.0, 1.0, fitsx[k]["nb"]] + fitsx[k]["b"] for k in masks}
    for r in range(a.nboot):
        w = rng.poisson(1.0, len(ue)).astype(float)[inv]
        for k, mk in masks.items():
            h = np.bincount(mbin[mk], weights=w[mk], minlength=FL.NB)
            rep[k][r] = FL.fit_fast(h, shapes[k], start[k])
            repx[k][r] = FL.fit_fast(h, shapes[k], startx[k], free_shape=False)
        if r % 20 == 0:
            print("boot", r, flush=True)
    # połowy danych (parzyste/nieparzyste zdarzenia)
    half = {}
    for par in (0, 1):
        hs = (ev[li] % 2) == par
        for k, mk in masks.items():
            h = np.bincount(mbin[mk & hs], minlength=FL.NB)
            half[(par,) + k] = FL.fit(h, shapes[k])
            half[(par, "fix") + k] = FL.fit(h, shapes[k], free_shape=False)
    # --- wielkości pochodne
    def derive(FT, RP, HP):
        der = {}
        for b in BINS:
            for wp in WPS:
                kH = ("H_s12345", bn(b), wp)
                for t in tags[1:3]:
                    kD = (t, bn(b), wp)
                    rH, rD = FT[kH], FT[kD]
                    relH, relD = rH["N_err"] / max(rH["N"], 1e-9), rD["N_err"] / max(rD["N"], 1e-9)
                    ok = (RP[kH][:, 0] > 0) & (RP[kD][:, 0] > 0) & (RP[kH][:, 4] > 0) & (RP[kD][:, 4] > 0)
                    g_rep = 1 - (RP[kD][ok, 1] / RP[kD][ok, 0]) / (RP[kH][ok, 1] / RP[kH][ok, 0])
                    # błąd empiryczny N z rozrzutu replik
                    eH, eD = np.std(RP[kH][ok, 0]), np.std(RP[kD][ok, 0])
                    x = dict(rel_H=relH, rel_D=relD, gain=1 - relD / relH, gain_err=float(np.std(g_rep)), n_ok=int(ok.sum()),
                             relboot_H=float(eH / max(rH["N"], 1e-9)), relboot_D=float(eD / max(rD["N"], 1e-9)),
                             gain_boot=float(1 - (eD / max(rD["N"], 1e-9)) / (eH / max(rH["N"], 1e-9))),
                             Nratio=rD["N"] / max(rH["N"], 1e-9), Nratio_err=float(np.std(RP[kD][ok, 0] / RP[kH][ok, 0])),
                             ddmu=rD["dmu"] - rH["dmu"], ddmu_err=float(np.std(RP[kD][ok, 2] - RP[kH][ok, 2])),
                             scratio=rD["sc"] / rH["sc"], scratio_err=float(np.std(RP[kD][ok, 3] / RP[kH][ok, 3])))
                    for ref in (0.90, 0.95):
                        if wp >= ref:
                            continue
                        kHr, kDr = ("H_s12345", bn(b), ref), (t, bn(b), ref)
                        emH = shapes[kH]["eps_mc"] / shapes[kHr]["eps_mc"]
                        emD = shapes[kD]["eps_mc"] / shapes[kDr]["eps_mc"]
                        sfH = FT[kH]["N"] / max(FT[kHr]["N"], 1e-9) / emH
                        sfD = FT[kD]["N"] / max(FT[kDr]["N"], 1e-9) / emD
                        o2 = ok & (RP[kHr][:, 0] > 0) & (RP[kDr][:, 0] > 0)
                        sfH_r = RP[kH][o2, 0] / RP[kHr][o2, 0] / emH
                        sfD_r = RP[kD][o2, 0] / RP[kDr][o2, 0] / emD
                        tg = f"ref{int(ref*100)}"
                        x.update({f"SFrel_H_{tg}": sfH, f"SFrel_H_{tg}_err": float(np.std(sfH_r)), f"SFrel_D_{tg}": sfD, f"SFrel_D_{tg}_err": float(np.std(sfD_r)),
                                  f"dSF_{tg}": sfD / sfH, f"dSF_{tg}_err": float(np.std(sfD_r / sfH_r))})
                        # połowy
                        hh = []
                        for par in (0, 1):
                            a1 = half[HP(par) + kH]["N"] / max(half[HP(par) + kHr]["N"], 1e-9) / emH
                            a2 = half[HP(par) + kD]["N"] / max(half[HP(par) + kDr]["N"], 1e-9) / emD
                            hh.append(a2 / max(a1, 1e-9))
                        x[f"dSF_{tg}_halves"] = hh
                    der[(t, bn(b), wp)] = x
        return der
    der = derive(fits, rep, lambda par: (par,))
    derx = derive(fitsx, repx, lambda par: (par, "fix"))
    # --- rzeźbienie w danych: ε_B(m) w SB, R w paśmie 45–70 MeV z liniowego fitu w |Δm| > 70 MeV
    sculpt = {}
    c1 = FL.LO + 0.001 * (np.arange(500) + 0.5)
    i1 = np.clip(((m - FL.LO) / 0.001).astype(int), 0, 499)
    c20 = c1.reshape(25, 20).mean(1)
    for t in tags[:3]:
        for bi, b in enumerate(BINS + [(3, 30)]):
            inb = (ib == bi) if bi < 3 else np.ones(len(m), bool)
            n0_1 = (H1[bi] if bi < 3 else H1.sum(0)).astype(float)
            for wp in (0.8, 0.7):
                ok = inb & (sc[t] >= tarr(t, wp))
                n1_1 = np.bincount(i1[ok], minlength=500).astype(float)
                n0, n1 = n0_1.reshape(25, 20).sum(1), n1_1.reshape(25, 20).sum(1)
                e = n1 / np.maximum(n0, 1)
                err = np.sqrt(np.maximum(e * (1 - e), 1e-12) / np.maximum(n0, 1))
                sbm = (np.abs(c20 - M0) > 0.07) & (n0 > 20)
                A = np.c_[np.ones(sbm.sum()), c20[sbm] - M0]
                W = 1 / err[sbm] ** 2
                cov = np.linalg.inv(A.T @ (A * W[:, None]))
                beta = cov @ (A.T @ (W * e[sbm]))
                chi2 = float(np.sum(W * (e[sbm] - A @ beta) ** 2))
                band = (np.abs(c1 - M0) >= 0.045) & (np.abs(c1 - M0) < 0.07)
                nb0, nb1 = n0_1[band].sum(), n1_1[band].sum()
                obs = nb1 / nb0
                xm = np.array([1.0, np.sum(n0_1[band] * (c1[band] - M0)) / nb0])
                pred = float(xm @ beta)
                vp = float(xm @ cov @ xm)
                eo = np.sqrt(obs * (1 - obs) / nb0)
                sculpt[(t, bn(b), wp)] = dict(R=float(obs / pred), R_err=float(obs / pred * np.hypot(eo / obs, np.sqrt(vp) / pred)),
                                              slope=float(beta[1] * 0.1 / beta[0]), slope_err=float(np.sqrt(cov[1, 1]) * 0.1 / beta[0]),
                                              chi2_lin=chi2, ndf_lin=int(sbm.sum() - 2), eps=e.tolist(), eps_err=err.tolist(), n_band=int(nb0))
    # --- sPlot przy H WP90: ε warunkowe, cechy, C2ST sygnału
    splot = {}
    feat = d["feat"]
    fi = d["fi"]
    fnames = [str(x) for x in d["feat_names"]]
    pos = np.full(len(d["m"]), -1, np.int64)  # wiersz luźny → wiersz macierzy cech
    pos[fi] = np.arange(len(fi))
    pos = pos[sel]
    p90 = sc["H_s12345"] >= tarr("H_s12345", 0.90)
    sw = np.zeros(len(m))
    for bi, b in enumerate(BINS):
        k = ("H_s12345", bn(b), 0.90)
        q = p90 & (ib == bi)
        sw[q] = sweights(m[q], fits[k], shapes[k])
        splot[bn(b)] = {"sum_sw": float(sw[q].sum()), "N_fit": fits[k]["N"]}
    mc90 = scS["H_s12345"] >= np.array([thr["H_s12345"][(lo, hi, 0.90)] for lo, hi in BINS])[np.digitize(mc["pt"], [3, 5, 8, 30]) - 1]
    ibm = np.digitize(mc["pt"], [3, 5, 8, 30]) - 1
    qi = np.flatnonzero(p90)
    uq, invq = np.unique(ev[qi], return_inverse=True)
    bw = np.random.default_rng(7).poisson(1.0, (a.nboot, len(uq))).astype(np.float32)
    bwm_u, bwm_inv = np.unique(mc["key"][mc90], return_inverse=True)
    bwm = np.random.default_rng(8).poisson(1.0, (a.nboot, len(bwm_u))).astype(np.float32)
    for bi, b in enumerate(BINS):
        qb = ib[qi] == bi
        mb_ = mc90 & (ibm == bi)
        mbl = (ibm[mc90] == bi)
        for t in tags[:3]:
            for wp in (0.8, 0.7, 0.5):
                pd_ = sc[t][qi] >= thr[t][(b[0], b[1], wp)]
                ed = (sw[qi] * pd_)[qb].sum() / sw[qi][qb].sum()
                edr = np.array([(bw[r][invq] * sw[qi] * pd_)[qb].sum() / (bw[r][invq] * sw[qi])[qb].sum() for r in range(a.nboot)])
                pm = scS[t] >= thr[t][(b[0], b[1], wp)]
                em = (mc["w"] * pm)[mb_].sum() / mc["w"][mb_].sum()
                wm90 = mc["w"][mc90]
                pm90 = pm[mc90]
                emr = np.array([(bwm[r][bwm_inv] * wm90 * pm90)[mbl].sum() / (bwm[r][bwm_inv] * wm90)[mbl].sum() for r in range(a.nboot)])
                splot[bn(b)][f"{t}_WP{int(wp*100)}"] = dict(eps_data=float(ed), eps_data_err=float(np.std(edr)), eps_mc=float(em), eps_mc_err=float(np.std(emr)),
                                                           sf=float(ed / em), sf_err=float(ed / em * np.hypot(np.std(edr) / ed, np.std(emr) / em)))
    # w binach dMinPt (pT-integralnie, progi per bin pT), H+DLZ WP70 | H WP90
    j = fnames.index("dMinPt_cn")
    S_feat = None
    fS = os.path.join(BASE, "oof", "S_feat.npz")
    if os.path.exists(fS):
        S_feat = np.load(fS)["feat"]
    if S_feat is not None:
        xd = feat[pos[qi], j]
        xm = S_feat[mc90, j]
        edges = np.quantile(xm, [0, 1 / 3, 2 / 3, 1])
        edges[0], edges[-1] = -np.inf, np.inf
        splot["dMinPt_bins"] = {"edges": [float(x) for x in np.quantile(xm, [1 / 3, 2 / 3])]}
        for t in ("HDLZ_s12345", "H_s12345"):
            if t not in sc:
                continue
            pd_ = sc[t][qi] >= tarr(t, 0.70)[qi]
            pm_ = (scS[t] >= np.array([thr[t][(lo, hi, 0.70)] for lo, hi in BINS])[ibm])[mc90]
            # MC przeważone do plonów sWeights per bin pT
            wm = mc["w"][mc90].copy()
            for bi in range(3):
                wm[ibm[mc90] == bi] *= max(sw[qi][ib[qi] == bi].sum(), 0) / wm[ibm[mc90] == bi].sum()
            rows = []
            for q in range(3):
                dq = (xd >= edges[q]) & (xd < edges[q + 1])
                mq = (xm >= edges[q]) & (xm < edges[q + 1])
                ed = (sw[qi] * pd_)[dq].sum() / sw[qi][dq].sum()
                edr = np.array([(bw[r][invq] * sw[qi] * pd_)[dq].sum() / (bw[r][invq] * sw[qi])[dq].sum() for r in range(a.nboot)])
                em = (wm * pm_)[mq].sum() / wm[mq].sum()
                rows.append(dict(eps_data=float(ed), eps_data_err=float(np.std(edr)), eps_mc=float(em), sf=float(ed / em), sf_err=float(np.std(edr) / em),
                                 frac_data=float(sw[qi][dq].sum() / sw[qi].sum()), frac_mc=float(wm[mq].sum() / wm.sum())))
            splot["dMinPt_bins"][t] = rows
        # rozkłady 1D 8 cech + score: χ² sWeighted dane vs MC (10 kwantyli MC), MC przeważone do plonów per bin pT
        wm = mc["w"][mc90].copy()
        for bi in range(3):
            wm[ibm[mc90] == bi] *= max(sw[qi][ib[qi] == bi].sum(), 0) / wm[ibm[mc90] == bi].sum()
        splot["feat1d"] = {}
        for f in fnames:
            jj = fnames.index(f)
            xm_, xd_ = S_feat[mc90, jj], feat[pos[qi], jj]
            ed_ = np.quantile(xm_, np.linspace(0, 1, 11))
            ed_[0], ed_[-1] = -np.inf, np.inf
            hd = np.array([sw[qi][(xd_ >= ed_[i]) & (xd_ < ed_[i + 1])].sum() for i in range(10)])
            hde = np.sqrt(np.array([(sw[qi] ** 2)[(xd_ >= ed_[i]) & (xd_ < ed_[i + 1])].sum() for i in range(10)]))
            hm = np.array([wm[(xm_ >= ed_[i]) & (xm_ < ed_[i + 1])].sum() for i in range(10)])
            hme = np.sqrt(np.array([(wm ** 2)[(xm_ >= ed_[i]) & (xm_ < ed_[i + 1])].sum() for i in range(10)]))
            hm_s = hm * hd.sum() / hm.sum()
            hme_s = hme * hd.sum() / hm.sum()
            chi2 = float(np.sum((hd - hm_s) ** 2 / (hde ** 2 + hme_s ** 2)))
            splot["feat1d"][f] = dict(chi2=chi2, ndf=9, ratio=(hd / hm_s).tolist(), ratio_err=(np.sqrt(hde ** 2 + hme_s ** 2) / hm_s).tolist(),
                                      mean_data=float((sw[qi] * xd_).sum() / sw[qi].sum()), mean_mc=float((wm * xm_).sum() / wm.sum()))
        # C2ST sygnału: dane (sWeights, H WP90) vs MC (H WP90), zestaw H vs H+DLZ; 5 foldów po zdarzeniach
        import xgboost as xgb
        res_c = {}
        Xd, Xm = feat[pos[qi]], S_feat[mc90]
        fd_ = (ev[qi] % 5).astype(int)
        fm_ = (mc["key"][mc90] % 5).astype(int)
        swq = sw[qi]
        for nm, cols in (("H", list(range(18))), ("HDLZ", list(range(26))), ("Hdal", list(range(20)))):
            aucs = []
            for seed in (0, 1, 2):
                sd, sm = np.zeros(len(Xd), np.float32), np.zeros(len(Xm), np.float32)
                for k in range(5):
                    trd, trm = fd_ != k, fm_ != k
                    X = np.vstack([Xd[trd][:, cols], Xm[trm][:, cols]])
                    y = np.r_[np.ones(trd.sum()), np.zeros(trm.sum())]
                    swt = np.clip(swq[trd], 0, None)  # XGB nie przyjmuje wag ujemnych: trening na części dodatniej sWeights, ocena AUC z pełnymi sWeights
                    w = np.r_[swt / swt.sum(), wm[trm] / wm[trm].sum()] * len(y) / 2
                    cl = xgb.XGBClassifier(n_estimators=150, max_depth=3, learning_rate=0.1, subsample=0.8, colsample_bytree=0.8, min_child_weight=20,
                                           n_jobs=int(os.environ.get("DLZ_NJ", "2")), random_state=seed, tree_method="hist")
                    cl.fit(X, y, sample_weight=w)
                    sd[~trd] = cl.predict_proba(Xd[~trd][:, cols])[:, 1]
                    sm[~trm] = cl.predict_proba(Xm[~trm][:, cols])[:, 1]
                aucs.append(wauc(sd, swq, sm, wm))
                if seed == 0:
                    br = [wauc(sd, swq * bw[r][invq], sm, wm * bwm[r][bwm_inv]) for r in range(min(a.nboot, 100))]
                    res_c[nm + "_booterr"] = float(np.std(br))
                    res_c[nm + "_bootrep"] = br
            res_c[nm] = aucs
            print("C2ST sygnału", nm, aucs, res_c[nm + "_booterr"], flush=True)
        if "H_bootrep" in res_c:
            for nm in ("HDLZ", "Hdal"):
                res_c[f"d_{nm}_minus_H"] = float(np.mean(res_c[nm]) - np.mean(res_c["H"]))
                res_c[f"d_{nm}_minus_H_err"] = float(np.std(np.array(res_c[nm + "_bootrep"]) - np.array(res_c["H_bootrep"])))
        splot["c2st_signal"] = {k: v for k, v in res_c.items() if not k.endswith("bootrep")}
    # --- zapis
    ks = lambda k: "|".join(str(x) for x in k)  # noqa: E731
    out["fits"] = {ks(k): {q: v for q, v in r.items() if q != "lam"} for k, r in fits.items()}
    out["lam"] = {ks(k): r["lam"] for k, r in fits.items()}
    out["hists"] = {ks(k): h.tolist() for k, h in hists.items()}
    out["shapes"] = {ks(k): v for k, v in shapes.items()}
    out["alt"] = {ks(k): {n: {q: r[q] for q in ("N", "N_err", "chi2_ndf", "valid", "Z")} for n, r in v.items()} for k, v in alt.items()}
    out["at_limit"] = {ks(k): r["at_limit"] for k, r in fits.items()}
    out["derived"] = {ks(k): v for k, v in der.items()}
    out["derived_fixshape"] = {ks(k): v for k, v in derx.items()}
    out["fits_fixshape"] = {ks(k): {q: v for q, v in r.items() if q != "lam"} for k, r in fitsx.items()}
    out["sculpt"] = {ks(k): v for k, v in sculpt.items()}
    out["splot"] = splot
    out["half"] = {ks(k): {q: r[q] for q in ("N", "N_err", "valid")} for k, r in half.items()}
    out["thr"] = {t: {ks(k): v for k, v in thr[t].items()} for t in thr}
    out["boot_ok_frac"] = {ks(k): float(np.mean(v[:, 4])) for k, v in rep.items()}
    hall = {bn(b): H1[bi].reshape(100, 5).sum(1).tolist() for bi, b in enumerate(BINS)}
    out["hist_all"] = hall
    os.makedirs(os.path.join(BASE, "out"), exist_ok=True)
    json.dump(out, open(os.path.join(BASE, "out", f"res_{run}.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(BASE, "out", f"rep_{run}.npz"), **{ks(k): v for k, v in rep.items()}, **{"fix|" + ks(k): v for k, v in repx.items()})
    print("zapisano", flush=True)


if __name__ == "__main__":
    main()
