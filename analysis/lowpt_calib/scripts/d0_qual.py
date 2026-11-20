"""LOWPT_calib pkt 4: JAKOŚCIOWY skan cięcia H dla inkluzywnego D0->Kpi (bez tagu D*, bez AFP) z odbiciem K<->pi.
Wejście: B_trig/runs/lxplus/bt_{data_*,fd_nd}.npz (pola all_*: m 1,75-2,00 GeV, score H, truth; truth 0/1 = sygnał, 100/101 = odbicie [zamiana hipotezy]).
Fit: 2 Gaussy (kształt z MC per WP) + szablon odbicia z MC (R = odbicie/sygnał ustalone z MC, warianty: R_MC, 0, 1,3 R_MC, R swobodne) + exp(pol2|pol3).
Błąd rho: wzór zagnieżdżony (Barlow) z hesjanów - bez bootstrapu (ocena jakościowa).  -> out/d0_qual.json"""
import json
import os
import sys

import numpy as np
from iminuit import Minuit
from scipy.special import ndtr

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
BT = os.path.join(os.path.dirname(TOP), "B_trig", "runs", "lxplus")
M0 = 1.86484
LO, HI, NB = 1.75, 2.00, 50
E = np.linspace(LO, HI, NB + 1)
X = np.linspace(LO, HI, 1001)
T = (X - M0) / 0.1
TP = np.stack([T ** (i + 1) for i in range(4)])
DX = np.diff(X)
BINS = [(3, 5), (5, 8), (8, 30)]
BN = ["3-5", "5-8", "8-30"]
WPS = [0.5, 0.6, 0.7, 0.8, 0.9]


def shape(m, w):
    e = np.linspace(M0 - 0.08, M0 + 0.08, 81)
    h, _ = np.histogram(m, e, weights=w)

    def nll(n, mu, s, k, f):
        lam = np.clip(n * np.diff(f * ndtr((e - mu) / s) + (1 - f) * ndtr((e - mu) / (s * k))), 1e-9, None)
        return float(np.sum(lam - h * np.log(lam)))
    mm = Minuit(nll, n=h.sum(), mu=M0, s=0.012, k=2.0, f=0.7)
    mm.errordef = 0.5
    mm.limits["s"] = (0.003, 0.04)
    mm.limits["k"] = (1.2, 6)
    mm.limits["f"] = (0.05, 1)
    mm.migrad()
    return {p: float(mm.values[p]) for p in ("mu", "s", "k", "f")}


def sigb(sh, dmu, sc):
    mu, s = sh["mu"] + dmu, sh["s"] * sc
    return np.diff(sh["f"] * ndtr((E - mu) / s) + (1 - sh["f"]) * ndtr((E - mu) / (s * sh["k"])))


def bkgb(b):
    ex = np.dot(np.asarray(b, float), TP[:len(b)])
    fb = np.exp(ex - ex.max())
    s = (0.5 * (fb[1:] + fb[:-1]) * DX).reshape(NB, 20).sum(1)
    return s / s.sum()


def fit(h, sh, refl, R, order, free_shape, free_R=False):
    names = ["N", "dmu", "sc", "nb", "R"] + [f"b{i+1}" for i in range(order)]

    def lam_of(p):
        return np.clip(p[0] * (sigb(sh, p[1], p[2]) + p[4] * refl) + p[3] * bkgb(p[5:]), 1e-9, None)

    def nll(*p):
        lam = lam_of(p)
        return float(np.sum(lam - h * np.log(lam)))
    c = 0.5 * (E[1:] + E[:-1])
    sr = np.abs(c - M0) < 0.03
    sb = np.abs(c - M0) > 0.07
    st = [max(h[sr].sum() - h[sb].mean() * sr.sum(), 10.0), 0.0, 1.0, h.sum(), R, -0.5] + [0.0] * (order - 1)
    m = Minuit(nll, *st, name=names)
    m.errordef = 0.5
    m.limits["nb"] = (0, None)
    m.limits["dmu"] = (-0.006, 0.006)
    m.limits["sc"] = (0.7, 1.5)
    m.limits["R"] = (0, 5)
    m.fixed["R"] = not free_R
    if not free_shape:
        m.fixed["dmu"] = True
        m.fixed["sc"] = True
    m.migrad(ncall=20000)
    if not m.valid:
        m.migrad(ncall=20000)
    m.hesse()
    p = [m.values[n] for n in names]
    lam = lam_of(p)
    return dict(N=float(p[0]), N_err=float(m.errors["N"]), R=float(p[4]), R_err=float(m.errors["R"]) if free_R else 0.0, valid=bool(m.valid),
                chi2=float(((h - lam) ** 2 / np.maximum(lam, 1)).sum()), ndf=int(NB - m.nfit), sc=float(p[2]))


def main():
    mc = np.load(os.path.join(BT, "bt_fd_nd.npz"))
    m, pt, s, tr, w = mc["all_m"].astype(float), mc["all_pt"].astype(float), mc["all_score"], mc["all_truth"], mc["all_w"].astype(float)
    sig, rf = tr <= 1, tr >= 100
    out = dict(mc={}, data={})
    thr, SH, RF, RR, EPS = {}, {}, {}, {}, {}
    for (lo, hi), bn in zip(BINS, BN):
        a = (pt >= lo) & (pt < hi)
        ss, ww = s[a & sig], w[a & sig]
        o = np.argsort(-ss)
        c = np.cumsum(ww[o]) / ww.sum()
        for wp in WPS:
            t = float(ss[o][min(np.searchsorted(c, wp), len(o) - 1)])
            thr[(bn, wp)] = t
            ps, pr = a & sig & (s >= t), a & rf & (s >= t)
            SH[(bn, wp)] = shape(m[ps], w[ps])
            hr, _ = np.histogram(m[pr], E, weights=w[pr])
            hs, _ = np.histogram(m[ps], E, weights=w[ps])
            RF[(bn, wp)] = hr / max(hr.sum(), 1e-9)
            RR[(bn, wp)] = float(hr.sum() / w[ps].sum())      # odbicie w oknie fitu / sygnał (cały)
            EPS[(bn, wp)] = float(w[ps].sum() / ww.sum())
            c0 = 0.5 * (E[1:] + E[:-1])
            win = np.abs(c0 - M0) < 0.04
            out["mc"][f"{bn}|{wp}"] = dict(thr=t, eps=EPS[(bn, wp)], R=RR[(bn, wp)], n_sig=int(ps.sum()), n_refl=int(pr.sum()),
                                          eps_refl=float(w[pr].sum() / max(w[a & rf].sum(), 1e-9)),
                                          refl_win_over_sig_win=float(hr[win].sum() / max(hs[win].sum(), 1e-9)),
                                          refl_rms=float(np.sqrt(np.average((m[pr] - M0) ** 2, weights=w[pr]))) if pr.sum() else 0.0)
    for run in ("427929", "428770"):
        d = np.load(os.path.join(BT, f"bt_data_{run}.npz"))
        sel = d["all_e_grl"] == 1
        if run == "427929":
            sel &= d["all_e_tp_HLT_mb_mbts_L1MBTS_1"] == 1
        dm, dpt, ds = d["all_m"][sel].astype(float), d["all_pt"][sel], d["all_score"][sel]
        for (lo, hi), bn in zip(BINS, BN):
            a = (dpt >= lo) & (dpt < hi)
            for wp in WPS:
                h, _ = np.histogram(dm[a & (ds >= thr[(bn, wp)])], E)
                h = h.astype(float)
                R = RR[(bn, wp)]
                res = {}
                for lab, RRv, order, fR in (("R_MC|pol3", R, 3, False), ("R_MC|pol2", R, 2, False), ("R0|pol3", 0.0, 3, False), ("R1.3|pol3", 1.3 * R, 3, False),
                                            ("Rfree|pol2", R, 2, True)):
                    res[lab] = fit(h, SH[(bn, wp)], RF[(bn, wp)], RRv, order, free_shape=(bn != "3-5"), free_R=fR)
                out["data"][f"{run}|{bn}|{wp}"] = dict(n=float(h.sum()), fits=res)
            print(run, bn, {wp: round(out["data"][f"{run}|{bn}|{wp}"]["fits"]["R_MC|pol3"]["N"] / EPS[(bn, wp)]) for wp in WPS}, flush=True)
    json.dump(out, open(os.path.join(TOP, "out", "d0_qual.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
