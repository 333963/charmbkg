"""K1c etap 0: progi, kształty MC, histogramy danych po cięciu, szablony sygnału MC; WARUNEK WSTĘPNY: odtworzenie res_RUN przy WP70
(fit pol2, oba warianty kształtu) i zgodność fastfit ↔ fitlib. Zapis: out/prep.json.  Użycie: python prep.py"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K1C = os.path.dirname(HERE)
SCAN = os.path.dirname(K1C)
DATA = os.path.join(os.path.dirname(SCAN), "DATA_2026-11-17")
CAND = os.path.join(SCAN, "cand")
sys.path.insert(0, os.path.join(DATA, "scripts"))
sys.path.insert(0, HERE)
import fitlib as FL  # noqa: E402
import fastfit as FF  # noqa: E402

RUNS = ["427929", "428770"]
BINS = [(3, 5), (5, 8), (8, 30)]
BN = ["3-5", "5-8", "8-30"]
TAGS = {"H": "H_s12345", "DLZ": "HDLZ_s12345"}
WPS = [0.5, 0.6, 0.7, 0.8, 0.9]


def thresholds(s, w, pt, wps):
    """jak analyze.thresholds (kwantyle ważone OOF sygnału MC per bin pT)."""
    t = {}
    for lo, hi in BINS:
        a = (pt >= lo) & (pt < hi)
        ss, ww = s[a], w[a]
        o = np.argsort(-ss)
        c = np.cumsum(ww[o]) / ww.sum()
        for wp in wps:
            t[(lo, hi, wp)] = float(ss[o][min(np.searchsorted(c, wp), len(o) - 1)])
    return t


def load_mc():
    b = np.load(os.path.join(CAND, "train", "oof", "base.npz"))
    mc = dict(m=b["S_m"].astype(float), pt=b["S_pt"].astype(float), w=b["S_w"].astype(float))
    sc = {t: np.load(os.path.join(CAND, "train", "oof", TAGS[t] + ".npz"))["S"] for t in TAGS}
    thr = {t: thresholds(sc[t], mc["w"], mc["pt"], WPS + [0.95]) for t in TAGS}
    return mc, sc, thr


def load_data(run):
    d = np.load(os.path.join(CAND, f"dlz_data_{run}.npz"))
    sel = np.ones(len(d["m"]), bool)
    H1 = d["hall_any"]
    if run == "427929":
        sel = (d["grl"] == 1) & (d["mbts"] == 1)
        H1 = d["hall_mb"]
    return dict(m=d["m"][sel].astype(float), pt=d["pt"][sel].astype(float), ev=d["event"][sel], s={t: d["s_" + TAGS[t]][sel] for t in TAGS},
                hall=H1.reshape(3, 100, 5).sum(2).astype(float))


def main():
    t0 = time.time()
    mc, scS, thr = load_mc()
    tj = json.load(open(os.path.join(CAND, "train", "thr.json")))
    out = {"thr": {}, "shapes": {}, "sig": {}, "hists": {}, "hall": {}, "check": {}}
    dmax = 0.0
    for t in TAGS:
        for (lo, hi), bn in zip(BINS, BN):
            for wp in WPS + [0.95]:
                out["thr"][f"{t}|{bn}|{wp}"] = thr[t][(lo, hi, wp)]
                kk = f"{wp:.2f}"
                if kk in tj[TAGS[t]][bn]:
                    dmax = max(dmax, abs(tj[TAGS[t]][bn][kk] - thr[t][(lo, hi, wp)]))
    print(f"progi vs thr.json: max |Δ| = {dmax:.2e}", flush=True)
    out["check"]["thr_maxdiff"] = dmax
    for t in TAGS:
        for (lo, hi), bn in zip(BINS, BN):
            a = (mc["pt"] >= lo) & (mc["pt"] < hi)
            W = mc["w"][a].sum()
            for wp in WPS:
                ps = a & (scS[t] >= thr[t][(lo, hi, wp)])
                sh = FL.mc_shape(mc["m"][ps], mc["w"][ps])
                sh["eps_mc"] = float(mc["w"][ps].sum() / W)
                sh["n_mc"] = int(ps.sum())
                out["shapes"][f"{t}|{bn}|{wp}"] = sh
                hs, _ = np.histogram(mc["m"][ps], FL.E, weights=mc["w"][ps])
                out["sig"][f"{t}|{bn}|{wp}"] = (hs / W).tolist()   # oczekiwany histogram sygnału na jednostkę S0 (suma ≈ ε)
    ok_all = True
    for run in RUNS:
        d = load_data(run)
        res = json.load(open(os.path.join(DATA, "runs", "lxplus", "ana", f"res_{run}.json")))
        ib = np.digitize(d["pt"], [3, 5, 8, 30]) - 1
        print(f"run {run}: kandydatów {len(d['m'])} (res: n_loose {res['n_loose']}), zdarzeń {len(np.unique(d['ev']))} (res: {res['n_events_loose']})", flush=True)
        out["check"][f"n_{run}"] = [int(len(d["m"])), int(res["n_loose"])]
        for bi, bn in enumerate(BN):
            out["hall"][f"{run}|{bn}"] = d["hall"][bi].tolist()
            for t in TAGS:
                for wp in WPS:
                    lo, hi = BINS[bi]
                    h, _ = np.histogram(d["m"][(ib == bi) & (d["s"][t] >= thr[t][(lo, hi, wp)])], FL.E)
                    k = f"{run}|{t}|{bn}|{wp}"
                    out["hists"][k] = h.tolist()
                    kr = f"{TAGS[t]}|{bn}|{wp}"
                    if kr in res["hists"]:
                        same = bool(np.array_equal(h, np.asarray(res["hists"][kr])))
                        out["check"][f"hist_{k}"] = same
                        if not same:
                            ok_all = False
                            print(f"  HISTOGRAM RÓŻNY {k}: suma {h.sum()} vs {int(np.sum(res['hists'][kr]))}", flush=True)
                    if wp == 0.7:
                        sh = out["shapes"][f"{t}|{bn}|{wp}"]
                        for fv, free, F in (("free", True, res["fits"]), ("fix", False, res["fits_fixshape"])):
                            r0 = FL.fit(h, sh, free_shape=free)
                            r1 = FF.fit(h, sh, order=2, free_shape=free)
                            ref = F[kr]
                            dN, dE = r0["N"] / ref["N"] - 1, r0["N_err"] / ref["N_err"] - 1
                            fN, fE = r1["N"] / r0["N"] - 1, r1["N_err"] / r0["N_err"] - 1
                            good = abs(dN) < 0.005 and abs(dE) < 0.02
                            goodf = abs(fN) < 0.002 and abs(fE) < 0.02
                            ok_all &= good and goodf
                            out["check"][f"fit_{k}_{fv}"] = dict(N=r0["N"], N_err=r0["N_err"], N_res=ref["N"], N_err_res=ref["N_err"], dN=dN, dErr=dE, fast_dN=fN, fast_dErr=fE, ok=bool(good), fast_ok=bool(goodf))
                            print(f"  {k} {fv}: N={r0['N']:.1f}±{r0['N_err']:.1f} (res {ref['N']:.1f}±{ref['N_err']:.1f}; ΔN {100*dN:+.3f} %, Δσ {100*dE:+.2f} %)  fastfit ΔN {100*fN:+.4f} %, Δσ {100*fE:+.3f} %  {'OK' if good and goodf else 'NIEZGODNE'}", flush=True)
    out["check"]["ALL_OK"] = bool(ok_all)
    json.dump(out, open(os.path.join(K1C, "out", "prep.json"), "w"))
    print("WARUNEK WSTĘPNY:", "SPEŁNIONY" if ok_all else "NIESPEŁNIONY", f"({time.time()-t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
