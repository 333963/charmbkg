"""HIGHSCORE krok 1b — ile deficytu > WP25 tłumaczy samo widmo pT(D0) wewnątrz binu 8–30 (MC przeważone do widma danych),
a ile zostaje. EKSPLORACYJNE, PO FAKCIE. Widmo danych: wydajki zliczeniem dla score >= WP_x w binach pT / wydajność MC tego cięcia
(założenie: SF(WP_x) płaskie w pT). Bootstrap po zdarzeniach obejmuje wagi pT i wydajki; MC bootstrap po parach.
python hs_ptw.py [NBOOT] -> out/ptw.json"""
import json, os, sys
import multiprocessing as mp
import numpy as np
import hs_common as C
R = C.R
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
EDG = {"drobne": [8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 23, 30], "grube": [8, 10, 12, 15, 20, 30]}
JOBS = [(428770, "0.70", "drobne"), (428770, "0.50", "drobne"), (428770, "0.80", "drobne"), (428770, "0.70", "grube"),
        (427929, "0.70", "grube"), (427929, "0.50", "grube")]
SUB = {"8-30": (8, 30), "8-12": (8, 12), "12-30": (12, 30)}


def job(a):
    run, wp, ename = a
    d = C.load(f"data_{run}"); M = C.load("fd_nd")
    mu = C.MU["8-30"]; W = 3.0
    ed = np.array(EDG[ename], float); nb = len(ed) - 1
    thr = C.thr_for("8-30", R.WP_MAIN); tw = C.THR["8-30"][wp]
    sd = C.base_mask(d, "8-30", isdata=True, run=run)
    ib, ok = C.dm_bin(d["dm"][sd]); sdi = np.where(sd)[0][ok]; ib = ib[ok]
    inv = d["inv"][sdi]; pt = d["pt"][sdi]; sc = d["score"][sdi]
    pb = np.clip(np.searchsorted(ed, pt, side="right") - 1, 0, nb - 1)
    selw = sc >= tw
    flat_w = (pb * R.NB + ib)[selw]
    cat = R.cat_index(sc, thr)
    flat_s = {k: ((cat * R.NB + ib)[(pt >= lo) & (pt < hi)], (pt >= lo) & (pt < hi)) for k, (lo, hi) in SUB.items()}
    sm = C.base_mask(M, "8-30", isdata=False) & (np.abs(M["dm"] - mu) < W)
    mpt = M["pt"][sm]; msc = M["score"][sm]
    mpb = np.clip(np.searchsorted(ed, mpt, side="right") - 1, 0, nb - 1); mcat = R.cat_index(msc, thr); mpass = msc >= tw

    def one(wd, wm, st):
        Hw = np.bincount(flat_w, weights=None if wd is None else wd[selw], minlength=nb * R.NB).reshape(nb, R.NB).astype(float)
        Sw, Pw, ncw = C.count_yields(Hw, mu, W, 3.0, st.get("ncw"), st.get("Pw"))
        st.setdefault("ncw", ncw); st.setdefault("Pw", Pw)
        n_all = np.bincount(mpb, weights=wm, minlength=nb); n_p = np.bincount(mpb[mpass], weights=wm[mpass], minlength=nb)
        Sd = np.clip(Sw, 0, None) / (n_p / n_all)
        wbin = (Sd / Sd.sum()) / (n_all / n_all.sum())
        wmc = wm * wbin[mpb]
        o = {"wbin": wbin}
        for k, (lo, hi) in SUB.items():
            fl, msk = flat_s[k]
            H = np.bincount(fl, weights=None if wd is None else wd[msk], minlength=6 * R.NB).reshape(6, R.NB).astype(float)
            S, P, nc = C.count_yields(H, mu, W, 3.0, st.get("nc" + k), st.get("P" + k))
            st.setdefault("nc" + k, nc); st.setdefault("P" + k, P)
            mm = (mpt >= lo) & (mpt < hi)
            f0 = np.bincount(mcat[mm], weights=wm[mm], minlength=6); f0 /= f0.sum()
            f1 = np.bincount(mcat[mm], weights=wmc[mm], minlength=6); f1 /= f1.sum()
            d0 = C.derived(S, f0); d1 = C.derived(S, f1)
            for q in ("rho50", "rho70", "SF50"):
                o[f"{k}|{q}|MC"] = d0[q]; o[f"{k}|{q}|MCptw"] = d1[q]; o[f"{k}|{q}|pred"] = d0[q] / d1[q]
            o[f"{k}|Rtop|MC"] = d0["R"][5]; o[f"{k}|Rtop|MCptw"] = d1["R"][5]; o[f"{k}|Rtop|pred"] = f1[5] / f0[5]
            o[f"{k}|R4|MC"] = d0["R"][4]; o[f"{k}|R4|MCptw"] = d1["R"][4]; o[f"{k}|R4|pred"] = f1[4] / f0[4]
        return o
    st = {}
    nom = one(None, np.ones(len(mpt)), st)
    rng = np.random.default_rng(31 + run)
    reps = [one(rng.poisson(1.0, d["nev"]).astype(float)[inv], rng.poisson(1.0, len(mpt)).astype(float), st) for _ in range(NB)]
    res = {"run": run, "wp": wp, "edges": EDG[ename], "ename": ename}
    for k in nom:
        v = np.array([r[k] for r in reps])
        res[k] = [np.asarray(nom[k]).tolist(), v.std(0, ddof=1).tolist()]
    return res


if __name__ == "__main__":
    C.load("fd_nd"); [C.load(f"data_{r}") for r in R.RUNS]
    with mp.get_context("fork").Pool(6) as P:
        out = P.map(job, JOBS, chunksize=1)
    json.dump(out, open(os.path.join(C.OUT, "ptw.json"), "w"), indent=1)
    for r in out:
        print(f"\n== {r['run']} wagi pT z >=WP{r['wp']} biny {r['ename']}: w = {np.round(r['wbin'][0],2).tolist()} ± {np.round(r['wbin'][1],2).tolist()}")
        for k in SUB:
            for q in ("Rtop", "R4", "rho50", "rho70", "SF50"):
                a, b, c = r[f"{k}|{q}|MC"], r[f"{k}|{q}|pred"], r[f"{k}|{q}|MCptw"]
                print(f"   pT {k:6s} {q:6s}: dane/MC = {a[0]:.3f}±{a[1]:.3f} | samo widmo pT (MCptw/MC) = {b[0]:.3f}±{b[1]:.3f} | dane/MCptw = {c[0]:.3f}±{c[1]:.3f}")
