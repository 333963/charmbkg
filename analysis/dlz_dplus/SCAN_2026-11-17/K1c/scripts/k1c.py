"""K1c: fity nominalne (pol2/3/4, F-test), sparowany bootstrap po zdarzeniach, pseudo-dane (kandydaci MC + tło n_all·ε_B). KRYTERIA_K1c.md.
Wymaga out/prep.json (prep.py).  Użycie: python k1c.py [--nrep 200] [--ntoy 200] [--nproc 7]"""
import argparse
import json
import os
import sys
import time
import multiprocessing as _mp
Pool = _mp.get_context("fork").Pool

import numpy as np
from scipy.stats import f as FD

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prep as PR  # noqa: E402
FF = PR.FF
K1C = PR.K1C
RUNS, BINS, BN, TAGS, WPS = PR.RUNS, PR.BINS, PR.BN, PR.TAGS, PR.WPS
ORD = [2, 3, 4]
NB = FF.NB
G = {}   # globalne dla workerów (fork)


def chain(h, sh, free, st0=None):
    """fity rzędu 2,3,4 po kolei (start z niższego rzędu)."""
    out, prev = {}, st0
    for o in ORD:
        r = FF.fit(h, sh, order=o, free_shape=free, start=None if prev is None else prev)
        out[o] = r
        prev = r["par"]
    return out


def ftest(fits):
    ft = {}
    for o1, o2 in ((2, 3), (3, 4), (2, 4)):
        d = fits[o1]["chi2"] - fits[o2]["chi2"]
        F = d / (o2 - o1) / (fits[o2]["chi2"] / fits[o2]["ndf"])
        ft[f"{o1}->{o2}"] = dict(F=float(F), p=float(FD.sf(F, o2 - o1, fits[o2]["ndf"])) if F > 0 else 1.0)
    if ft["2->3"]["p"] < 0.05:
        sel = 4 if ft["3->4"]["p"] < 0.05 else 3
    else:
        sel = 4 if ft["2->4"]["p"] < 0.05 else 2
    return ft, sel


def nominal(a):
    run, bn, t, wp = a
    P = G["P"]
    h = np.asarray(P["hists"][f"{run}|{t}|{bn}|{wp}"], float)
    sh = P["shapes"][f"{t}|{bn}|{wp}"]
    out = {"lab": f"{run}|{bn}|{t}|{wp}"}
    for fv, free in (("free", True), ("fix", False)):
        fits = chain(h, sh, free)
        ft, sel = ftest(fits)
        for o in ORD:
            fits[o]["rel"] = fits[o]["N_err"] / max(fits[o]["N"], 1e-9)
        out[fv] = dict(fits={str(o): fits[o] for o in ORD}, ftest=ft, sel=sel)
    return out


def boot_job(a):
    run, bi, t, wp, fv, starts = a
    sh = G["P"]["shapes"][f"{t}|{BN[bi]}|{wp}"]
    hb = G["boot"][run][:, bi, list(TAGS).index(t), WPS.index(wp), :]
    res = {}
    for o in ORD:
        res[str(o)] = [FF.fit_fast(hb[i], sh, starts[str(o)], order=o, free_shape=(fv == "free")) for i in range(hb.shape[0])]
    return f"{run}|{BN[bi]}|{t}|{wp}", res


def eps_smooth(npass, nall, tails, gap, h=0.025):
    """ε_B(m): regresja lokalnie liniowa (jądro Gaussa h) udziału przechodzących w pasmach bocznych |Δm| > gap, po odjęciu ogonów sygnału."""
    c = FF.C
    sb = (np.abs(c - FF.M0) > gap) & (nall > 0)
    e = (npass - tails)[sb] / nall[sb]
    var = np.maximum(npass[sb], 1.0) / nall[sb] ** 2
    x = c[sb]
    out = np.zeros(NB)
    for i, x0 in enumerate(c):
        w = np.exp(-0.5 * ((x - x0) / h) ** 2) / var
        A = np.c_[np.ones(len(x)), x - x0]
        M = A.T @ (A * w[:, None])
        out[i] = np.linalg.solve(M, A.T @ (w * e))[0]
    return np.clip(out, 0, None)


def toy_job(a):
    run, bi, t, wp, fv, S0, ntoy, seed_sig, seed_bkg, do_v2_toys = a
    bn = BN[bi]
    P = G["P"]
    free = fv == "free"
    sh = P["shapes"][f"{t}|{bn}|{wp}"]
    sigexp = S0 * np.asarray(P["sig"][f"{t}|{bn}|{wp}"])
    ntrue = S0 * sh["eps_mc"]
    npass = np.asarray(P["hists"][f"{run}|{t}|{bn}|{wp}"], float)
    nall = np.asarray(P["hall"][f"{run}|{bn}"], float)
    out = {"lab": f"{run}|{bn}|{t}|{wp}", "N_true": float(ntrue), "sig_in_window_frac": float(sigexp.sum() / ntrue)}
    sbm = np.abs(FF.C - FF.M0) > 0.045
    win = ~sbm
    # sygnał: prawdziwi kandydaci MC, te same losowania dla wszystkich (klasyfikator, WP) w danym (run, bin)
    mcb = G["mc"][bi]
    rs = np.random.default_rng(seed_sig)
    ps = mcb["pass"][(t, wp)]
    lam = S0 * mcb["w"]
    sig_toys = np.zeros((ntoy, NB))
    for i in range(ntoy):
        cnt = rs.poisson(lam)
        sig_toys[i] = np.bincount(mcb["mb"][ps], weights=cnt[ps], minlength=NB + 2)[1:NB + 1]
    out["sig_toy_mean_over_true"] = float(sig_toys.sum(1).mean() / ntrue)
    rb = np.random.default_rng(seed_bkg)
    for vn, gap in (("V1", 0.045), ("V2", 0.070)):
        eb = eps_smooth(npass, nall, sigexp, gap)
        B = nall * eb
        d = npass - sigexp - B
        r = dict(B_win=float(B[win].sum()), SB_true=float(sigexp[win].sum() / max(B[win].sum(), 1e-9)),
                 chi2_sb=float(np.sum((d ** 2 / np.maximum(npass, 1))[sbm]) / sbm.sum()),
                 data_win_minus_SplusB=float((npass[win] - sigexp[win] - B[win]).sum()), fit={})
        asim = chain(sigexp + B, sh, free, st0=[ntrue, 0.0, 1.0, B.sum(), -0.5])
        dotoys = vn == "V1" or do_v2_toys
        if dotoys:
            toys = sig_toys + rb.poisson(B, (ntoy, NB))
        for o in ORD:
            q = dict(asimov=float(asim[o]["N"] / ntrue), asimov_valid=asim[o]["valid"])
            if dotoys:
                x = np.array([FF.fit_fast(toys[i], sh, asim[o]["par"], order=o, free_shape=free) for i in range(ntoy)])
                ok = (x[:, 4] > 0) & np.isfinite(x[:, 0]) & (x[:, 1] > 0)
                n = int(ok.sum())
                if n > 1:
                    rr = x[ok, 0] / ntrue
                    pull = (x[ok, 0] - ntrue) / x[ok, 1]
                    q.update(n_ok=n, f_ok=float(ok.mean()), b=float(rr.mean()), b_err=float(rr.std(ddof=1) / np.sqrt(n)), med=float(np.median(rr)), rms=float(rr.std(ddof=1)),
                             pull_mean=float(pull.mean()), pull_w=float(pull.std(ddof=1)), pull_w_err=float(pull.std(ddof=1) / np.sqrt(2 * n)),
                             relerr_hesse=float(np.mean(x[ok, 1]) / ntrue))
                else:
                    q.update(n_ok=n, f_ok=float(ok.mean()))
            r["fit"][str(o)] = q
        out[vn] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=200)
    ap.add_argument("--ntoy", type=int, default=200)
    ap.add_argument("--nproc", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()
    P = json.load(open(os.path.join(K1C, "out", "prep.json")))
    assert P["check"]["ALL_OK"], "warunek wstępny niespełniony"
    G["P"] = P
    cells = [(run, bn, t, wp) for run in RUNS for bn in BN for t in TAGS for wp in WPS]
    with Pool(a.nproc) as p:
        nom = {r["lab"]: r for r in p.imap_unordered(nominal, cells)}
    print(f"fity nominalne: {time.time()-t0:.0f} s", flush=True)
    fitvar = {}
    for run in RUNS:
        for bn in BN:
            fv = "fix" if bn == "3-5" else "free"
            if fv == "free":
                for t in TAGS:
                    for wp in WPS:
                        n = nom[f"{run}|{bn}|{t}|{wp}"]["free"]
                        if n["fits"][str(n["sel"])]["at_limit"]:
                            fv = "fix"
            fitvar[f"{run}|{bn}"] = fv
    out = {"nrep": a.nrep, "ntoy": a.ntoy, "nominal": nom, "fitvar": fitvar, "boot": {}, "toys": {}, "S0": {}}
    json.dump(out, open(os.path.join(K1C, "out", "k1c_nominal.json"), "w"))
    # --- sparowany bootstrap po zdarzeniach
    mc, scS, thr = PR.load_mc()
    G["boot"] = {}
    for run in RUNS:
        d = PR.load_data(run)
        ib = np.digitize(d["pt"], [3, 5, 8, 30]) - 1
        ue, inv = np.unique(d["ev"], return_inverse=True)
        mbin = np.clip(((d["m"] - FF.LO) / 0.005).astype(int), 0, NB - 1)
        idx = {}
        for bi, (lo, hi) in enumerate(BINS):
            for ti, t in enumerate(TAGS):
                for wi, wp in enumerate(WPS):
                    mk = np.flatnonzero((ib == bi) & (d["s"][t] >= thr[t][(lo, hi, wp)]))
                    idx[(bi, ti, wi)] = (mbin[mk], inv[mk])
                    h0 = np.bincount(mbin[mk], minlength=NB)
                    assert np.array_equal(h0, np.asarray(P["hists"][f"{run}|{t}|{BN[bi]}|{wp}"])), "histogram bootstrapu ≠ prep"
        rng = np.random.default_rng(20261117 + int(run))
        hb = np.zeros((a.nrep, 3, 2, len(WPS), NB))
        for r in range(a.nrep):
            w = rng.poisson(1.0, len(ue)).astype(float)
            for (bi, ti, wi), (mb_, ev_) in idx.items():
                hb[r, bi, ti, wi] = np.bincount(mb_, weights=w[ev_], minlength=NB)
        G["boot"][run] = hb
        del d
        print(f"repliki {run}: {time.time()-t0:.0f} s", flush=True)
    jobs = []
    for run in RUNS:
        for bi, bn in enumerate(BN):
            fv = fitvar[f"{run}|{bn}"]
            for t in TAGS:
                for wp in WPS:
                    n = nom[f"{run}|{bn}|{t}|{wp}"][fv]["fits"]
                    jobs.append((run, bi, t, wp, fv, {str(o): n[str(o)]["par"] for o in ORD}))
    # WP70 najpierw
    jobs.sort(key=lambda j: j[3] != 0.7)
    with Pool(a.nproc) as p:
        for lab, r in p.imap_unordered(boot_job, jobs, chunksize=1):
            out["boot"][lab] = r
    print(f"bootstrap: {time.time()-t0:.0f} s", flush=True)
    json.dump(out, open(os.path.join(K1C, "out", "k1c_boot.json"), "w"))
    G["boot"] = {}
    # --- pseudo-dane
    G["mc"] = {}
    mbm = np.digitize(mc["m"], FF.E)   # 0 = poniżej okna, NB+1 = powyżej
    for bi, (lo, hi) in enumerate(BINS):
        sel = (mc["pt"] >= lo) & (mc["pt"] < hi)
        G["mc"][bi] = dict(mb=mbm[sel], w=mc["w"][sel] / mc["w"][sel].sum(),
                           **{"pass": {(t, wp): (scS[t][sel] >= thr[t][(lo, hi, wp)]) for t in TAGS for wp in WPS}})
    jobs = []
    for ri, run in enumerate(RUNS):
        for bi, bn in enumerate(BN):
            fv = fitvar[f"{run}|{bn}"]
            nH = nom[f"{run}|{bn}|H|0.7"][fv]
            S0 = nH["fits"][str(nH["sel"])]["N"] / P["shapes"][f"H|{bn}|0.7"]["eps_mc"]
            out["S0"][f"{run}|{bn}"] = S0
            for t in TAGS:
                for wp in WPS:
                    jobs.append((run, bi, t, wp, fv, S0, a.ntoy, 777000 + 10 * ri + bi, 888000 + len(jobs), wp == 0.7))
    jobs.sort(key=lambda j: j[3] != 0.7)
    with Pool(a.nproc) as p:
        for r in p.imap_unordered(toy_job, jobs, chunksize=1):
            out["toys"][r["lab"]] = r
    print(f"pseudo-dane: {time.time()-t0:.0f} s", flush=True)
    json.dump(out, open(os.path.join(K1C, "out", "k1c.json"), "w"))
    print("zapisano", flush=True)


if __name__ == "__main__":
    main()
