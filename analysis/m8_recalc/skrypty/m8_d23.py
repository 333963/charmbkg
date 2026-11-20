#!/usr/bin/env python3
"""M8_przeliczenie pkt 2 — pik D0 w 2–3 GeV: stary model etapu 4 (kształt stały, exp2, 1,75–2,00) vs M8. POST HOC.
Wejście: D23_D0_weryfikacja_2026-11-19/in/*.npz (tylko odczyt). W MC NIE MA tła dla 2–3 GeV → szablonu DMM+CM dla tego binu nie ma.
Warianty: S4 (etap 4), M0 (kształt wolny), M8n (M8 bez składowej DMM+CM: Gauss + zakres 1,70–2,05),
          M8a (szablon i κ z 3–5 GeV wg przepisu M8: κ_T(ε ≤ 0,5) stałe, × 1,364), M8b (szablon 3–5; κ z prawdy mb_nd 3–5 dla tego samego klasyfikatora i ε, × 1,364).
  m8_d23.py → out/d23_m8.json"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, "..", "..", "FITBIAS_2026-11-19", "skrypty"))
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

IN = os.path.join(HERE, "..", "..", "D23_D0_weryfikacja_2026-11-19", "in")
M0, WIN, KSC, GC = 1.86484, 0.055, 1.364, (0.0008, 0.10)
NM = ["drzewo", "siec"]
HEAD = {("427929", 0): 0.1, ("427929", 1): 0.2, ("428770", 0): 0.1, ("428770", 1): 0.1}
YDP = {"427929": (3441.0, 260.0), "428770": (3055.0, 363.0)}          # y(D+) 2–3 GeV, D23 (H WP20)
FE = np.arange(1.70, 2.05 + F.BW / 2, F.BW)
SETS = {"S4": dict(shape="fix", rng=(1.75, 2.00), pk=None), "M0": dict(shape="free", rng=(1.75, 2.00), pk=None),
        "M8n": dict(shape="gauss", rng=(1.70, 2.05), pk=None), "M8a": dict(shape="gauss", rng=(1.70, 2.05), pk="a"), "M8b": dict(shape="gauss", rng=(1.70, 2.05), pk="b")}


def load(n):
    z = np.load(os.path.join(IN, n + ".npz"))
    return {k: z[k] for k in z.files}


def mc_wp(Fd, j, wp):
    sw = (Fd["cat"] == 0) & (np.abs(Fd["m"] - M0) < WIN)
    return float(np.quantile(Fd["lg"][sw, j].astype(np.float16).astype(np.float64), 1 - wp))


def templ(Fd, j, lo, hi=np.inf):
    w = Fd["w"].astype(np.float64)
    sig, rf = Fd["cat"] == 0, Fd["cat"] == 1
    s = (Fd["lg"][:, j] >= lo) & (Fd["lg"][:, j] < hi)
    ps, pr = sig & s, rf & s
    return dict(sigT=F.Tmpl(Fd["m"][ps], w[ps], 1.70, 2.05), refT=F.Tmpl(Fd["m"][pr], w[pr], 1.70, 2.05) if pr.sum() >= 300 else None,
                frefl=float(w[pr].sum() / w[ps].sum()) if pr.sum() >= 300 else 0.0, eps=float(w[ps].sum() / w[sig].sum()))


def kappa_truth(MB, F35, j, wp_hi, wp_lo=None):
    """κ = N(truth CM {2,3} + DMM {6,7}) / N(sygnał, poprawna hipoteza) w mb_nd 3–5 GeV, klasyfikator j, wycinek [WP hi, WP lo)."""
    lo = mc_wp(F35, j, wp_hi)
    hi = mc_wp(F35, j, wp_lo) if wp_lo else np.inf
    s = (MB["lg"][:, j] >= lo) & (MB["lg"][:, j] < hi)
    w = MB["w"].astype(np.float64)
    pk = np.isin(MB["truth"], (2, 3, 6, 7)) & s
    sg = (MB["cat"] == 0) & s
    k = float(w[pk].sum() / w[sg].sum())
    return k, int(pk.sum()), int(sg.sum())


def hist(m, w=None):
    return np.histogram(m, FE, weights=w)[0].astype(np.float64)


def cut(h, rng):
    i0, i1 = int(round((rng[0] - 1.70) / F.BW)), int(round((rng[1] - 1.70) / F.BW))
    return h[i0:i1], np.arange(rng[0], rng[1] + F.BW / 2, F.BW)


def fit(h, q, name, kap, peakT, start=None, hesse=True, fixN0=False):
    S = SETS[name]
    H, e = cut(h, S["rng"])
    tied = [(q["refT"], q["frefl"])] if q["refT"] is not None else []
    if S["pk"]:
        tied.append((peakT, kap[S["pk"]]))
    return F.fit(H[None], [q["sigT"]], [tied], e, M0, 2, shape=S["shape"], gc=GC, win=WIN, start=start, hesse=hesse, fixN0=fixN0)


def measure(d, s, q, kap, peakT, nboot, seed):
    rng = np.random.default_rng(seed)
    u, inv = np.unique(d["uid"][s], return_inverse=True)
    m = d["m"][s]
    h = hist(m)
    W = [rng.poisson(1.0, len(u)).astype(np.float64)[inv] for _ in range(nboot)]
    HB = [hist(m, w) for w in W]
    o = {}
    for name in SETS:
        f = fit(h, q, name, kap, peakT)
        f0 = fit(h, q, name, kap, peakT, fixN0=True)
        bN = []
        for hb in HB:
            fb = fit(hb, q, name, kap, peakT, start=f["x"], hesse=False)
            if fb["valid"]:
                bN.append(fb["N"][0])
        err_h = (f.get("err") or [float("nan")])[0]
        rms = float(np.std(bN)) if len(bN) > 5 else float("nan")
        err = float(np.nanmax([err_h, rms]))
        ctr = 0.5 * (FE[1:] + FE[:-1])
        o[name] = dict(N=f["N"][0], err_hesse=err_h, err_boot=rms, n_boot=len(bN), err=err, Z=f["N"][0] / err, ZLR=math.sqrt(max(f0["fval"] - f["fval"], 0.0)), chi2=f["chi2"], ndf=f["ndf"],
                       p=F.pval(f["chi2"], f["ndf"]), dmu=f["dmu"], sc=f["sc"], valid=f["valid"], at_limit=f["at_limit"], B=f["bkg_win"][0], eps=q["eps"], y=f["N"][0] / q["eps"], yerr=err / q["eps"],
                       kap=(kap[SETS[name]["pk"]] if SETS[name]["pk"] else 0.0), boot=bN)
    return o


def main():
    t0 = time.time()
    V._init()
    F23, F35, MB = load("fd_d23"), load("fd_nd_d150"), load("mb_nd_d150")
    D = {r: load(f"data_{r}_d23") for r in ("427929", "428770")}
    peakT = V.peak_tmpl("D0", 0)                                     # kształt DMM+CM: T, 3–5 GeV, WP80 (jak M8)
    kT = V.kappa_eps(V.kappa_nodes("D0", 0), 0.5)                    # κ_T(ε ≤ 0,5) = stałe poza węzłami (przepis M8)
    out = {"kappa_T_3-5_eps0.5": kT, "kscale": KSC, "head": {}, "wp": {}, "slice": {}, "kappa_truth_3-5": {}}
    print(f"kappa_T(3-5, eps<=0.5) = {kT:.3f} (x{KSC} = {kT * KSC:.3f})", flush=True)
    # punkty pracy skumulowane: WP30, WP20, WP10
    for run in D:
        d = D[run]
        for j in (0, 1):
            for wp in (0.3, 0.2, 0.1):
                th = mc_wp(F23, j, wp)
                q = templ(F23, j, th)
                kb, npk, nsg = kappa_truth(MB, F35, j, wp)
                out["kappa_truth_3-5"][f"{NM[j]}|{wp}"] = dict(kap=kb, n_peak=npk, n_sig=nsg)
                kap = {"a": kT * KSC, "b": kb * KSC}
                s = d["lg"][:, j] >= th
                head = HEAD[(run, j)] == wp
                o = measure(d, s, q, kap, peakT, 100 if head else 40, 11 + int(wp * 100) + 7 * j)
                out["wp"][f"{run}|{NM[j]}|{wp}"] = o
                print(f"[{time.time() - t0:4.0f}s] {run} {NM[j]} WP{int(wp * 100)}{' *' if head else ''} (eps {q['eps']:.4f}, kap_b {kb:.3f} [{npk}/{nsg}]): " +
                      "; ".join(f"{n}: N {v['N']:.0f} ± {v['err']:.0f} ZLR {v['ZLR']:.1f} p {v['p']:.2g} y {v['y']:.0f}" for n, v in o.items()), flush=True)
            # wycinki rozłączne (test 5 D23: średnia ważona 3 ciasnych wycinków ε < 0,3)
            for hi_wp, lo_wp in ((0.3, 0.2), (0.2, 0.1), (0.1, None)):
                lo = mc_wp(F23, j, hi_wp)
                hi = mc_wp(F23, j, lo_wp) if lo_wp else np.inf
                q = templ(F23, j, lo, hi)
                kb, npk, nsg = kappa_truth(MB, F35, j, hi_wp, lo_wp)
                kap = {"a": kT * KSC, "b": kb * KSC}
                s = (d["lg"][:, j] >= lo) & (d["lg"][:, j] < hi)
                o = measure(d, s, q, kap, peakT, 40, 5 + int(hi_wp * 100) + 7 * j)
                out["slice"][f"{run}|{NM[j]}|{hi_wp}-{lo_wp or 0.0}"] = o
            for name in SETS:
                ys = np.array([out["slice"][f"{run}|{NM[j]}|{a}-{b}"][name]["y"] for a, b in ((0.3, 0.2), (0.2, 0.1), (0.1, 0.0))])
                es = np.array([out["slice"][f"{run}|{NM[j]}|{a}-{b}"][name]["yerr"] for a, b in ((0.3, 0.2), (0.2, 0.1), (0.1, 0.0))])
                w = 1 / es ** 2
                mean, err = float((ys * w).sum() / w.sum()), float(1 / math.sqrt(w.sum()))
                c2 = float((((ys - mean) / es) ** 2).sum())
                yd, ed = YDP[run]
                r = mean / yd
                out["slice"][f"{run}|{NM[j]}|mean|{name}"] = dict(y=mean, yerr=err, chi2=c2, ndf=2, ratio_Dp=r, ratio_Dp_err=r * math.hypot(err / mean, ed / yd), ys=ys.tolist(), es=es.tolist())
                print(f"   <y> {run} {NM[j]} {name}: {mean:.0f} ± {err:.0f} (chi2 {c2:.1f}/2); y(D0)/y(D+) = {r:.2f} ± {r * math.hypot(err / mean, ed / yd):.2f}", flush=True)
    json.dump(out, open(os.path.join(HERE, "..", "out", "d23_m8.json"), "w"), indent=1, default=float)
    print(f"zapisano, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
