#!/usr/bin/env python3
"""M8_przeliczenie pkt 3 — metoda R (sieć × drzewo, kategorie rozłączne b/n/t) dla D0: domknięcie na MC starym modelem (etap 4) i M8. POST HOC.
Wejście: in/app_d150/app_{fd_nd,mb_nd}_*.npz (pobrane z EOS /eos/user/k/klyczak/charm_q/dnet4/app/d150, tylko MC; bez jobów).
dn4_ana (Q_lxplus/dnet4) i fb_fit / fb_var (FITBIAS) importowane tylko do odczytu.
M8 dla kategorii łącznych (ROZSZERZENIE — przepis M8 nie definiuje κ dla kategorii sieć × drzewo): ograniczenie Gaussa (0,8 MeV; 0,10) wspólne dla 3 kategorii,
zakres 1,70–2,05, szablon DMM+CM = T (WP80 danego binu, jak M8), κ_c = N(truth CM {2,3} + DMM {6,7})/N(sygnał) w kategorii c z prawdy mb_nd (na MC bez skali BR).
Wariant M8n: to samo bez składowej DMM+CM. Bramka etapu 4: fit określony, |R_fit − R_prawda| ≤ 2,5σ i ≤ 0,04, σ ≤ 0,05; σ = max(HESSE, bootstrap po zdarzeniach).
  m8_R.py [--nboot 100] [--tags s1,s2,s3] → out/R_m8.json"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, "..", "..", "FITBIAS_2026-11-19", "skrypty"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "Q_lxplus", "dnet4"))
import dn4_ana as A  # noqa: E402
import fb_fit as F  # noqa: E402
import fb_var as V  # noqa: E402

APP = os.path.join(HERE, "..", "in", "app_d150")
SP, REF, GC = "D0", "D0.xgb.H.sb.n.s1", (0.0008, 0.10)
G = {}


def init():
    V._init()
    G["F"] = A.load_app(APP, "fd_nd", SP)
    G["M"] = A.load_app(APP, "mb_nd", SP)


def fit_m(H, mc, kap, peakT, model, start=None, hesse=True):
    c = A.SP[SP]
    if model == "M0":
        rng = c["fit"]
        tied = [[(mc.refT[j], mc.frefl[j])] if mc.refT[j] is not None else [] for j in range(3)]
        shape = "free"
    else:
        rng = c["alt"]
        tied = []
        for j in range(3):
            t = [(mc.refT[j], mc.frefl[j])] if mc.refT[j] is not None else []
            if model == "M8":
                t.append((peakT, kap[j]))
            tied.append(t)
        shape = "gauss"
    i0, i1 = int(round((rng[0] - 1.70) / F.BW)), int(round((rng[1] - 1.70) / F.BW))
    e = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    f = F.fit(H[:, i0:i1], mc.sigT, tied, e, c["m0"], c["nbkg"], shape=shape, gc=GC, win=c["win"], start=start, hesse=hesse)
    f.pop("mu", None)
    return f


def task(t):
    tag, b, wp, nboot = t
    Fd, M = G["F"], G["M"]
    c = A.SP[SP]
    bn = A.bname(b)
    mc = A.MC(SP, Fd, tag, REF, b, wp, (1.70, 2.05))
    inb = M["ok"] & (M["truth"] < 100) & (M["pt"] >= b[0]) & (M["pt"] < b[1])
    cat = mc.cat(M["lg"][inb, M["tags"].index(tag)], M["lg"][inb, M["tags"].index(REF)])
    m, uid, tr, sw = M["m"][inb].astype(np.float64), M["uid"][inb], M["truth"][inb], M["swap"][inb]
    fe = np.arange(1.70, 2.05 + F.BW / 2, F.BW)
    H, ok, idx = A.hist_cats(m, cat, fe)
    ts = np.isin(tr, (0, 1)) & (sw != 1)
    pk = np.isin(tr, (2, 3, 6, 7))
    Tn = [int((ts & (cat == j)).sum()) for j in range(3)]
    Pn = [int((pk & (cat == j)).sum()) for j in range(3)]
    kap = [Pn[j] / max(Tn[j], 1) for j in range(3)]
    peakT = V.peak_tmpl(SP, A.BINS.index(b))
    rt = A.r_of(Tn, None, mc.e)
    Rt = rt["R"]
    Rt_err = Rt * math.sqrt(Tn[1] / max(Tn[0] + Tn[1], 1) ** 2 + Tn[2] / max(Tn[0] + Tn[2], 1) ** 2 + mc.eps_ratio_relerr ** 2)
    res = {"tag": tag, "bin": bn, "wp": wp, "eps": mc.e, "eps_net": mc.eps_net, "eps_tree": mc.eps_tree, "thr": mc.thr, "truth_N": Tn, "peak_N": Pn, "kap": kap, "frefl": mc.frefl,
           "R_truth": Rt, "R_truth_err": Rt_err, "n_rows": [float(x) for x in H.sum(1)],
           "kap_cum": {"net": (Pn[0] + Pn[1]) / (Tn[0] + Tn[1]), "tree": (Pn[0] + Pn[2]) / (Tn[0] + Tn[2])}}
    rng = np.random.default_rng(4 + 17 * A.BINS.index(b) + int(wp * 100))
    u, inv = np.unique(uid[ok], return_inverse=True)
    HB = []
    for _ in range(nboot):
        wv = rng.poisson(1.0, len(u)).astype(np.float64)[inv]
        HB.append(np.bincount(idx, weights=wv, minlength=H.size).reshape(H.shape))
    for model in ("M0", "M8n", "M8"):
        f = fit_m(H, mc, kap, peakT, model)
        o = {"fit": {k: f[k] for k in ("N", "dmu", "sc", "valid", "chi2", "ndf", "at_limit")}, "p": F.pval(f["chi2"], f["ndf"])}
        r = A.r_of(f["N"], f.get("cov_N"), mc.e) if f["valid"] else None
        if r is not None and "R_err" in r:
            o.update(r)
            o["R_err_hesse"] = r["R_err"]
            Rb, Gb = [], []
            for Hb in HB:
                fb = fit_m(Hb, mc, kap, peakT, model, start=f["x"], hesse=False)
                rb = A.r_of(fb["N"], None, mc.e) if fb["valid"] else None
                if rb is not None:
                    Rb.append(rb["R"])
            o["boot_rms"] = float(np.std(Rb)) if len(Rb) > 5 else None
            o["boot_ok"] = len(Rb)
            if o["boot_rms"]:
                o["R_err"] = max(o["R_err_hesse"], o["boot_rms"])
            o["bias"] = o["R"] - Rt
            o["N_fit_over_truth"] = [f["N"][j] / Tn[j] if Tn[j] else None for j in range(3)]
            o["pass"] = bool(f["valid"] and not f["at_limit"] and abs(o["bias"]) <= 2.5 * o["R_err"] and abs(o["bias"]) <= 0.04 and o["R_err"] <= 0.05)
        else:
            o["pass"] = False
        res[model] = o
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nboot", type=int, default=100)
    ap.add_argument("--tags", default="s1,s2,s3")
    ap.add_argument("--nproc", type=int, default=3)
    a = ap.parse_args()
    t0 = time.time()
    tasks = []
    for s in a.tags.split(","):
        tag = f"D0.tr.H+e04.mc.d5.{s}"
        for b in A.BINS:
            for wp in ((0.7, 0.5) if s == "s1" else (0.7,)):
                tasks.append((tag, b, wp, a.nboot if wp == 0.7 else 0))
    out = []
    with Pool(a.nproc, initializer=init) as p:
        for r in p.imap_unordered(task, tasks):
            out.append(r)
            print(f"[{time.time() - t0:4.0f}s] {r['tag']} {r['bin']} WP{int(r['wp'] * 100)}: R_prawda {r['R_truth']:.4f}; kap b/n/t {r['kap'][0]:.2f}/{r['kap'][1]:.2f}/{r['kap'][2]:.2f}; " +
                  "; ".join(f"{m}: R {r[m].get('R', float('nan')):.4f} ± {r[m].get('R_err', float('nan')):.4f} (H {r[m].get('R_err_hesse', float('nan')):.4f}) p {r[m]['p']:.2g} G {r[m].get('G_fit') or float('nan'):.3f} {'OK' if r[m]['pass'] else 'NIE'}"
                            for m in ("M0", "M8n", "M8")), flush=True)
    out.sort(key=lambda r: (r["tag"], [A.bname(b) for b in A.BINS].index(r["bin"]), -r["wp"]))
    json.dump({"res": out, "nboot": a.nboot}, open(os.path.join(HERE, "..", "out", "R_m8.json"), "w"), indent=1, default=float)
    print(f"zapisano, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
