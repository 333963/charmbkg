#!/usr/bin/env python3
"""FITBIAS diagnoza W — to samo co diag_T, ale na histogramach WPSCAN (rzeczywiste drzewo H, WP70; prawda tylko: sygnał / odbicie / reszta).
Wyjście: out/diag_W.json.  KRYTERIA §1.2."""
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fb_fit as F  # noqa: E402

WP = os.path.join(HERE, "..", "..", "WPSCAN_2026-11-19", "runs", "lxplus", "wps_hist.npz")
TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
Z = INFO = None


def _init():
    global Z, INFO
    Z = np.load(WP, allow_pickle=False)
    INFO = json.loads(str(Z["meta"]))


def sel(h, jlo, jhi=None):
    """przedziały score: ε_lo ≥ ... > ε_hi  ⇔ indeksy jlo+1 .. jhi (jhi=None → do końca)."""
    return h[..., jlo + 1:(None if jhi is None else jhi + 1), :].sum(-2)


def catdef(name):
    e = INFO["eps_desc"].index
    return {"WP90": (e(0.9), None), "WP70": (e(0.7), None), "A50": (e(0.5), None), "A40": (e(0.4), None), "A35": (e(0.35), None),
            "B50": (e(0.7), e(0.5)), "B40": (e(0.7), e(0.4)), "B35": (e(0.7), e(0.35)), "L90-70": (e(0.9), e(0.7))}[name]


def build(sp, bn, name):
    c = F.SP[sp]
    lo, hi = c["mr"]
    key = f"{sp}|{TAG[sp]}|{bn}"
    j = catdef(name)
    sw = sel(Z[key + "|fd|sig_w"], *j)
    o = {"sigT": F.tmpl1(sw, lo), "refT": None, "frefl_fd": 0.0, "frefl_mb": 0.0, "e": float(sw.sum() / Z[key + "|fd|sig_w"].sum())}
    comp = {"all": sel(Z[key + "|mb|all"], *j).astype(float), "Sig": sel(Z[key + "|mb|sig"], *j).astype(float)}
    comp["Refl"] = sel(Z[key + "|mb|refl"], *j).astype(float) if c["refl"] else np.zeros_like(comp["Sig"])
    comp["rest"] = comp["all"] - comp["Sig"] - comp["Refl"]
    if c["refl"]:
        rw = sel(Z[key + "|fd|refl_w"], *j)
        o["refT"], o["frefl_fd"] = F.tmpl1(rw, lo), float(rw.sum() / sw.sum())
        o["frefl_mb"] = float(comp["Refl"].sum() / comp["Sig"].sum())
        o["refT_mb"] = F.tmpl1(comp["Refl"], lo)
    o["comp"] = comp
    o["N_true"] = float(comp["Sig"].sum())
    return o


def brief(f, ntrue, i=0):
    return {"N": f["N"][i], "err": (f.get("err") or [None] * 9)[i], "b": f["N"][i] / ntrue - 1, "valid": f["valid"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"],
            "at_limit": f["at_limit"]}


def task(t):
    sp, bn = t
    c = F.SP[sp]
    lo, hi = c["mr"]
    flo, fhi = c["fit"]
    edges = np.arange(flo, fhi + F.BW / 2, F.BW)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    wm = np.abs(ctr - c["m0"]) < c["win"]
    res = {"sp": sp, "bin": bn, "single": {}, "joint": {}}
    B = {n: build(sp, bn, n) for n in ("WP90", "WP70", "A50", "A40", "A35", "B50", "B40", "B35", "L90-70")}
    rb = lambda h, r=None: F.rebin1(h, lo, *(r or c["fit"]))  # noqa: E731
    for n, b in B.items():
        nt = b["N_true"]
        H = {g: rb(h) for g, h in b["comp"].items()}
        o = {"N_true": nt, "frefl_fd": b["frefl_fd"], "frefl_mb": b["frefl_mb"], "eps": b["e"],
             "win": {g: float(h[wm].sum()) for g, h in H.items()}, "SB_win": float(H["Sig"][wm].sum() / max(H["all"][wm].sum() - H["Sig"][wm].sum(), 1e-9))}
        T_fd = [[(b["refT"], b["frefl_fd"])]] if c["refl"] else [[]]
        T_mb = [[(b["refT"], b["frefl_mb"])]] if c["refl"] else [[]]
        T_mbs = [[(b["refT_mb"], b["frefl_mb"])]] if c["refl"] else [[]]
        ft = lambda h, shape, tied, nbkg=None, **kw: F.fit(h[None], [b["sigT"]], tied, edges, c["m0"], c["nbkg"] if nbkg is None else nbkg, shape=shape, win=c["win"], **kw)  # noqa: E731
        for shape in ("fix", "free"):
            f = ft(H["Sig"], shape, [[]])
            o[f"Sig|{shape}"] = brief(f, nt)
            o[f"Sig|{shape}"]["b_err"] = math.sqrt(max(f["err"][0] ** 2 - nt, 0)) / nt if "err" in f else None
            if c["refl"]:
                o[f"Sig+Refl|{shape}|frefl_fd"] = brief(ft(H["Sig"] + H["Refl"], shape, T_fd), nt)
                o[f"Sig+Refl|{shape}|frefl_mb"] = brief(ft(H["Sig"] + H["Refl"], shape, T_mb), nt)
                o[f"all|{shape}|frefl_mb"] = brief(ft(H["all"], shape, T_mb), nt)
                o[f"all|{shape}|refl_mb_shape"] = brief(ft(H["all"], shape, T_mbs), nt)
                o[f"Sig+rest|{shape}"] = brief(ft(H["Sig"] + H["rest"], shape, [[]]), nt)
            o[f"all|{shape}|frefl_fd"] = brief(ft(H["all"], shape, T_fd), nt)
            for k in (c["nbkg"] + 1, c["nbkg"] + 2):
                o[f"all|{shape}|frefl_fd|exp{k}"] = brief(ft(H["all"], shape, T_fd, nbkg=k), nt)
        for k in range(c["nbkg"] - 1, c["nbkg"] + 3):
            f0 = ft(H["rest"], "fix", [[]], nbkg=k, fixN0=True)
            f1 = ft(H["rest"], "fix", [[]], nbkg=k)
            o[f"rest|exp{k}"] = {"chi2_0": f0["chi2"], "ndf_0": f0["ndf"], "p_0": F.pval(f0["chi2"], f0["ndf"]), "N_spur": f1["N"][0], "err": (f1.get("err") or [None])[0],
                                 "rel": f1["N"][0] / nt, "rel_err": f1["err"][0] / nt if "err" in f1 else None}
        res["single"][n] = o
    # fit łączny A + B
    for e in ("50", "40", "35"):
        a, b = B["A" + e], B["B" + e]
        nA, nB = a["N_true"], b["N_true"]
        f_true = nA / (nA + nB)
        f_mc = a["e"] / (a["e"] + b["e"])
        o = {"f_true": f_true, "f_mc": f_mc, "NA_true": nA, "NB_true": nB}
        sigT = [a["sigT"], b["sigT"]]

        def jf(HA, HB, shape, tied):
            f = F.fit(np.stack([HA, HB]), sigT, tied, edges, c["m0"], c["nbkg"], shape=shape, win=c["win"])
            NA, NB = f["N"]
            r = {"NA": NA, "NB": NB, "bA": NA / nA - 1, "bB": NB / nB - 1, "f": NA / (NA + NB), "SF": NA / (NA + NB) / f_mc, "SF_true": f_true / f_mc, "valid": f["valid"], "sc": f["sc"],
                 "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"]}
            if "cov_N" in f:
                C = np.array(f["cov_N"]); E = NA + NB
                g = np.array([NB / E ** 2, -NA / E ** 2])
                r["f_err"] = math.sqrt(max(g @ C @ g, 0)); r["SF_err"] = r["f_err"] / f_mc
                r["NA_err"], r["NB_err"] = math.sqrt(C[0, 0]), math.sqrt(C[1, 1])
            return r
        HA = {g: rb(h) for g, h in a["comp"].items()}
        HB = {g: rb(h) for g, h in b["comp"].items()}
        Tfd = [[(a["refT"], a["frefl_fd"])], [(b["refT"], b["frefl_fd"])]] if c["refl"] else [[], []]
        Tmb = [[(a["refT"], a["frefl_mb"])], [(b["refT"], b["frefl_mb"])]] if c["refl"] else [[], []]
        for shape in ("fix", "free"):
            o[f"all|{shape}|frefl_fd"] = jf(HA["all"], HB["all"], shape, Tfd)
            o[f"Sig|{shape}"] = jf(HA["Sig"], HB["Sig"], shape, [[], []])
            if c["refl"]:
                o[f"all|{shape}|frefl_mb"] = jf(HA["all"], HB["all"], shape, Tmb)
                o[f"Sig+Refl|{shape}|frefl_fd"] = jf(HA["Sig"] + HA["Refl"], HB["Sig"] + HB["Refl"], shape, Tfd)
                o[f"Sig+rest|{shape}"] = jf(HA["Sig"] + HA["rest"], HB["Sig"] + HB["rest"], shape, [[], []])
        res["joint"][e] = o
    return res


def main():
    t0 = time.time()
    tasks = [(sp, bn) for sp in TAG for bn in F.BINS]
    with Pool(6, initializer=_init) as p:
        out = p.map(task, tasks)
    json.dump({"res": out}, open(os.path.join(HERE, "..", "out", "diag_W.json"), "w"))
    print(f"zapisano diag_W.json, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
