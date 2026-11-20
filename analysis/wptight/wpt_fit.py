#!/usr/bin/env python3
"""WPTIGHT — fit pass/fail (A = ciasny, B = WP70 i nie ciasny) inkluzywnego piku; kryteria: KRYTERIA.md.
Wejście: WPSCAN wps_hist.npz (masa 1 MeV × 19 przedziałów score, 30 replik bootstrapu po zdarzeniach). Model: dn4_ana.fit_cats (bez zmian).
  wpt_fit.py mc   --hist wps_hist.npz --out wpt_mc.json      domknięcie na mb_nd (pseudo-dane), repliki Poissona po binach
  wpt_fit.py data --hist wps_hist.npz --out wpt_data.json    dane: 427929|MB427929, 428770|all, 428770|J12only, 428770|notJ12only
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "Q_lxplus", "dnet4"))
import dn4_ana as A  # noqa: E402

TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
EPS_T = [0.5, 0.4, 0.35]
Z = INFO = None
NREP = 30


def _init(path, nrep):
    global Z, INFO, NREP
    Z = np.load(path, allow_pickle=False)
    INFO = json.loads(str(Z["meta"]))
    NREP = nrep


def tmpl(h, lo, hi):
    ctr = lo + (np.arange(len(h)) + 0.5) * INFO["mbw"]
    return A.Tmpl(ctr, h, lo, hi)


def rebin(h, lo, flo, fhi):
    k = int(round(A.BW / INFO["mbw"]))
    i0, i1 = int(round((flo - lo) / INFO["mbw"])), int(round((fhi - lo) / INFO["mbw"]))
    x = h[..., i0:i1]
    return x.reshape(x.shape[:-1] + ((i1 - i0) // k, k)).sum(-1).astype(np.float64)


def cats(h, j70, jt):
    """h [..., 19, nm] → [..., 2, nm]: A = przedziały > jt, B = przedziały j70+1..jt."""
    return np.stack([h[..., jt + 1:, :].sum(-2), h[..., j70 + 1:jt + 1, :].sum(-2)], axis=-2)


def pval(chi2, ndf):
    from scipy.stats import chi2 as C2
    return float(C2.sf(chi2, ndf))


def lowmu_frac(f, H):
    return float((H < 10).mean())


def jfit(H2, mc, edges, c, nbkg=None, kind="exp", start=None, common=False):
    """Fit łączny A+B. → dict z N_A, N_B, f, E, błędy HESSE, korelacja."""
    sig = [mc["allT"], mc["allT"]] if common else mc["sigT"]
    f = A.fit_cats(H2, sig, mc["refT"], mc["frefl"], edges, c["m0"], c["nbkg"] if nbkg is None else nbkg, kind=kind, start=start, win=c["win"])
    NA, NB = f["N"]
    o = {"NA": NA, "NB": NB, "E": NA + NB, "f": NA / (NA + NB) if NA + NB > 0 else None, "valid": f["valid"], "at_limit": f["at_limit"],
         "chi2": f["chi2"], "ndf": f["ndf"], "dmu": f["dmu"], "sc": f["sc"], "x": f["x"], "bkg_win": f["bkg_win"]}
    if "cov_N" in f and NA + NB > 0 and NA > 0:
        C = np.array(f["cov_N"])
        E = NA + NB
        o["NA_err"] = math.sqrt(max(C[0, 0], 0)); o["NB_err"] = math.sqrt(max(C[1, 1], 0))
        o["E_err"] = math.sqrt(max(C.sum(), 0))
        g = np.array([NB / E ** 2, -NA / E ** 2])
        o["f_err"] = math.sqrt(max(g @ C @ g, 0))
        o["cov_AB"] = float(C[0, 1])
        cf = (C[0] @ g)                                   # cov(N_A, f)
        o["corr_NA_f"] = float(cf / (o["NA_err"] * o["f_err"])) if o["NA_err"] > 0 and o["f_err"] > 0 else None
    return o


def sfit(H, sigT, refT, frefl, edges, c, nbkg=None, kind="exp", start=None):
    f = A.fit_cats(H[None, :], [sigT], [refT], [frefl], edges, c["m0"], c["nbkg"] if nbkg is None else nbkg, kind=kind, start=start, win=c["win"])
    return {"N": f["N"][0], "err": math.sqrt(max(f["cov_N"][0][0], 0)) if "cov_N" in f else None, "valid": f["valid"], "at_limit": f["at_limit"],
            "chi2": f["chi2"], "ndf": f["ndf"], "dmu": f["dmu"], "sc": f["sc"], "x": f["x"], "bkg_win": f["bkg_win"][0]}


def build_mc(sp, bn, jt):
    c = A.SP[sp]
    lo, hi = INFO["sp"][sp]["mr"]
    key = f"{sp}|{TAG[sp]}|{bn}"
    j70 = INFO["eps_desc"].index(0.7)
    sw, sn = Z[key + "|fd|sig_w"], Z[key + "|fd|sig_n"]
    W = float(sw.sum())
    cw, cn = cats(sw, j70, jt), cats(sn, j70, jt)
    eA, eB = float(cw[0].sum() / W), float(cw[1].sum() / W)
    nA, nB = int(cn[0].sum()), int(cn[1].sum())
    fmc = eA / (eA + eB)
    mc = {"eA": eA, "eB": eB, "f_mc": fmc, "f_mc_err": math.sqrt(fmc * (1 - fmc) / max(nA + nB, 1)), "nA": nA, "nB": nB,
          "sigT": [tmpl(cw[0], lo, hi), tmpl(cw[1], lo, hi)], "allT": tmpl(cw.sum(0), lo, hi), "refT": [None, None], "frefl": [0.0, 0.0],
          "refAll": None, "freflAll": 0.0, "refA": None}
    if c["refl"]:
        rw = cats(Z[key + "|fd|refl_w"], j70, jt)
        mc["refT"] = [tmpl(rw[0], lo, hi), tmpl(rw[1], lo, hi)]
        mc["frefl"] = [float(rw[0].sum() / cw[0].sum()), float(rw[1].sum() / cw[1].sum())]
        mc["refAll"], mc["freflAll"] = tmpl(rw.sum(0), lo, hi), float(rw.sum() / cw.sum())
    return mc, j70, key, c, lo, hi


def scal(o, keys):
    return {k: o.get(k) for k in keys}


JK = ("NA", "NB", "E", "f", "NA_err", "NB_err", "E_err", "f_err", "corr_NA_f", "valid", "at_limit", "chi2", "ndf", "dmu", "sc")
SK = ("N", "err", "valid", "at_limit", "chi2", "ndf", "dmu", "sc", "bkg_win")


def task(t):
    sp, bn, src, eps_t = t
    jt = INFO["eps_desc"].index(eps_t)
    mc, j70, key, c, lo, hi = build_mc(sp, bn, jt)
    flo, fhi = c["fit"]
    edges = np.arange(flo, fhi + A.BW / 2, A.BW)
    res = {"sp": sp, "bin": bn, "src": src, "eps_t": eps_t, "eps_t_mc": float(Z[key + "|eps"][jt]), "eps70_mc": float(Z[key + "|eps"][j70]),
           **scal(mc, ("eA", "eB", "f_mc", "f_mc_err", "nA", "nB", "frefl", "freflAll"))}
    is_mc = src == "mb|all"
    D = Z[f"{key}|{src}"]
    if is_mc:                                            # pseudo-dane mb_nd; repliki Poissona po binach (5 MeV)
        H0 = rebin(cats(D, j70, jt), lo, flo, fhi)
        rng = np.random.default_rng(1000 * jt + hash((sp, bn)) % 1000)
        Hs = np.concatenate([H0[None], rng.poisson(H0, size=(NREP,) + H0.shape).astype(np.float64)])
        ts = cats(Z[f"{key}|mb|sig"], j70, jt)
        res["truth"] = {"NA": int(ts[0].sum()), "NB": int(ts[1].sum())}
        if c["refl"]:
            tr = cats(Z[f"{key}|mb|refl"], j70, jt)
            res["truth"]["refl_A"], res["truth"]["refl_B"] = int(tr[0].sum()), int(tr[1].sum())
    else:
        Hs = rebin(cats(D, j70, jt), lo, flo, fhi)        # [nrep+1, 2, nb]
    res["n_rows"] = [float(Hs[0][0].sum()), float(Hs[0][1].sum())]
    res["lowmu_frac"] = float((Hs[0] < 10).mean())
    # --- fit łączny
    f = jfit(Hs[0], mc, edges, c)
    start = f.pop("x")
    f["p"] = pval(f["chi2"], f["ndf"])
    res["joint"] = f
    var = {}
    for nm, kw in ((f"exp{c['nbkg'] - 1}", dict(nbkg=c["nbkg"] - 1)), (f"exp{c['nbkg'] + 1}", dict(nbkg=c["nbkg"] + 1)),
                   (f"pol{c['nbkg']}", dict(kind="pol")), ("common_tmpl", dict(common=True))):
        v = jfit(Hs[0], mc, edges, c, **kw); v.pop("x"); var[nm] = scal(v, JK)
    alo, ahi = max(c["alt"][0], lo), min(c["alt"][1], hi)
    e2 = np.arange(alo, ahi + A.BW / 2, A.BW)
    Halt = rebin(cats(D if is_mc else D[0], j70, jt), lo, alo, ahi)
    v = jfit(Halt, mc, e2, c); v.pop("x"); var["range"] = scal(v, JK)
    res["joint_var"] = var
    # --- fit jednej kategorii: WP70 (A+B) i ciasny (A)
    H70 = Hs[:, 0] + Hs[:, 1]
    s70 = sfit(H70[0], mc["allT"], mc["refAll"], mc["freflAll"], edges, c)
    st70 = s70.pop("x"); s70["p"] = pval(s70["chi2"], s70["ndf"])
    res["wp70"] = s70
    v70 = {}
    for nm, kw in ((f"exp{c['nbkg'] - 1}", dict(nbkg=c["nbkg"] - 1)), (f"exp{c['nbkg'] + 1}", dict(nbkg=c["nbkg"] + 1)), (f"pol{c['nbkg']}", dict(kind="pol"))):
        v = sfit(H70[0], mc["allT"], mc["refAll"], mc["freflAll"], edges, c, **kw); v.pop("x"); v70[nm] = scal(v, SK)
    v = sfit(Halt.sum(0), mc["allT"], mc["refAll"], mc["freflAll"], e2, c); v.pop("x"); v70["range"] = scal(v, SK)
    res["wp70_var"] = v70
    sA = sfit(Hs[0, 0], mc["sigT"][0], mc["refT"][0], mc["frefl"][0], edges, c)
    stA = sA.pop("x"); sA["p"] = pval(sA["chi2"], sA["ndf"])
    res["tight"] = sA
    # --- repliki
    bj, b70, bA = [], [], []
    for r in range(1, min(NREP, Hs.shape[0] - 1) + 1):
        fb = jfit(Hs[r], mc, edges, c, start=start)
        bj.append(scal(fb, ("NA", "NB", "E", "f", "E_err", "f_err", "NA_err")) if fb["valid"] and "E_err" in fb else None)
        fb = sfit(H70[r], mc["allT"], mc["refAll"], mc["freflAll"], edges, c, start=st70)
        b70.append(scal(fb, ("N", "err")) if fb["valid"] and fb["err"] is not None else None)
        fb = sfit(Hs[r, 0], mc["sigT"][0], mc["refT"][0], mc["frefl"][0], edges, c, start=stA)
        bA.append(scal(fb, ("N", "err")) if fb["valid"] and fb["err"] is not None else None)
    res["boot_joint"], res["boot_wp70"], res["boot_tight"] = bj, b70, bA
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["mc", "data"])
    ap.add_argument("--hist", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--nproc", type=int, default=6); ap.add_argument("--nrep", type=int, default=30)
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    _init(a.hist, a.nrep)
    srcs = ["mb|all"] if a.cmd == "mc" else ["427929|MB427929", "428770|all", "428770|J12only", "428770|notJ12only"]
    tasks = [(sp, A.bname(b), src, e) for sp in TAG for b in INFO["bins"] for src in srcs for e in EPS_T]
    if a.only:
        tasks = [t for t in tasks if a.only in "|".join(map(str, t))]
    print(f"zadań: {len(tasks)}", flush=True)
    t0 = time.time()
    out = []
    with Pool(a.nproc, initializer=_init, initargs=(a.hist, a.nrep)) as p:
        for i, r in enumerate(p.imap_unordered(task, tasks)):
            out.append(r)
            if (i + 1) % 6 == 0:
                print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)}", flush=True)
    out.sort(key=lambda r: (r["sp"], r["bin"], r["src"], -r["eps_t"]))
    json.dump({"info": INFO, "res": out}, open(a.out, "w"))
    print(f"zapisano {a.out}, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
