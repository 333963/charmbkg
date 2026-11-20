#!/usr/bin/env python3
"""WPSCAN — fity piku wzdłuż skanu punktu pracy (kryteria: KRYTERIA.md). Wejście: wps_hist.npz (Q_lxplus/wpscan/wps_hist.py).
Model = dn4_ana.fit_cats z jedną kategorią (bez zmian). Na punkt (gatunek, model, bin, próbka|klasa, ε): fit nominalny, warianty tła
(rząd ±1, wielomian bez exp, zakres alternatywny), repliki bootstrapu (fit nominalny). Domknięcie: pseudo-dane mb_nd.
  wps_fit.py fit --hist wps_hist.npz --out wps_fit.json [--nproc 6] [--nrep 30] [--classes MB427929,all]
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

Z = None
INFO = None
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
    """hist 1 MeV od lo → biny 5 MeV w [flo, fhi]."""
    k = int(round(A.BW / INFO["mbw"]))
    i0, i1 = int(round((flo - lo) / INFO["mbw"])), int(round((fhi - lo) / INFO["mbw"]))
    x = h[..., i0:i1]
    return x.reshape(x.shape[:-1] + ((i1 - i0) // k, k)).sum(-1).astype(np.float64)


def model_mu(f, sigT, refT, frefl, edges, m0, nbkg, kind="exp"):
    v = np.array(f["x"])
    ctr = 0.5 * (edges[1:] + edges[:-1])
    x = 2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1
    s = sigT.bins(edges, m0, v[0], v[1])
    if refT is not None and frefl > 0:
        s = s + frefl * refT.bins(edges, m0, v[0], 1.0)
    sh = A.bshape(x, v[4:4 + nbkg], kind)
    return np.maximum(v[2] * s + v[3] * sh / sh.sum(), 1e-9)


def one(H, sigT, refT, frefl, edges, m0, nbkg, win, kind="exp", start=None, detail=False):
    f = A.fit_cats(H[None, :], [sigT], [refT], [frefl], edges, m0, nbkg, kind=kind, start=start, win=win)
    o = {"N": f["N"][0], "valid": f["valid"], "at_limit": f["at_limit"], "chi2": f["chi2"], "ndf": f["ndf"], "dmu": f["dmu"], "sc": f["sc"],
         "err": math.sqrt(max(f["cov_N"][0][0], 0)) if "cov_N" in f else None, "x": f["x"]}
    if detail:
        from scipy.stats import chi2 as C2
        mu = model_mu(f, sigT, refT, frefl, edges, m0, nbkg, kind)
        pull = np.where(mu >= 10, (H - mu) / np.sqrt(mu), 0.0)      # pull Pearsona tylko dla binów z μ ≥ 10 (KRYTERIA §2, doprecyzowanie)
        ctr = 0.5 * (edges[1:] + edges[:-1])
        o.update({"p": float(C2.sf(f["chi2"], f["ndf"])), "maxpull": float(np.abs(pull).max()), "m_maxpull": float(ctr[np.abs(pull).argmax()]),
                  "n_lowmu": int((mu < 10).sum()), "bkg_win": f["bkg_win"][0], "NB": f["NB"][0], "n_rows": float(H.sum()),
                  "pull_peak_chi2": float((pull[np.abs(ctr - m0) < win] ** 2).sum()), "n_peak_bins": int((np.abs(ctr - m0) < win).sum())})
    return o


def task(t):
    sp, tag, bn, src, j = t
    c = A.SP[sp]
    lo, hi = INFO["sp"][sp]["mr"]
    key = f"{sp}|{tag}|{bn}"
    eps = float(Z[key + "|eps"][j])
    cum = lambda h: h[..., j + 1:, :].sum(-2)  # noqa: E731
    sw, sn = cum(Z[key + "|fd|sig_w"]), cum(Z[key + "|fd|sig_n"])
    res = {"sp": sp, "tag": tag, "bin": bn, "src": src, "j": j, "eps_nom": INFO["eps_desc"][j], "eps": eps, "thr": float(Z[key + "|thr"][j]),
           "n_mc_pass": int(sn.sum()), "n_mc_tot": int(Z[key + "|nsig_tot"][0])}
    if sn.sum() < 300:
        res["skip"] = "za mało sygnału MC"; return res
    sigT = tmpl(sw, lo, hi)
    refT, frefl = None, 0.0
    if c["refl"]:
        rw, rn = cum(Z[key + "|fd|refl_w"]), cum(Z[key + "|fd|refl_n"])
        if rn.sum() >= 300:
            refT, frefl = tmpl(rw, lo, hi), float(rw.sum() / sw.sum())
    res["frefl"] = frefl
    flo, fhi = c["fit"]
    edges = np.arange(flo, fhi + A.BW / 2, A.BW)
    if src == "mb|all":
        Hall = rebin(cum(Z[f"{key}|mb|all"]), lo, flo, fhi)
        f = one(Hall, sigT, refT, frefl, edges, c["m0"], c["nbkg"], c["win"], detail=True)
        f.pop("x")
        res["nom"] = f
        res["N_truth"] = int(cum(Z[f"{key}|mb|sig"]).sum())
        res["N_truth_refl"] = int(cum(Z[f"{key}|mb|refl"]).sum()) if c["refl"] else 0
        return res
    D = cum(Z[f"{key}|{src}"])                    # [nrep + 1, nm]
    Hs = rebin(D, lo, flo, fhi)
    f = one(Hs[0], sigT, refT, frefl, edges, c["m0"], c["nbkg"], c["win"], detail=True)
    start = f.pop("x")
    res["nom"] = f
    var = {}
    for nm_, kw in ((f"exp{c['nbkg'] - 1}", dict(nbkg=c["nbkg"] - 1)), (f"exp{c['nbkg'] + 1}", dict(nbkg=c["nbkg"] + 1)), (f"pol{c['nbkg']}", dict(nbkg=c["nbkg"], kind="pol"))):
        v = one(Hs[0], sigT, refT, frefl, edges, c["m0"], kw["nbkg"], c["win"], kind=kw.get("kind", "exp"), detail=True)
        v.pop("x"); var[nm_] = v
    alo, ahi = max(c["alt"][0], lo), min(c["alt"][1], hi)
    e2 = np.arange(alo, ahi + A.BW / 2, A.BW)
    v = one(rebin(D[0], lo, alo, ahi), sigT, refT, frefl, e2, c["m0"], c["nbkg"], c["win"], detail=True)
    v.pop("x"); var["range"] = v
    res["var"] = var
    bN, bE = [], []
    for r in range(1, min(NREP, Hs.shape[0] - 1) + 1):
        fb = one(Hs[r], sigT, refT, frefl, edges, c["m0"], c["nbkg"], c["win"], start=start)
        ok = fb["valid"] and fb["err"] is not None
        bN.append(fb["N"] if ok else None); bE.append(fb["err"] if ok else None)
    res["boot_N"], res["boot_err"] = bN, bE
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["fit"])
    ap.add_argument("--hist", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--nproc", type=int, default=6); ap.add_argument("--nrep", type=int, default=30)
    ap.add_argument("--classes", default="427929|MB427929,428770|all,mb|all")
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    _init(a.hist, a.nrep)
    t0 = time.time()
    tasks = []
    for sp, tags in INFO["tags"].items():
        for tag in tags:
            for b in INFO["bins"]:
                bn = A.bname(b)
                for src in a.classes.split(","):
                    if f"{sp}|{tag}|{bn}|{src}" not in Z.files:
                        continue
                    for j in range(len(INFO["eps_desc"])):
                        tasks.append((sp, tag, bn, src, j))
    if a.only:
        tasks = [t for t in tasks if a.only in "|".join(map(str, t))]
    print(f"zadań: {len(tasks)}", flush=True)
    out = []
    with Pool(a.nproc, initializer=_init, initargs=(a.hist, a.nrep)) as p:
        for i, r in enumerate(p.imap_unordered(task, tasks, chunksize=2)):
            out.append(r)
            if (i + 1) % 50 == 0:
                print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)}", flush=True)
    out.sort(key=lambda r: (r["sp"], r["tag"], r["bin"], r["src"], r["j"]))
    json.dump({"info": INFO, "res": out}, open(a.out, "w"))
    print(f"zapisano {a.out}, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
