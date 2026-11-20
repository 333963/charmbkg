#!/usr/bin/env python3
"""FITBIAS diagnoza T — rozkład obciążenia fitu na składowe prawdy (MC_templates: mb_nd, fd_nd; 24 kategorie prawdy × 4 przedziały score).
KRYTERIA §1.2 pkt 1–5. Wyjście: out/diag_T.json.  Użycie: fb_diag_T.py [--ntoy 200] [--nproc 7]"""
import argparse
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

TD = os.path.join(HERE, "..", "..", "MC_templates", "runs", "lxplus")
GROUPS = {"Sig": [0, 1], "Refl": [2, 3], "CM": [4, 5, 6, 7, 8, 9, 10], "DMM": [11, 12], "PC": [13, 14, 15, 16, 17, 18, 19], "Comb": [20, 21, 22, 23]}
CATS = {"WP95": [1, 2, 3], "WP80": [2, 3], "WP50": [3], "L95-80": [1], "L80-50": [2]}
M = FD = None


def _init():
    global M, FD
    M = np.load(os.path.join(TD, "mctmpl_mb_nd.npz"))
    FD = np.load(os.path.join(TD, "mctmpl_fd_nd.npz"))


def comp(Z, sp, ib, sbs):
    """→ dict grupa → hist [nm] (5 MeV, pełny zakres ekstrakcji)."""
    h = Z["d0_m" if sp == "D0" else "dp_m"][ib]            # [sb, K, cell, (veto), m]
    h = h.sum(2) if sp == "D0" else h.sum(2)[:, :, 1]
    h = h[sbs].sum(0).astype(float)                        # [K, m]
    o = {g: h[ix].sum(0) for g, ix in GROUPS.items()}
    o["fine"] = h
    return o


def cut(h, sp, rng=None):
    c = F.SP[sp]
    lo, hi = c["mr"]
    flo, fhi = rng or c["fit"]
    i0, i1 = int(round((flo - lo) / F.BW)), int(round((fhi - lo) / F.BW))
    return h[..., i0:i1]


def rel_err(f, ntrue):
    """błąd N_fit/N_prawda po odjęciu fluktuacji Poissona samego sygnału (te same zdarzenia w liczniku i mianowniku)."""
    if "err" not in f:
        return None
    return math.sqrt(max(f["err"][0] ** 2 - ntrue, 0.0)) / ntrue


def task(t):
    sp, ib, cat, ntoy = t
    c = F.SP[sp]
    lo, hi = c["mr"]
    full_edges = np.arange(lo, hi + F.BW / 2, F.BW)
    mb, fd = comp(M, sp, ib, CATS[cat]), comp(FD, sp, ib, CATS[cat])
    ntrue = float(mb["Sig"].sum())
    sT = F.TmplE(full_edges, fd["Sig"])
    frefl_fd = float(fd["Refl"].sum() / fd["Sig"].sum()) if c["refl"] else 0.0
    frefl_mb = float(mb["Refl"].sum() / mb["Sig"].sum()) if c["refl"] else 0.0
    rT = F.TmplE(full_edges, fd["Refl"]) if c["refl"] else None
    tied = [[(rT, frefl_fd)]] if c["refl"] else [[]]
    flo, fhi = c["fit"]
    edges = np.arange(flo, fhi + F.BW / 2, F.BW)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    wm = np.abs(ctr - c["m0"]) < c["win"]
    res = {"sp": sp, "bin": F.BINS[ib], "cat": cat, "N_true": ntrue, "frefl_fd": frefl_fd, "frefl_mb": frefl_mb,
           "n_fd_sig": float(fd["Sig"].sum()),
           "win": {g: float(cut(mb[g], sp)[wm].sum()) for g in GROUPS}, "tot": {g: float(cut(mb[g], sp).sum()) for g in GROUPS},
           "kappa_fd": {g: float(fd[g].sum() / fd["Sig"].sum()) for g in GROUPS}, "kappa_mb": {g: float(mb[g].sum() / mb["Sig"].sum()) for g in GROUPS}}
    res["SB_win"] = res["win"]["Sig"] / max(sum(v for g, v in res["win"].items() if g != "Sig"), 1e-9)

    def fit1(H, shape="free", nbkg=None, tied_=None, rng=None, **kw):
        e = edges if rng is None else np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
        f = F.fit(H[None], [sT], tied if tied_ is None else tied_, e, c["m0"], c["nbkg"] if nbkg is None else nbkg, shape=shape, win=c["win"], **kw)
        f.pop("mu", None)
        return f

    def brief(f):
        return {"N": f["N"][0], "err": (f.get("err") or [None])[0], "b": f["N"][0] / ntrue - 1, "valid": f["valid"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"],
                "at_limit": f["at_limit"]}
    # (1)+(2)+(3) dokładanie składowych — dwie kolejności, kształt stały / wolny
    H = {g: cut(mb[g], sp) for g in GROUPS}
    seqs = {"I": ["Sig", "Refl", "DMM", "CM", "PC", "Comb"], "II": ["Sig", "Comb", "PC", "CM", "DMM", "Refl"]}
    res["seq"] = {}
    for sn, order in seqs.items():
        for shape in ("fix", "free"):
            acc = np.zeros_like(H["Sig"])
            rows = []
            has_refl = False
            for g in order:
                if g == "Refl" and not c["refl"]:
                    continue
                acc = acc + H[g]
                has_refl = has_refl or g == "Refl"
                f = fit1(acc, shape, tied_=None if has_refl else [[]])     # odbicie w modelu dopiero, gdy jest w pseudo-danych
                r = brief(f); r["add"] = g
                if g == "Sig":
                    r["b_err"] = rel_err(f, ntrue)
                rows.append(r)
            res["seq"][f"{sn}|{shape}"] = rows
    # odbicie: wariant z f_refl z prawdy mb (kształt fd) — pełne pseudo-dane
    allH = sum(H.values())
    if c["refl"]:
        for shape in ("fix", "free"):
            res[f"all_freflmb|{shape}"] = brief(fit1(allH, shape, tied_=[[(rT, frefl_mb)]]))
            res[f"all_norefl|{shape}"] = brief(fit1(allH, shape, tied_=[[]]))
    # (3) sygnał pozorny każdej składowej tła osobno (sygnał stały, bez odbicia w modelu)
    res["spur"] = {}
    fine = cut(mb["fine"], sp)
    names = ["SigC", "SigB", "ReflSwap", "ReflOther", "CM_411", "CM_421", "CM_431", "CM_413", "CM_423", "CM_4122", "CM_other", "DMM_C", "DMM_B",
             "PC_411", "PC_421", "PC_423", "PC_413", "PC_431", "PC_4122", "PC_other", "MatchedBL", "OtherNoB", "OtherB", "Unknown"]
    items = {g: H[g] for g in ("CM", "DMM", "PC", "Comb")}
    items.update({names[k]: fine[k] for k in range(4, 20) if fine[k].sum() >= 50})
    items["PC+Comb"] = H["PC"] + H["Comb"]
    items["bkg_all"] = H["CM"] + H["DMM"] + H["PC"] + H["Comb"]
    for g, h in items.items():
        f = fit1(h, "fix", tied_=[[]])
        res["spur"][g] = {"N": f["N"][0], "err": (f.get("err") or [None])[0], "rel": f["N"][0] / ntrue, "rel_err": (f["err"][0] / ntrue if "err" in f else None),
                          "n_tot": float(h.sum()), "n_win": float(h[wm].sum()), "chi2": f["chi2"], "ndf": f["ndf"], "valid": f["valid"]}
    # (4) gładkie tło: rząd i zakres; samo tło (N = 0) oraz sygnał pozorny
    res["smooth"] = {}
    for gname, hfull in (("PC+Comb", mb["PC"] + mb["Comb"]), ("bkg_all", mb["CM"] + mb["DMM"] + mb["PC"] + mb["Comb"])):
        for rname, rng in (("nom", c["fit"]), ("alt", c["alt"])):
            h = cut(hfull, sp, rng)
            for k in range(c["nbkg"] - 1, c["nbkg"] + 3):
                f0 = fit1(h, "fix", nbkg=k, tied_=[[]], rng=rng, fixN0=True)
                f1 = fit1(h, "fix", nbkg=k, tied_=[[]], rng=rng)
                res["smooth"][f"{gname}|{rname}|exp{k}"] = {"chi2_0": f0["chi2"], "ndf_0": f0["ndf"], "p_0": F.pval(f0["chi2"], f0["ndf"]), "N_spur": f1["N"][0],
                                                            "err": (f1.get("err") or [None])[0], "rel": f1["N"][0] / ntrue, "rel_err": (f1["err"][0] / ntrue if "err" in f1 else None),
                                                            "valid": f1["valid"]}
    # (4b) Asimov bez szumu: gładkie tło = dopasowanie exp(nbkg+2) do bkg_all; sygnał = N_true·(szablon fd + f_refl_fd·odbicie fd)
    hb = H["CM"] + H["DMM"] + H["PC"] + H["Comb"]
    fa = F.fit(hb[None], [sT], [[]], edges, c["m0"], c["nbkg"] + 2, shape="fix", win=c["win"], fixN0=True)
    smooth = fa["mu"][0]
    s_as = ntrue * (sT.bins(edges, c["m0"]) + (frefl_fd * rT.bins(edges, c["m0"]) if c["refl"] else 0.0))
    res["asimov"] = {}
    for shape in ("fix", "free"):
        res["asimov"][f"model_sig+smooth{c['nbkg'] + 2}|{shape}"] = brief(fit1(s_as + smooth, shape))          # tylko rząd tła
        res["asimov"][f"true_sig+smooth{c['nbkg'] + 2}|{shape}"] = brief(fit1(H["Sig"] + (H["Refl"] if c["refl"] else 0) + smooth, shape))  # + kształt sygnału/odbicia
        res["asimov"][f"all|{shape}"] = brief(fit1(allH, shape))
        res["asimov"][f"all|{shape}|exp{c['nbkg'] + 1}"] = brief(fit1(allH, shape, nbkg=c["nbkg"] + 1))
        res["asimov"][f"all|{shape}|exp{c['nbkg'] + 2}"] = brief(fit1(allH, shape, nbkg=c["nbkg"] + 2))
    # (5) obciążenie własne estymatora: zabawki z modelu nominalnego dopasowanego do pełnych pseudo-danych (prawda = model)
    f0 = F.fit(allH[None], [sT], tied, edges, c["m0"], c["nbkg"], shape="fix", win=c["win"])
    # generator: sygnał o N = N_true (kształt szablonu), tło z fitu o stałym kształcie przeskalowane tak, by suma = dane
    sh = F.bshape(2 * (ctr - edges[0]) / (edges[-1] - edges[0]) - 1, np.array(f0["x"][4:4 + c["nbkg"]]))
    sg = sT.bins(edges, c["m0"]) + (frefl_fd * rT.bins(edges, c["m0"]) if c["refl"] else 0.0)
    gen = ntrue * sg + max(allH.sum() - ntrue * sg.sum(), 1.0) * sh / sh.sum()
    rng_ = np.random.default_rng(12345 + 100 * ib + 7 * list(CATS).index(cat) + (1000 if sp == "D0" else 0))
    toys = {s: [] for s in ("fix", "free")}
    for _ in range(ntoy):
        Ht = rng_.poisson(gen).astype(float)
        for shape in ("fix", "free"):
            f = fit1(Ht, shape)
            toys[shape].append([f["N"][0], (f.get("err") or [float("nan")])[0], f["sc"], f["dmu"], float(f["valid"]), float(f["at_limit"])])
    res["toys"] = {}
    for shape, a in toys.items():
        a = np.array(a)
        ok = (a[:, 4] > 0) & np.isfinite(a[:, 1]) & (a[:, 1] > 0)
        pull = (a[ok, 0] - ntrue) / a[ok, 1]
        res["toys"][shape] = {"n": int(len(a)), "ok": int(ok.sum()), "bias": float(a[ok, 0].mean() / ntrue - 1), "bias_err": float(a[ok, 0].std() / ntrue / math.sqrt(max(ok.sum(), 1))),
                              "bias_med": float(np.median(a[ok, 0]) / ntrue - 1), "rms_rel": float(a[ok, 0].std() / ntrue), "err_rel": float(np.median(a[ok, 1]) / ntrue),
                              "pull_mean": float(pull.mean()), "pull_rms": float(pull.std()), "cov68": float((np.abs(pull) < 1).mean()), "sc_mean": float(a[ok, 2].mean()),
                              "sc_rms": float(a[ok, 2].std()), "at_limit": float(a[:, 5].mean())}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ntoy", type=int, default=200); ap.add_argument("--nproc", type=int, default=7)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "out", "diag_T.json")); ap.add_argument("--only", default="")
    a = ap.parse_args()
    tasks = [(sp, ib, cat, a.ntoy) for sp in ("D0", "Dplus") for ib in range(3) for cat in CATS]
    if a.only:
        tasks = [t for t in tasks if a.only in f"{t[0]}|{F.BINS[t[1]]}|{t[2]}"]
    t0 = time.time()
    out = []
    with Pool(a.nproc, initializer=_init) as p:
        for i, r in enumerate(p.imap_unordered(task, tasks)):
            out.append(r)
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)} {r['sp']} {r['bin']} {r['cat']}", flush=True)
    out.sort(key=lambda r: (r["sp"], r["bin"], list(CATS).index(r["cat"])))
    json.dump({"res": out}, open(a.out, "w"))
    print(f"zapisano {a.out}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
