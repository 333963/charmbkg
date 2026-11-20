#!/usr/bin/env python3
"""FITBIAS — warianty modelu M0–M9 (KRYTERIA §2.1) i próby wyboru (§2.2): asimov (T), mc (W, mb_nd), toys (pasma danych + wstrzyknięcie), data (po wyborze).
  fb_var.py asimov|mc|toys|data [--ntoy 150] [--nproc 7] [--chosen D0=M4,Dplus=M4]
Wyjście: out/var_<tryb>.json"""
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

WP = os.path.join(HERE, "..", "..", "WPSCAN_2026-11-19", "runs", "lxplus", "wps_hist.npz")
TD = os.path.join(HERE, "..", "..", "MC_templates", "runs", "lxplus")
TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
SET = {"M0": dict(shape="free", peak=False, dn=0, rng="fit"), "M1": dict(shape="gauss", peak=False, dn=0, rng="fit"), "M2": dict(shape="fix", peak=False, dn=0, rng="fit"),
       "M3": dict(shape="free", peak=True, dn=0, rng="fit"), "M4": dict(shape="gauss", peak=True, dn=0, rng="fit"), "M5": dict(shape="fix", peak=True, dn=0, rng="fit"),
       "M6": dict(shape="free", peak=False, dn=1, rng="fit"), "M7": dict(shape="free", peak=False, dn=0, rng="alt"), "M8": dict(shape="gauss", peak=True, dn=0, rng="alt"),
       "M9": dict(shape="gauss", peak=True, dn=1, rng="fit")}
GC = {"D0": (0.0008, 0.10), "Dplus": (0.0005, 0.10)}
GROUPS = {"Sig": [0, 1], "Refl": [2, 3], "CM": [4, 5, 6, 7, 8, 9, 10], "DMM": [11, 12], "PC": [13, 14, 15, 16, 17, 18, 19], "Comb": [20, 21, 22, 23]}
TCATS = {"WP95": [1, 2, 3], "WP80": [2, 3], "WP50": [3], "L80-50": [2]}
Z = INFO = M = FD = None
SRC = {"427929": "427929|MB427929", "428770": "428770|all"}
KSCALE_RETRACTED = "D0=1.364,Dplus=1.068"  # WYCOFANE 20.11.2026 (cver/RAPORT.md pkt 5); nominał = bez skali κ (1.0) — patrz main()


def _init():
    global Z, INFO, M, FD
    Z = np.load(WP, allow_pickle=False)
    INFO = json.loads(str(Z["meta"]))
    M = np.load(os.path.join(TD, "mctmpl_mb_nd.npz"))
    FD = np.load(os.path.join(TD, "mctmpl_fd_nd.npz"))


# ------------------------------------------------------------------------------------------------ histogramy
def full_edges(sp):
    lo, hi = F.SP[sp]["mr"]
    return np.arange(lo, hi + F.BW / 2, F.BW)


def full5(h1):
    """1 MeV, pełny zakres → 5 MeV, pełny zakres."""
    return h1.reshape(h1.shape[:-1] + (h1.shape[-1] // 5, 5)).sum(-1).astype(np.float64)


def cutr(h5, sp, rng):
    lo = F.SP[sp]["mr"][0]
    i0, i1 = int(round((rng[0] - lo) / F.BW)), int(round((rng[1] - lo) / F.BW))
    return h5[..., i0:i1]


def tcomp(Zt, sp, ib, sbs):
    h = Zt["d0_m" if sp == "D0" else "dp_m"][ib]
    h = h.sum(2) if sp == "D0" else h.sum(2)[:, :, 1]
    h = h[sbs].sum(0).astype(float)
    return {g: h[ix].sum(0) for g, ix in GROUPS.items()}


def kappa_nodes(sp, ib):
    """κ(ε) = N_(DMM+CM)/N_Sig skumulowane w węzłach ε = 0,95; 0,80; 0,50 (mb_nd, T)."""
    o = {}
    for e, cat in ((0.95, "WP95"), (0.80, "WP80"), (0.50, "WP50")):
        c = tcomp(M, sp, ib, TCATS[cat])
        o[e] = float((c["DMM"].sum() + c["CM"].sum()) / c["Sig"].sum())
    return o


def kappa_eps(nodes, e):
    xs = [0.50, 0.80, 0.95]
    return float(np.interp(e, xs, [nodes[x] for x in xs]))


def peak_tmpl(sp, ib):
    c = tcomp(M, sp, ib, TCATS["WP80"])
    return F.TmplE(full_edges(sp), c["DMM"] + c["CM"])


def wsel(h, jlo, jhi=None):
    return h[..., jlo + 1:(None if jhi is None else jhi + 1), :].sum(-2)


def wcat(sp, bn, elo, ehi=None):
    """Kategoria W: ε_lo ≥ score-WP > ε_hi (ehi=None → skumulowana). → składniki modelu."""
    c = F.SP[sp]
    lo = c["mr"][0]
    key = f"{sp}|{TAG[sp]}|{bn}"
    ed = INFO["eps_desc"]
    jlo, jhi = ed.index(elo), (None if ehi is None else ed.index(ehi))
    eps = Z[key + "|eps"]
    sw = wsel(Z[key + "|fd|sig_w"], jlo, jhi)
    o = {"sigT": F.tmpl1(sw, lo), "refT": None, "frefl": 0.0, "j": (jlo, jhi), "key": key, "e": float(sw.sum() / Z[key + "|fd|sig_w"].sum())}
    if c["refl"]:
        rw = wsel(Z[key + "|fd|refl_w"], jlo, jhi)
        o["refT"], o["frefl"] = F.tmpl1(rw, lo), float(rw.sum() / sw.sum())
    ib = F.BINS.index(bn)
    nodes = kappa_nodes(sp, ib)
    e_lo = float(eps[jlo])
    if jhi is None:
        o["kap"] = kappa_eps(nodes, elo)
    else:
        e_hi = float(eps[jhi])
        o["kap"] = (kappa_eps(nodes, elo) * e_lo - kappa_eps(nodes, ehi) * e_hi) / (e_lo - e_hi)
    o["peakT"] = peak_tmpl(sp, ib)
    return o


def wdata(cat, src):
    """hist 5 MeV pełnego zakresu; src = 'mb|all', 'mb|sig', '427929|MB427929', ... (dane: [nrep+1, nfull])."""
    return full5(wsel(Z[f"{cat['key']}|{src}"], *cat["j"]))


def fit_set(H5, cats, sp, mname, start=None, hesse=True):
    """H5 [ncat, nfull] (pełny zakres 5 MeV); cats = lista kategorii (sigT, refT, frefl, kap, peakT)."""
    c, S = F.SP[sp], SET[mname]
    rng = c[S["rng"]]
    edges = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    tied = []
    for k in cats:
        t = [(k["refT"], k["frefl"])] if k["refT"] is not None else []
        if S["peak"]:
            t.append((k["peakT"], k["kap"]))
        tied.append(t)
    f = F.fit(cutr(np.asarray(H5, float), sp, rng), [k["sigT"] for k in cats], tied, edges, c["m0"], c["nbkg"] + S["dn"], shape=S["shape"], gc=GC[sp], win=c["win"], start=start, hesse=hesse)
    f.pop("mu", None)
    return f


def sres(f):
    return {"N": f["N"][0], "err": (f.get("err") or [None])[0], "valid": f["valid"], "at_limit": f["at_limit"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"],
            "bkg_win": f["bkg_win"][0]}


def jres(f):
    NA, NB = f["N"]
    E = NA + NB
    o = {"NA": NA, "NB": NB, "E": E, "f": NA / E if E > 0 else None, "valid": f["valid"], "at_limit": f["at_limit"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"]}
    if "cov_N" in f and E > 0:
        C = np.array(f["cov_N"])
        g = np.array([NB / E ** 2, -NA / E ** 2])
        o.update({"NA_err": math.sqrt(max(C[0, 0], 0)), "NB_err": math.sqrt(max(C[1, 1], 0)), "E_err": math.sqrt(max(C.sum(), 0)), "f_err": math.sqrt(max(g @ C @ g, 0))})
    return o


# ------------------------------------------------------------------------------------------------ 1) Asimov T
def t_asimov(t):
    sp, ib = t
    c = F.SP[sp]
    fe = full_edges(sp)
    ctr = 0.5 * (fe[1:] + fe[:-1])
    cats, pseudo, truth = {}, {}, {}
    pk = peak_tmpl(sp, ib)
    for cn in ("WP80", "WP50", "L80-50"):
        mb, fd = tcomp(M, sp, ib, TCATS[cn]), tcomp(FD, sp, ib, TCATS[cn])
        k = {"sigT": F.TmplE(fe, fd["Sig"]), "refT": F.TmplE(fe, fd["Refl"]) if c["refl"] else None, "frefl": float(fd["Refl"].sum() / fd["Sig"].sum()) if c["refl"] else 0.0,
             "kap": float((mb["DMM"].sum() + mb["CM"].sum()) / mb["Sig"].sum()), "peakT": pk}
        hb = mb["PC"] + mb["Comb"]
        fa = F.fit(hb[None], [k["sigT"]], [[]], fe, c["m0"], c["nbkg"] + 2, shape="fix", win=c["win"], fixN0=True)
        cats[cn] = k
        truth[cn] = float(mb["Sig"].sum())
        pseudo[cn] = mb["Sig"] + mb["Refl"] + mb["DMM"] + mb["CM"] + fa["mu"][0]
    res = {"sp": sp, "bin": F.BINS[ib], "truth": truth, "kap": {cn: cats[cn]["kap"] for cn in cats}, "S": {}, "J": {}}
    for m in SET:
        res["S"][m] = {cn: sres(fit_set(pseudo[cn][None], [cats[cn]], sp, m)) for cn in cats}
        res["J"][m] = jres(fit_set(np.stack([pseudo["WP50"], pseudo["L80-50"]]), [cats["WP50"], cats["L80-50"]], sp, m))
    return res


# ------------------------------------------------------------------------------------------------ 2) domknięcie mb_nd (W)
def t_mc(t):
    sp, bn, nrep = t
    ed = INFO["eps_desc"]
    res = {"sp": sp, "bin": bn, "S": {}, "J": {}}
    single = {f"WP{int(round(e * 100))}": (e, None) for e in ed}
    single.update({"B50": (0.7, 0.5), "B40": (0.7, 0.4), "B35": (0.7, 0.35), "L90-70": (0.9, 0.7)})
    for cn, (elo, ehi) in single.items():
        k = wcat(sp, bn, elo, ehi)
        H = wdata(k, "mb|all")
        nt = float(wdata(k, "mb|sig").sum())
        o = {"N_true": nt, "eps": k["e"], "kap": k["kap"], "frefl": k["frefl"]}
        for m in SET:
            o[m] = sres(fit_set(H[None], [k], sp, m))
        res["S"][cn] = o
    for et in (0.5, 0.4, 0.35):
        a, b, w70 = wcat(sp, bn, et), wcat(sp, bn, 0.7, et), wcat(sp, bn, 0.7)
        HA, HB = wdata(a, "mb|all"), wdata(b, "mb|all")
        nA, nB = float(wdata(a, "mb|sig").sum()), float(wdata(b, "mb|sig").sum())
        o = {"NA_true": nA, "NB_true": nB, "f_true": nA / (nA + nB), "f_mc": a["e"] / (a["e"] + b["e"])}
        rng = np.random.default_rng(777 + int(et * 100) + 13 * F.BINS.index(bn) + (500 if sp == "D0" else 0))
        reps = [(rng.poisson(HA).astype(float), rng.poisson(HB).astype(float)) for _ in range(nrep)]
        for m in SET:
            f = fit_set(np.stack([HA, HB]), [a, b], sp, m)
            r = jres(f)
            s = sres(fit_set((HA + HB)[None], [w70], sp, m))
            fr, Er, Sr = [], [], []
            for ra, rb in reps:
                fj = jres(fit_set(np.stack([ra, rb]), [a, b], sp, m, hesse=False))
                if fj["valid"] and fj["f"] is not None:
                    fr.append(fj["f"]); Er.append(fj["E"])
                fs = fit_set((ra + rb)[None], [w70], sp, m, hesse=False)
                if fs["valid"]:
                    Sr.append(fs["N"][0])
            r.update({"f_rms": float(np.std(fr)) if len(fr) > 5 else None, "E_rms": float(np.std(Er)) if len(Er) > 5 else None, "N70": s["N"], "N70_err": s["err"],
                      "N70_rms": float(np.std(Sr)) if len(Sr) > 5 else None, "n_rep_ok": len(fr)})
            o[m] = r
        res["J"][str(et)] = o
    return res


# ------------------------------------------------------------------------------------------------ 3) zabawki: pasma danych + wstrzyknięcie
def t_toys(t):
    sp, bn, run, level, ntoy = t
    c = F.SP[sp]
    fe = full_edges(sp)
    ctr = 0.5 * (fe[1:] + fe[:-1])
    blind = np.abs(ctr - c["m0"]) < c["win"]
    a, b, w70 = wcat(sp, bn, 0.5), wcat(sp, bn, 0.7, 0.5), wcat(sp, bn, 0.7)
    gen, Ninj, info = [], [], {}
    for nm, k in (("A", a), ("B", b)):
        Hd = wdata(k, SRC[run])[0]                                   # replika 0 = nominał; UŻYWANE TYLKO PASMA (okno zaślepione)
        fb = F.fit(Hd[None], [k["sigT"]], [[]], fe, c["m0"], c["nbkg"], shape="fix", win=c["win"], blind=blind, hesse=False)
        bkg = fb["mu"][0]
        Hmb, Smb = wdata(k, "mb|all"), wdata(k, "mb|sig")
        n = level * float(Smb.sum()) * float(Hd[~blind].sum()) / float(Hmb[~blind].sum())
        gen.append(bkg); Ninj.append(n)
        info[nm] = {"N_inj": n, "bkg_win": float(bkg[blind].sum()), "sb_data": float(Hd[~blind].sum()), "sb_mb": float(Hmb[~blind].sum()), "bkg_fit_valid": fb["valid"],
                    "chi2_sb": fb["chi2"], "ndf_sb": fb["ndf"]}
    sw = lambda k, dm, s: k["sigT"].bins(fe, c["m0"], dm, s) + (k["frefl"] * k["refT"].bins(fe, c["m0"], dm, 1.0) if k["refT"] is not None else 0.0) + k["kap"] * k["peakT"].bins(fe, c["m0"], dm, 1.0)  # noqa: E731
    s0 = a["sigT"].bins(fe, c["m0"])
    info["SB_win_A"] = Ninj[0] * float(s0[blind].sum()) / max(info["A"]["bkg_win"], 1e-9)
    s0b = b["sigT"].bins(fe, c["m0"])
    info["SB_win_70"] = (Ninj[0] * float(s0[blind].sum()) + Ninj[1] * float(s0b[blind].sum())) / max(info["A"]["bkg_win"] + info["B"]["bkg_win"], 1e-9)
    info["SB_win_B"] = Ninj[1] * float(s0b[blind].sum()) / max(info["B"]["bkg_win"], 1e-9)
    rng = np.random.default_rng(424242 + 1000 * F.BINS.index(bn) + (5000 if sp == "D0" else 0) + (70000 if run == "428770" else 0) + int(level * 10))
    rows = {m: [] for m in SET}
    for _ in range(ntoy):
        dm, s = rng.normal(0, GC[sp][0]), rng.normal(1, GC[sp][1])
        s = min(max(s, 0.75), 1.4)
        HA = rng.poisson(gen[0] + Ninj[0] * sw(a, dm, s)).astype(float)
        HB = rng.poisson(gen[1] + Ninj[1] * sw(b, dm, s)).astype(float)
        for m in SET:
            fs = sres(fit_set((HA + HB)[None], [w70], sp, m))
            fj = jres(fit_set(np.stack([HA, HB]), [a, b], sp, m))
            nan = float("nan")
            rows[m].append([fs["N"], fs["err"] if fs["err"] is not None else nan, float(fs["valid"]), float(fs["at_limit"]), fs["sc"], fs["dmu"],
                            fj["NA"], fj["NB"], fj.get("E_err", nan), fj["f"] if fj["f"] is not None else nan, fj.get("f_err", nan), float(fj["valid"]), float(fj["at_limit"]), dm, s])
    return {"sp": sp, "bin": bn, "run": run, "level": level, "info": info, "N_inj": Ninj, "rows": {m: np.array(v).tolist() for m, v in rows.items()}}


# ------------------------------------------------------------------------------------------------ 4) dane (po wyborze)
def t_data(t):
    sp, bn, run, chosen, nrep, ksc = t
    ed = INFO["eps_desc"]
    res = {"sp": sp, "bin": bn, "run": run, "kscale": ksc[sp], "S": {}, "J": {}}
    _w = wcat

    def wcat_s(*a_, **k_):                                           # κ_dane = κ_MC × r_BR (KRYTERIA cz. 3)
        k = _w(*a_, **k_)
        k["kap"] *= ksc[sp]
        return k
    rep_m = ["M0"] + ([chosen] if chosen != "M0" else [])
    for e in ed:
        cn = f"WP{int(round(e * 100))}"
        k = wcat_s(sp, bn, e)
        D = wdata(k, SRC[run])                                       # [nrep+1, nfull]
        o = {"eps": k["e"], "n_rows": float(D[0].sum())}
        for m in (SET if e == 0.7 else rep_m):
            f = fit_set(D[0][None], [k], sp, m)
            r = sres(f)
            if m in rep_m and e in (0.7, 0.5, 0.35):
                bN = []
                for i in range(1, min(nrep, D.shape[0] - 1) + 1):
                    fb = fit_set(D[i][None], [k], sp, m, start=f["x"], hesse=False)
                    bN.append(fb["N"][0] if fb["valid"] else None)
                r["boot"] = bN
            o[m] = r
        res["S"][cn] = o
    for et in (0.5, 0.4, 0.35):
        a, b = wcat_s(sp, bn, et), wcat_s(sp, bn, 0.7, et)
        DA, DB = wdata(a, SRC[run]), wdata(b, SRC[run])
        o = {"f_mc": a["e"] / (a["e"] + b["e"])}
        for m in SET:
            f = fit_set(np.stack([DA[0], DB[0]]), [a, b], sp, m)
            r = jres(f)
            if m in rep_m:
                bf, bE = [], []
                for i in range(1, min(nrep, DA.shape[0] - 1) + 1):
                    fb = jres(fit_set(np.stack([DA[i], DB[i]]), [a, b], sp, m, start=f["x"], hesse=False))
                    ok = fb["valid"] and fb["f"] is not None
                    bf.append(fb["f"] if ok else None); bE.append(fb["E"] if ok else None)
                r["boot_f"], r["boot_E"] = bf, bE
            o[m] = r
        res["J"][str(et)] = o
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["asimov", "mc", "toys", "data"])
    ap.add_argument("--ntoy", type=int, default=150); ap.add_argument("--nproc", type=int, default=7); ap.add_argument("--nrep", type=int, default=30)
    # UWAGA (oznaczone 20.11.2026): wartość domyślna --kscale poniżej (KSCALE_RETRACTED) jest WYCOFANA — BR sygnału w MC = PDG,
    # więc κ_dane = κ_MC, r = 1 (docs/results/cver/RAPORT.md pkt 5). Zostaje bez zmian tylko po to, by ten plik odtwarzał zapisany
    # out/var_data.json; wynik bez skali (out/var_data_k1.json) wymaga jawnego --kscale D0=1.0,Dplus=1.0. Uzyski D0 przekroju
    # końcowego to fity r = 1 z analysis/xsec_final/skrypty/kappa_fit.py (out/kappa_A.json).
    ap.add_argument("--chosen", default=""); ap.add_argument("--kscale", default=KSCALE_RETRACTED, help="skala κ per gatunek; wartość domyślna WYCOFANA 20.11.2026 — nominał: D0=1.0,Dplus=1.0"); ap.add_argument("--only", default=""); ap.add_argument("--out", default="")
    a = ap.parse_args()
    if a.cmd == "data" and a.kscale == KSCALE_RETRACTED:
        print(f"UWAGA: --kscale = {KSCALE_RETRACTED} to WYCOFANA wartość domyślna (cver/RAPORT.md pkt 5); nominał: --kscale D0=1.0,Dplus=1.0", file=sys.stderr)
    t0 = time.time()
    if a.cmd == "asimov":
        fn, tasks = t_asimov, [(sp, ib) for sp in TAG for ib in range(3)]
    elif a.cmd == "mc":
        fn, tasks = t_mc, [(sp, bn, a.nrep) for sp in TAG for bn in F.BINS]
    elif a.cmd == "toys":
        fn, tasks = t_toys, [(sp, bn, run, lv, a.ntoy) for sp in TAG for bn in F.BINS for run in SRC for lv in (1.0, 0.5)]
    else:
        ch = dict(x.split("=") for x in a.chosen.split(","))
        ks = {k: float(v) for k, v in (x.split("=") for x in a.kscale.split(","))}
        fn, tasks = t_data, [(sp, bn, run, ch[sp], a.nrep, ks) for sp in TAG for bn in F.BINS for run in SRC]
    if a.only:
        tasks = [t for t in tasks if a.only in "|".join(map(str, t))]
    out = []
    with Pool(a.nproc, initializer=_init) as p:
        for i, r in enumerate(p.imap_unordered(fn, tasks)):
            out.append(r)
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(tasks)}", flush=True)
    out.sort(key=lambda r: (r["sp"], F.BINS.index(r["bin"]), r.get("run", ""), r.get("level", 0)))
    fn_out = a.out or os.path.join(HERE, "..", "out", f"var_{a.cmd}.json")
    json.dump({"res": out, "chosen": a.chosen, "sec": time.time() - t0}, open(fn_out, "w"))
    print(f"zapisano {fn_out}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
