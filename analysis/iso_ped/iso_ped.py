#!/usr/bin/env python3
"""ISO_ped — izolacja D0 z odjętym piedestałem (R4). Wejście: B_trig/TRKJET2/runs/lxplus/tj2_*.npz (tylko odczyt).
Kryteria: KRYTERIA.md (zapisane przed uruchomieniem). Wynik: out/iso_ped.json."""
import json, os, sys, time
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "B_trig", "TRKJET2"))
import analyze as A  # noqa: E402
from iso_common import C, DR, fit_pf, sbs_w, _bkg_yield, E_DM, _WIN, chi2_const  # noqa: E402

IN = os.path.join(CW, "B_trig", "TRKJET2", "runs", "lxplus")
OUT = os.path.join(HERE, "out")
CLS = ["9:MB", "8:all", "8:J12only", "8:MBTSbased", "8:other", "8:notJ12"]
BINS = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30), "3-30": (3, 30)}
VARS = ["raw", "V1", "V2", "V3", "V1r"]
NE = [0, 10, 20, 30, 40, 50, 60, 80, 100, 150, 1e9]
PE = [0, 0.3, 0.6, 1.0, 1.5, 2.2, 3.2, 5.0, 1e9]
NB = int(os.environ.get("ISO_NB", 300))
NDM = len(E_DM) - 1


def ped_map(D, m):
    n = np.clip(np.digitize(D["ntrk05"][m], NE) - 1, 0, len(NE) - 2)
    p = np.minimum(D["perp_avg"][m], 10.0)
    s = np.bincount(n, weights=p, minlength=len(NE) - 1); c = np.bincount(n, minlength=len(NE) - 1)
    tot = p.mean()
    return np.where(c >= 30, s / np.maximum(c, 1), tot), c


def add_vars(D, pmap):
    n = np.clip(np.digitize(D["ntrk05"], NE) - 1, 0, len(NE) - 2)
    D["raw"] = D["cone"]
    D["V1"] = D["cone"] - D["perp_avg"]
    D["V2"] = D["cone"] - pmap[n]
    D["V3"] = D["cone"] / (D["perp_avg"] + 1.0)
    D["V1r"] = D["cone"] - D["rnd"].astype(np.float64)


def sbs_h(hp, hf):
    Bp, eBp, Np = _bkg_yield(hp, E_DM); Bf, eBf, Nf = _bkg_yield(hf, E_DM)
    Sp, Sf = Np - Bp, Nf - Bf
    vp, vf = Np + eBp ** 2, Nf + eBf ** 2
    T = Sp + Sf
    if T <= 0:
        return np.nan, np.nan, float(T), float(np.sqrt(vp + vf))
    return float(Sp / T), float(np.sqrt(Sf ** 2 * vp + Sp ** 2 * vf) / T ** 2), float(T), float(np.sqrt(vp + vf))


def _fit(job):
    key, dp, df_, shape = job
    try:
        return key, fit_pf(dp, df_, shape, interval=False)
    except Exception as ex:  # noqa: BLE001
        return key, {"eps": np.nan, "err": np.nan, "N": 0.0, "N_err": 0.0, "det": False, "valid": False, "method": "fit", "reason": str(ex)}


G = {}


def _boot(seed):
    """Jedna replika: Poisson po zdarzeniach (dane per run, MC fd, MC mb). Zwraca {(cls,bin): (eps_S[var], epsMC[var], epsB[var], epsBmc[var])}."""
    rng = np.random.default_rng(1000 + seed)
    wev = {k: rng.poisson(1.0, n).astype(float) for k, n in G["nev"].items()}
    o = {}
    for key, c in G["cells"].items():
        w = wev[c["src"]][c["ev"]]
        eS = []
        for v in VARS:
            f = c["flag"][v]
            hp = np.bincount(c["ib"][f], weights=w[f], minlength=NDM); hf = np.bincount(c["ib"][~f], weights=w[~f], minlength=NDM)
            eS.append(sbs_h(hp, hf)[0])
        wm = wev["fd"][c["mev"]] * c["mw"]
        eM = [float((wm * c["mflag"][v]).sum() / wm.sum()) for v in VARS]
        wb = wev[c["src"]][c["bev"]]
        eB = [float((wb * c["bflag"][v]).sum() / max(wb.sum(), 1e-9)) for v in VARS]
        wq = wev["mb"][c["qev"]]
        eQ = [float((wq * c["qflag"][v]).sum() / max(wq.sum(), 1e-9)) for v in VARS]
        o[key] = (eS, eM, eB, eQ)
    return o


def evidx(D, m, mc=False):
    k = D["event"][m].astype(np.int64)
    if mc:
        k = k + D["file_idx"][m].astype(np.int64) * (1 << 40)
    return k


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    D9, _ = A.load(IN, "data_427929"); D8, _ = A.load(IN, "data_428770")
    M, _ = A.load(IN, "fd_nd"); MB, _ = A.load(IN, "mb_nd")
    S, _ = A.make_classes(D9, D8)
    c8 = A.classes(D8, 428770); S["8:other"] = (D8, c8["other"])
    mm0 = A.mc_sel(M)
    shape = DR.shape_from_mc(M["dm"][mm0] * 1e3)
    # mapy piedestału P(n)
    maps = {}
    for nm, D, cm in (("427929", D9, S["9:MB"][1]), ("428770", D8, S["8:all"][1])):
        maps[nm], cnt = ped_map(D, A.base(D) & cm & (D["dm"] > 0.1485)); maps[nm + "_n"] = cnt
    maps["MC"], cnt = ped_map(M, A.base(M)); maps["MC_n"] = cnt
    add_vars(D9, maps["427929"]); add_vars(D8, maps["428770"]); add_vars(M, maps["MC"]); add_vars(MB, maps["MC"])
    R = {"maps": {k: [float(x) for x in v] for k, v in maps.items()}, "NE": NE, "thr": {}, "ref": {}, "cells": {}, "ped": {}, "nboot": NB}
    # progi: równa eps_S^MC (bez wag) per bin pT
    for b, (lo, hi) in BINS.items():
        mm = A.mc_sel(M, lo, hi)
        e0 = float((M["raw"][mm] < 2.0).mean())
        R["thr"][b] = {"eps_mc": e0, "n_mc": int(mm.sum()), "raw": 2.0}
        for v in VARS[1:]:
            R["thr"][b][v] = float(np.quantile(M[v][mm], e0))
    print("progi", json.dumps(R["thr"], indent=0)[:1500], flush=True)
    # zdarzenia → indeksy bootstrapu
    def uq(keys):
        u, inv = np.unique(keys, return_inverse=True)
        return len(u), inv
    full = {}
    for src, D, mc in (("d9", D9, False), ("d8", D8, False), ("fd", M, True), ("mb", MB, True)):
        n, inv = uq(evidx(D, np.ones(len(D["m"]), bool), mc)); full[src] = inv; G.setdefault("nev", {})[src] = n
    jobs, cells = [], {}
    # krok 0: odniesienie (gałąź ptcone40 < 2, MC bez wag, pT 3–30)
    for name in ("9:MB", "8:all", "8:J12only"):
        D, cm = S[name]; m = A.base(D) & cm
        jobs.append((("ref", name), D["dm"][m & D["_pas"]], D["dm"][m & ~D["_pas"]], shape))
        R["ref"][name] = {"eps_mc": float(M["_pas"][mm0].mean())}
    for name in CLS:
        D, cm = S[name]; src = "d9" if name.startswith("9") else "d8"
        for b, (lo, hi) in BINS.items():
            m = A.base(D, lo, hi) & cm
            mm = A.mc_sel(M, lo, hi)
            wsig = sbs_w(D["dm"][m], D["pt"][m], D["_pas"][m])
            mw = A.mc_pt_w(wsig, D["pt"][m], M["pt"][mm])
            if mw.sum() <= 0:
                mw = np.ones(mm.sum())
            dm = D["dm"][m] * 1e3
            inr = (dm >= E_DM[0]) & (dm < E_DM[-1])
            ib = np.clip(np.digitize(dm, E_DM) - 1, 0, NDM - 1)
            # tło: pasma boczne masy D0, jeden wpis na D0
            sb = (np.abs(D["m"] - C.M0_D0) > 0.050) & (np.abs(D["m"] - C.M0_D0) < 0.075) & (D["score"] > A.STHR) & (D["pt"] >= lo) & (D["pt"] < hi) & D["dok"] & cm
            idx = np.flatnonzero(sb)
            _, first = np.unique(D["file_idx"][idx].astype(np.int64) * (1 << 40) + D["row"][idx].astype(np.int64), return_index=True)
            idx = idx[first]
            qb = (np.abs(MB["m"] - C.M0_D0) > 0.050) & (np.abs(MB["m"] - C.M0_D0) < 0.075) & (MB["score"] > A.STHR) & (MB["pt"] >= lo) & (MB["pt"] < hi) & MB["dok"] & (MB["truth"] >= 8)
            qidx = np.flatnonzero(qb)
            _, first = np.unique(MB["file_idx"][qidx].astype(np.int64) * (1 << 40) + MB["row"][qidx].astype(np.int64), return_index=True)
            qidx = qidx[first]
            c = {"src": src, "ev": full[src][m][inr], "ib": ib[inr], "flag": {}, "mev": full["fd"][mm], "mw": mw, "mflag": {},
                 "bev": full[src][idx], "bflag": {}, "qev": full["mb"][qidx], "qflag": {}}
            r = {"ncand": int(m.sum()), "nb": int(len(idx)), "nb_mc": int(len(qidx)), "var": {}}
            for v in VARS:
                t = R["thr"][b][v]
                f = D[v][m] < t
                c["flag"][v] = f[inr]; c["mflag"][v] = (M[v][mm] < t).astype(float)
                c["bflag"][v] = (D[v][idx] < t).astype(float); c["qflag"][v] = (MB[v][qidx] < t).astype(float)
                hp = np.histogram(dm[f], E_DM)[0]; hf = np.histogram(dm[~f], E_DM)[0]
                e, er, N, eN = sbs_h(hp, hf)
                emc = float((mw * c["mflag"][v]).sum() / mw.sum())
                r["var"][v] = {"thr": t, "sbs": {"eps": e, "err": er, "N": N, "N_err": eN}, "eps_mc": emc, "eps_mc_unw": float(c["mflag"][v].mean()),
                               "epsB": float(c["bflag"][v].mean()) if len(idx) else np.nan, "epsB_mc": float(c["qflag"][v].mean()) if len(qidx) else np.nan}
                jobs.append((("cell", name, b, v), D["dm"][m][f], D["dm"][m][~f], shape))
            cells[(name, b)] = c; R["cells"][f"{name}|{b}"] = r
            # --- rozkład piedestału sygnału (krok 3)
            pm = M["perp_avg"][mm]; pd_ = D["perp_avg"][m]
            def wstat(w, x):
                T = w.sum(); mu = (w * x).sum() / T
                return float(mu), float(np.sqrt((w ** 2 * (x - mu) ** 2).sum()) / abs(T))
            P = {}
            for lab, xd, xm in (("perp", np.minimum(pd_, 10.0), np.minimum(pm, 10.0)),
                                ("rnd", np.minimum(D["rnd"][m], 10.0), np.minimum(M["rnd"][mm], 10.0)),
                                ("near", np.clip(D["cone"][m] - pd_, -10, 30), np.clip(M["cone"][mm] - pm, -10, 30)),
                                ("cone", np.minimum(D["cone"][m], 30.0), np.minimum(M["cone"][mm], 30.0)),
                                ("ntrk05", np.minimum(D["ntrk05"][m], 300.0), np.minimum(M["ntrk05"][mm], 300.0))):
                if wsig.sum() > 0:
                    P[lab] = {"data": wstat(wsig, xd), "mc": wstat(mw, xm)}
            for t in (0.5, 1.0, 2.0):
                if wsig.sum() > 0:
                    P[f"f<{t}"] = {"data": wstat(wsig, (pd_ < t).astype(float)), "mc": wstat(mw, (pm < t).astype(float))}
            if wsig.sum() > 0:
                hd = np.histogram(pd_, PE, weights=wsig)[0]; vd = np.histogram(pd_, PE, weights=wsig ** 2)[0]
                hm = np.histogram(pm, PE, weights=mw)[0]; vm = np.histogram(pm, PE, weights=mw ** 2)[0]
                fd, fm = hd / hd.sum(), hm / hm.sum()
                ed, em = np.sqrt(vd) / hd.sum(), np.sqrt(vm) / hm.sum()
                P["hist"] = {"edges": PE, "data": fd.tolist(), "data_err": ed.tolist(), "mc": fm.tolist(), "mc_err": em.tolist(),
                             "chi2": float(((fd - fm) ** 2 / (ed ** 2 + em ** 2 + 1e-12)).sum()), "ndf": len(PE) - 2}
            R["ped"][f"{name}|{b}"] = P
    print(f"przygotowane: {len(jobs)} fitów, {len(cells)} komórek, {time.time()-t0:.0f}s", flush=True)
    G["cells"] = cells
    import multiprocessing as mp, pickle
    cache = os.path.join(OUT, "fits_cache.pkl")                 # fity są deterministyczne; cache tylko po to, by nie liczyć ich dwa razy
    with mp.get_context("fork").Pool(7) as Pl:                  # fork: repliki bootstrapu czytają G (Python 3.14 domyślnie forkserver)
        if os.path.exists(cache):
            fits = pickle.load(open(cache, "rb"))
            assert [k for k, _ in fits] == [j[0] for j in jobs]
        else:
            fits = Pl.map(_fit, jobs, chunksize=2)
            pickle.dump(fits, open(cache, "wb"))
        print(f"fity: {time.time()-t0:.0f}s", flush=True)
        boots = Pl.map(_boot, range(NB), chunksize=4)
    print(f"bootstrap: {time.time()-t0:.0f}s", flush=True)
    for key, f in fits:
        f = {k: (float(x) if isinstance(x, (float, np.floating)) else x) for k, x in f.items()}
        if key[0] == "ref":
            R["ref"][key[1]]["fit"] = f; R["ref"][key[1]]["SF"] = f["eps"] / R["ref"][key[1]]["eps_mc"]; R["ref"][key[1]]["SF_err"] = f["err"] / R["ref"][key[1]]["eps_mc"]
        else:
            R["cells"][f"{key[1]}|{key[2]}"]["var"][key[3]]["fit"] = f
    # --- złożenie: SF, D, bootstrap
    for (name, b) in cells:
        r = R["cells"][f"{name}|{b}"]
        bs = np.array([[x[(name, b)][i] for i in range(4)] for x in boots])          # (NB, 4, nvar)
        eS, eM, eB, eQ = bs[:, 0], bs[:, 1], bs[:, 2], bs[:, 3]
        sfb = eS / eM
        N, eN = r["var"]["raw"]["sbs"]["N"], r["var"]["raw"]["sbs"]["N_err"]
        r["N_sig"], r["N_err"] = N, eN
        r["quality"] = "ok" if (N >= 150 and eN > 0 and N / eN >= 3) else "niemiarodajna"
        for i, v in enumerate(VARS):
            x = r["var"][v]; f = x.get("fit", {})
            okf = bool(f.get("det") and f.get("valid") and np.isfinite(f.get("eps", np.nan)) and np.isfinite(f.get("err", np.nan)) and f.get("err", 0) > 0)
            x["fit_ok"] = okf
            e, er, meth = (f["eps"], f["err"], "fit") if okf else (x["sbs"]["eps"], x["sbs"]["err"], "SBS")
            x["eps"], x["eps_err"], x["method"] = e, er, meth
            x["SF"] = e / x["eps_mc"]; x["SF_sbs"] = x["sbs"]["eps"] / x["eps_mc"]
            x["SF_err_boot"] = float(np.nanstd(sfb[:, i])); x["SF_err"] = float(np.hypot(er / x["eps_mc"], x["SF"] * np.nanstd(eM[:, i]) / x["eps_mc"]))
            x["epsB_err"] = float(np.nanstd(eB[:, i])); x["epsB_mc_err"] = float(np.nanstd(eQ[:, i]))
            if v != "raw":
                x0 = r["var"]["raw"]
                Db = np.abs(sfb[:, 0] - 1) - np.abs(sfb[:, i] - 1)
                x["D"] = abs(x0["SF"] - 1) - abs(x["SF"] - 1); x["D_sbs"] = abs(x0["SF_sbs"] - 1) - abs(x["SF_sbs"] - 1)
                x["D_err"] = float(np.nanstd(Db)); x["dSF_err"] = float(np.nanstd(sfb[:, i] - sfb[:, 0]))
                x["rB"] = x["epsB"] / x0["epsB"] if x0["epsB"] > 0 else np.nan
                x["rB_err"] = float(np.nanstd(eB[:, i] / np.maximum(eB[:, 0], 1e-9)))
                x["rB_mc"] = x["epsB_mc"] / x0["epsB_mc"] if x0["epsB_mc"] > 0 else np.nan
                zg = (eS[:, i] / eS[:, 0]) / np.sqrt(eB[:, i] / np.maximum(eB[:, 0], 1e-9))
                x["Zgain"] = (x["eps"] / x0["eps"]) / np.sqrt(x["rB"]) if x["rB"] > 0 else np.nan
                x["Zgain_err"] = float(np.nanstd(zg))
                x["Zgain_mc"] = (x["eps_mc"] / x0["eps_mc"]) / np.sqrt(x["rB_mc"]) if x["rB_mc"] > 0 else np.nan
    # stałość SF po klasach rozłącznych
    R["const"] = {}
    for b in BINS:
        for v in VARS:
            xs = [R["cells"][f"{n}|{b}"]["var"][v] for n in ("9:MB", "8:J12only", "8:MBTSbased", "8:other")]
            R["const"][f"{b}|{v}"] = chi2_const([x["SF"] for x in xs], [x["SF_err"] for x in xs])
    json.dump(R, open(os.path.join(OUT, "iso_ped.json"), "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"koniec {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
