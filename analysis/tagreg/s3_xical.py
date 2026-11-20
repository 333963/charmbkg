#!/usr/bin/env python3
"""s3 (część C): rozkład Δ = ξ_cal − ξ_AFP w zdarzeniach tagowanych (T1), składowa przypadkowa z mieszania, ε_true cięcia Δ ≤ 0,02.
Bez żadnego fitu masy; próbki: wszystkie zdarzenia / z kandydatem D (slim) / z kandydatem w pasmach bocznych."""
import numpy as np, os, sys, common as C
rng = np.random.default_rng(20261118)
FE = np.unique(np.round(np.concatenate([np.arange(-0.10, 0.3001, 0.02), [0.045]]), 4))   # krawędzie (zawierają 0.02, 0.045, 0.10, 0.20)
MAXPOOL = 1_500_000
E = np.concatenate([[-np.inf], FE, [np.inf]]); NB = len(E) - 1
def reg(lo, hi): return (E[:-1] >= lo - 1e-9) & (E[1:] <= hi + 1e-9)
SR = reg(-np.inf, 0.02); NORM = {"n10": reg(0.10, np.inf), "n045": reg(0.045, np.inf), "n20": reg(0.20, np.inf)}
NCL = len(C.NSEL_EDGES) - 1
NSG = {"nSel<15": (0, 3), "15-30": (3, 6), "30-50": (6, 8), ">=50": (8, 10)}

def hist2(cell, x, ncell):
    k = np.searchsorted(E, x, side="left") - 1          # przedział (lo, hi]
    k = np.clip(k, 0, NB - 1)
    return np.bincount(cell * NB + k, minlength=ncell * NB).reshape(ncell, NB).astype(float)

def template(xc, gc, cc, xp, gp, ng, ncell):
    """T[cell, bin] = Σ_i∈cell ułamek protonów bloku z Δ w binie."""
    T = np.zeros((ncell, NB))
    for g in range(ng):
        mc = np.flatnonzero(gc == g); p = np.sort(xp[gp == g])
        if len(p) == 0 or len(mc) == 0: continue
        x = xc[mc]; c = cc[mc]
        prev = np.zeros(len(x))
        for k in range(NB):
            hi = E[k + 1]
            cur = np.ones(len(x)) if hi == np.inf else 1.0 - np.searchsorted(p, x - hi, side="left") / len(p)   # frakcja z x−ξ ≤ hi
            T[:, k] += np.bincount(c, cur - prev, minlength=ncell); prev = cur
    return T

def eps_from(N, T, cells=None, nboot=200):
    """N, T: [cell, bin]. → dict wariantów normalizacji: ε_raw, ε_true, f_acc(SR), N_acc itd. + błąd stat (Poisson danych, bootstrap)."""
    if cells is not None: N = N[cells]; T = T[cells]
    out = {"N_tot": float(N.sum()), "N_SR": float(N[:, SR].sum()), "eps_raw": float(N[:, SR].sum() / max(N.sum(), 1))}
    for nm, nr in NORM.items():
        def calc(Nx):
            tn = T[:, nr].sum(1); good = tn > 0
            sc = np.where(good, Nx[:, nr].sum(1) / np.where(good, tn, 1), 0.0)
            acc = T * sc[:, None]; a_sr = acc[:, SR].sum(); a_tot = acc.sum()
            n_sr = Nx[:, SR].sum(); n_tot = Nx.sum()
            return (n_sr - a_sr) / (n_tot - a_tot), a_sr / max(n_sr, 1), a_sr, a_tot, (acc.sum(0))
        e, f, a_sr, a_tot, acch = calc(N)
        bs = np.array([calc(rng.poisson(N).astype(float))[:2] for _ in range(nboot)])
        out[nm] = {"eps_true": float(e), "eps_err": float(bs[:, 0].std()), "f_acc_SR": float(f), "f_acc_err": float(bs[:, 1].std()), "N_acc_SR": float(a_sr), "N_acc_tot": float(a_tot),
                   "N_nopool_norm": float(N[T[:, nr].sum(1) <= 0][:, nr].sum())}
        if nm == "n10": out["acc_hist"] = acch
    out["data_hist"] = N.sum(0)
    return out

res = {"edges": E}
for r in C.RUNS:
    d = C.load(r); cl, b = C.classes(d); T0, T1, T2 = C.regions(d); R = {}
    base = d["grl"] & (d["npu"] == 0)
    ncls = np.clip(np.searchsorted(C.NSEL_EDGES, d["nsel"], side="right") - 1, 0, NCL - 1)
    streams = {"all": np.ones(len(base), bool)}
    if r == 428770:
        streams.update({"AFP": b["A"], "nonAFP": (b["J"] | b["M"]) & ~b["A"], "none": cl["none"], "AFPsptrk4_any": b["as_"], "AFPnoalg_any": b["an"]})
    nonafp = b["J"] | b["M"]
    pools = {"nonAFP": nonafp, "MBTS": b["M"]} if r == 428770 else {"all": np.ones(len(base), bool)}
    samples = {"wszystkie": np.ones(len(base), bool), "z kand. D": d["hasD"], "z kand. SB": d["hasSB"], "z kand. D0": d["hasD0"],
               "SB Dplus": d["hasSB_Dplus"], "SB D0": d["hasSB_D0"], "SB Dstar": d["hasSB_Dstar"]}
    for sname, sm in streams.items():
        for smp, pm in samples.items():
            if sname not in ("all", "AFP", "nonAFP", "none") and smp not in ("z kand. D",): continue
            for pname, pool in pools.items():
                if pname != list(pools)[0] and smp not in ("z kand. D",): continue
                key = f"{sname}|{smp}|pool={pname}"; acc = {}
                for s in C.SIDES:
                    tag = base & d[f"tag_{s}"] & sm
                    dat = tag & pm
                    if dat.sum() < 50: continue
                    xi = d[f"xi1_{s}"]; xc_all = d[f"xical_{s}"]
                    e = np.unique(np.quantile(d["lb"][tag], np.linspace(0, 1, 9)[1:-1])); g = np.searchsorted(e, d["lb"], side="right").astype(np.int64); ng = len(e) + 1
                    cell = g * NCL + ncls; ncell = ng * NCL
                    prot = tag & (xc_all > 0.14)                       # protony „niemożliwe” strumienia (z D i bez D)
                    cen = base & pool & pm & (d[f"np_{s}"] == 0)       # pula centralna bez protonu po stronie
                    if cen.sum() > MAXPOOL:
                        idx = np.flatnonzero(cen); keep = rng.choice(idx, MAXPOOL, replace=False); cen = np.zeros(len(base), bool); cen[keep] = True
                    N = hist2(cell[dat], (xc_all - xi)[dat], ncell)
                    T = template(xc_all[cen], g[cen], cell[cen], xi[prot], g[prot], ng, ncell)
                    o = {"all": eps_from(N, T), "n_prot": int(prot.sum()), "n_cen": int(cen.sum())}
                    if smp in ("z kand. D", "wszystkie") and pname == list(pools)[0] and sname in ("all", "AFP", "nonAFP", "none"):
                        cidx = np.arange(ncell)
                        for nm, (a, bb) in NSG.items():
                            o[nm] = eps_from(N, T, cells=np.flatnonzero((cidx % NCL >= a) & (cidx % NCL < bb)), nboot=100)
                        for k in range(ng):
                            o[f"LB{k}"] = eps_from(N, T, cells=np.flatnonzero(cidx // NCL == k), nboot=100)
                            o[f"LB{k}"]["mu"] = float(d["mu"][dat & (g == k)].mean())
                        Tn = template(xc_all[cen], g[cen], g[cen], xi[prot], g[prot], ng, ng)      # bez warstw nSelBL
                        o["nostrata"] = eps_from(hist2(g[dat], (xc_all - xi)[dat], ng), Tn, nboot=50)
                        for lo, hi in ((0.03, 0.05), (0.05, 0.075), (0.075, 0.1001)):
                            mx = (xi >= lo) & (xi < hi)
                            Tx = template(xc_all[cen], g[cen], cell[cen], xi[prot & mx], g[prot & mx], ng, ncell)
                            o[f"xi {lo:g}-{hi:g}"] = eps_from(hist2(cell[dat & mx], (xc_all - xi)[dat & mx], ncell), Tx, nboot=100)
                    acc[s] = o
                    print(r, key, s, "N", int(o["all"]["N_tot"]), "raw %.4f" % o["all"]["eps_raw"], {nm: "%.4f±%.4f f=%.4f" % (o["all"][nm]["eps_true"], o["all"][nm]["eps_err"], o["all"][nm]["f_acc_SR"]) for nm in NORM}, flush=True)
                R[key] = acc
    res[str(r)] = R
C.jdump(res, os.path.join(C.OUT, "s3.json"))
