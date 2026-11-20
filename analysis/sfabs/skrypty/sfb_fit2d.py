#!/usr/bin/env python3
"""SFABS metoda B: tag D*, fit 2D (m_Kpi, dm) laczny pass/fail WP70 z eps jako parametrem. Kryteria: ../KRYTERIA.md par. B.
  sfb_fit2d.py mc    CACHE OUT     domkniecie G_B1 (mb_nd, pelne tlo) + test samego sygnalu
  sfb_fit2d.py data  CACHE OUT     dane (oba runy, 3 biny): fit nominalny, warianty, bootstrap
  sfb_fit2d.py toys  CACHE OUT DATAJSON   G_B2 (sygnal fd_nd w tle z pasm masy danych)
Skladowe: S (D* prawdziwe), R (D0 prawdziwe x przypadkowy pion), C (kombinatoryka), P (pik w dm, gladkie w m)."""
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

M0, MPI, DM0 = 1.86484, 0.13957, 0.14543
ME = np.linspace(1.79, 1.94, 16)                 # 10 MeV
DE = np.linspace(0.1396, 0.1560, 83)             # 0,2 MeV (82 biny)
THR = {"3-5": 0.822, "5-8": 0.8719, "8-30": 0.8631}
PTB = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
RUNS = {"427929": ("data_427929", lambda c: (c & 3) == 3), "428770": ("data_428770", lambda c: (c & 1) == 1)}
GC_SC = 0.10
CACHE = None


class Cdf:
    """Dystrybuanta empiryczna na gestej siatce (PCHIP); bins(edges, x0, shift, scale)."""

    def __init__(self, x, lo, hi, n=600, w=None):
        from scipy.interpolate import PchipInterpolator
        g = np.linspace(lo, hi, n + 1)
        h = np.histogram(x, g, weights=w)[0].astype(float)
        self.g = g
        self.f = PchipInterpolator(g, np.concatenate([[0.0], np.cumsum(h) / max(h.sum(), 1e-300)]))

    def bins(self, edges, x0, d=0.0, s=1.0):
        return np.diff(self.f(np.clip((edges - x0 - d) / s + x0, self.g[0], self.g[-1])))


def load(s):
    z = np.load(os.path.join(CACHE, f"pairs_{s}.npz"))
    return {k: z[k] for k in z.files}


def sel(T, bn, extra=None):
    a, b = PTB[bn]
    w = (T["pt"] >= a) & (T["pt"] < b) & (T["m"] >= ME[0]) & (T["m"] < ME[-1]) & (T["ds_dm"] >= DE[0]) & (T["ds_dm"] < DE[-1])
    if extra is not None:
        w &= extra
    return w


def is_S(T):
    return np.isin(T["truth"], (0, 1)) & np.isin(T["ds_truth"], (0, 1)) & (T["swap"] == 0)


def is_R(T):
    return np.isin(T["truth"], (0, 1)) & ~np.isin(T["ds_truth"], (0, 1)) & (T["swap"] == 0)


def hist2(T, w, bn, wt=None):
    """-> [2 (pass, fail), nm, nd]"""
    p = T["score"][w] >= THR[bn]
    o = []
    for k in (p, ~p):
        o.append(np.histogram2d(T["m"][w][k], T["ds_dm"][w][k], [ME, DE], weights=None if wt is None else wt[w][k])[0])
    return np.array(o)


def templates(fd, bn):
    w = sel(fd, bn, is_S(fd))
    p = fd["score"][w] >= THR[bn]
    m, dm = fd["m"][w], fd["ds_dm"][w]
    return {"m": [Cdf(m[p], 1.74, 1.99), Cdf(m[~p], 1.74, 1.99)], "m_all": Cdf(m, 1.74, 1.99), "dm": Cdf(dm, 0.1380, 0.1580, 2000),
            "eps_mc": float(p.mean()), "n_mc": int(w.sum())}


def thr_shape(p, c):
    """Prog (dm - mpi)^p exp(-c (dm - mpi)) calkowany po binach (4 podpunkty), znormalizowany do 1 w zakresie."""
    k = 4
    x = (DE[:-1, None] + (np.arange(k)[None, :] + 0.5) / k * np.diff(DE)[:, None]) - MPI
    v = (np.clip(x, 1e-9, None) ** p * np.exp(-c * x)).mean(1)
    return v / v.sum()


def fit2d(H, tp, opt=None, start=None, hesse=True, minos=False):
    """H [2, nm, nd]. opt: P (bool), expthr (bool), common (bool), mcut (int: bez tylu skrajnych binow m z kazdej strony), mord (1|2)."""
    from iminuit import Minuit
    o = dict(P=True, expthr=True, common=False, mcut=0, mord=1, tieR=True, kR=1.0)
    o.update(opt or {})
    mc = o["mcut"]
    ms = slice(mc, len(ME) - 1 - mc)
    Hf = H[:, ms, :]
    xm = (0.5 * (ME[1:] + ME[:-1]) - M0) / 0.075
    names = ["dmu", "sc", "ddm", "dsc", "N", "eps"]
    p0 = [0.0, 1.0, 0.0, 1.0, 0.0, tp["eps_mc"]]
    lim = [(-0.006, 0.006), (0.7, 1.5), (-0.0004, 0.0004), (0.6, 2.5), (0, None), (0.0, 1.5)]
    ntot = H.sum(axis=(1, 2))
    wm = (np.abs(0.5 * (ME[1:] + ME[:-1]) - M0) < 0.03)
    wd = (0.5 * (DE[1:] + DE[:-1]) > 0.1445) & (0.5 * (DE[1:] + DE[:-1]) < 0.1465)
    wdo = (0.5 * (DE[1:] + DE[:-1]) > 0.1480) & (0.5 * (DE[1:] + DE[:-1]) < 0.1520)
    s_est = []
    for c in range(2):
        s_est.append(max(H[c][wm][:, wd].sum() - H[c][wm][:, wdo].sum() * wd.sum() / wdo.sum() * 0.93, 0.0) * 1.35)
    p0[4] = max(sum(s_est), 10.0)
    if s_est[0] + s_est[1] > 0:
        p0[5] = min(max(s_est[0] / (s_est[0] + s_est[1]), 0.05), 1.2)
    lim[4] = (0.0, 3 * p0[4] + 50 * math.sqrt(ntot.sum()))
    for c in "pf":
        nt = ntot[0 if c == "p" else 1]
        names += [f"NR{c}", f"NC{c}", f"a{c}", f"a2{c}", f"p{c}", f"c{c}", f"NP{c}"]
        p0 += [0.02 * nt, 0.95 * nt, 0.0, 0.0, 0.5, 0.0, 0.0]
        lim += [(0, nt + 50), (0, 1.5 * nt + 50), (-5, 5), (-5, 5), (0.05, 3.0), (-300, 300), (0, 0.5 * nt + 50)]
    if start is not None and len(start) == len(p0):
        p0 = list(start)

    def comps(q):
        dmu, sc, ddm, dsc, N, eps = q[:6]
        sd = tp["dm"].bins(DE, DM0, ddm, dsc)
        out = []
        for c in range(2):
            r = q[6 + 7 * c: 13 + 7 * c]
            T = tp["m_all"] if o["common"] else tp["m"][c]
            sm = T.bins(ME, M0, dmu, sc)
            th = thr_shape(r[4], r[5])
            bm = np.exp(np.clip(r[2] * xm + r[3] * (2 * xm * xm - 1), -30, 30))
            bm = bm / bm[ms].sum()
            ns = N * (eps if c == 0 else 1.0 - eps)
            nr = (q[6] * (min(o["kR"] * eps, 1.0) if c == 0 else 1.0 - min(o["kR"] * eps, 1.0))) if o["tieR"] else r[0]      # tieR: eps_R = eps (NRp = N_R calkowite)
            out.append((ns * np.outer(sm, sd), nr * np.outer(sm, th), r[1] * np.outer(bm, th), r[6] * np.outer(bm, sd)))
        return out

    def model(q):
        return np.maximum(np.array([sum(cc) for cc in comps(q)])[:, ms, :], 1e-9)

    def cost(q):
        mu = model(q)
        v = 2.0 * float(np.sum(mu - Hf) + np.sum(np.where(Hf > 0, Hf * np.log(np.where(Hf > 0, Hf, 1.0) / mu), 0.0)))
        return v + ((q[1] - 1.0) / GC_SC) ** 2
    m = Minuit(cost, p0, name=names)
    m.errordef = 1.0
    m.strategy = 1
    for n_, l_ in zip(names, lim):
        m.limits[n_] = l_
    for c in "pf":
        if o["mord"] < 2:
            m.fixed[f"a2{c}"] = True
        if not o["expthr"]:
            m.fixed[f"c{c}"] = True
        if not o["P"]:
            m.fixed[f"NP{c}"] = True
    if o["tieR"]:
        m.values["NRf"] = 0.0; m.fixed["NRf"] = True
        m.limits["NRp"] = (0, 0.2 * ntot.sum() + 100)
    for n_ in ("dmu", "sc", "ddm", "dsc"):
        m.fixed[n_] = True
    m.migrad(ncall=30000)
    for n_ in ("dmu", "sc", "ddm", "dsc"):
        m.fixed[n_] = False
    m.migrad(ncall=30000)
    for it in range(2):
        if m.valid:
            break
        m.simplex(ncall=30000)
        m.migrad(ncall=60000)
    v = np.array(m.values)
    res = {"x": v.tolist(), "N": float(v[4]), "eps": float(v[5]), "valid": bool(m.valid), "fval": float(m.fval), "dmu": float(v[0]), "sc": float(v[1]), "ddm": float(v[2]), "dsc": float(v[3]),
           "NR": [float(v[6]), float(v[13])], "NC": [float(v[7]), float(v[14])], "NP": [float(v[12]), float(v[19])], "pthr": [float(v[10]), float(v[17])]}
    if hesse:
        try:
            m.hesse()
            res["eps_err"] = float(m.errors["eps"]); res["N_err"] = float(m.errors["N"])
            res["acc"] = bool(m.fmin.has_accurate_covar)
        except Exception:  # noqa: BLE001
            res["eps_err"] = None
    if minos:
        try:
            m.minos("eps")
            res["eps_minos"] = [float(m.merrors["eps"].lower), float(m.merrors["eps"].upper)]
        except Exception:  # noqa: BLE001
            pass
    mu = model(v)
    # chi2 po zgrubieniu (dm x4; przy malej statystyce: m x3, dm x8), tylko komorki z mu >= 5
    from scipy.stats import chi2 as C2
    for km, kd in ((1, 4), (3, 8), (5, 16)):
        nd = (len(DE) - 1) // kd * kd; nm = Hf.shape[1] // km * km
        Hc = Hf[:, :nm, :nd].reshape(2, nm // km, km, nd // kd, kd).sum((2, 4)); Mc = mu[:, :nm, :nd].reshape(2, nm // km, km, nd // kd, kd).sum((2, 4))
        ok = Mc >= 5
        res["chi2"] = float((((Hc - Mc) ** 2) / Mc)[ok].sum()); res["ndf"] = int(ok.sum() - m.nfit); res["gof_grid"] = [km, kd]
        if res["ndf"] >= 15:
            break
    res["p"] = float(C2.sf(res["chi2"], res["ndf"])) if res["ndf"] > 0 else None
    res["at_limit"] = bool(v[5] < 0.01 or v[5] > 1.49 or v[3] < 0.61 or v[3] > 2.49 or abs(v[0]) > 0.0059 or v[1] < 0.71 or v[1] > 1.49)
    res["_model"] = model
    return res


def clean(r):
    return {k: v for k, v in r.items() if not k.startswith("_")}


VARS = {"freeR": dict(tieR=False), "kR1.05": dict(kR=1.05), "kR0.95": dict(kR=0.95), "noP": dict(P=False), "noexp": dict(expthr=False), "common": dict(common=True), "mcut1": dict(mcut=1), "mord2": dict(mord=2)}


# ------------------------------------------------------------------ mc
def cmd_mc(out):
    fd, mb = load("fd_nd"), load("mb_nd")
    R = {}
    for bn in PTB:
        tp = templates(fd, bn)
        w = sel(mb, bn)
        H = hist2(mb, w, bn)
        ws = sel(mb, bn, is_S(mb)); tr_p = int((mb["score"][ws] >= THR[bn]).sum()); tr_n = int(ws.sum())
        wr = sel(mb, bn, is_R(mb))
        f = fit2d(H, tp, minos=True)
        r = {"eps_mc_fd": tp["eps_mc"], "truth": {"N": tr_n, "eps": tr_p / tr_n, "eps_err": math.sqrt(tr_p * (tr_n - tr_p) / tr_n ** 3), "NR": int(wr.sum()),
                                                  "NR_pass": int((mb["score"][wr] >= THR[bn]).sum())}, "ntot": [float(H[0].sum()), float(H[1].sum())], "fit": clean(f)}
        r["var"] = {k: clean(fit2d(H, tp, opt=v)) for k, v in VARS.items()}
        # sam sygnal fd (test faktoryzacji szablonu): polowa fd jako pseudo-dane + tlo mb bez S
        e = f["eps"]; s = f.get("eps_err") or float("nan")
        print(f"[mc] {bn}: prawda eps {r['truth']['eps']:.4f} (N {tr_n}, R {int(wr.sum())}); fit eps {e:.4f} +- {s:.4f} N {f['N']:.0f} +- {f.get('N_err', 0):.0f} NR {f['NR']} NP {f['NP']} p(chi2) {f['p']:.3f} "
              f"pull {(e - r['truth']['eps']) / s:+.2f}; warianty " + " ".join(f"{k}:{v['eps']:.3f}" for k, v in r["var"].items()), flush=True)
        # bootstrap Poissona po binach: rozrzut vs HESSE
        rng = np.random.default_rng(7)
        bs = []
        for i in range(20):
            fb = fit2d(rng.poisson(H).astype(float), tp, start=f["x"], hesse=False)
            if fb["valid"]:
                bs.append(fb["eps"])
        r["boot_rms"] = float(np.std(bs, ddof=1)); r["boot_n"] = len(bs)
        print(f"      rozrzut Poissona eps: {r['boot_rms']:.4f} (n {len(bs)})", flush=True)
        R[bn] = r
    json.dump(R, open(out, "w"), indent=1)


# ------------------------------------------------------------------ data
def _boot_task(a):
    H, tp, start = a
    f = fit2d(H, tp, start=start, hesse=False)
    return f["eps"] if f["valid"] else None


def cmd_data(out, nboot=30):
    fd = load("fd_nd")
    R = {}
    for run, (smp, cf) in RUNS.items():
        T = load(smp)
        okc = cf(T["cls"])
        ev = T["uid"] >> np.uint64(8)                      # przyblizenie zdarzenia: uid kandydata bez dolnych bitow (bootstrap po kandydatach-grupach)
        for bn in PTB:
            tp = templates(fd, bn)
            w = sel(T, bn, okc)
            H = hist2(T, w, bn)
            f = fit2d(H, tp, minos=True)
            r = {"eps_mc": tp["eps_mc"], "ntot": [float(H[0].sum()), float(H[1].sum())], "fit": clean(f)}
            r["var"] = {k: clean(fit2d(H, tp, opt=v)) for k, v in VARS.items()}
            rng = np.random.default_rng(11)
            uq, inv = np.unique(ev[w], return_inverse=True)
            tasks = []
            for i in range(nboot):
                wt = np.zeros(len(T["m"])); wt[w] = rng.poisson(1.0, len(uq))[inv]
                tasks.append((hist2(T, w, bn, wt), tp, f["x"]))
            with Pool(6) as p:
                bs = [b for b in p.map(_boot_task, tasks) if b is not None]
            r["boot_rms"] = float(np.std(bs, ddof=1)) if len(bs) > 3 else None; r["boot_n"] = len(bs)
            e, s = f["eps"], f.get("eps_err") or float("nan")
            r["SF"] = e / tp["eps_mc"]; r["SF_hesse"] = s / tp["eps_mc"]; r["SF_boot"] = (r["boot_rms"] or float("nan")) / tp["eps_mc"]
            vmax = max(abs(v["eps"] - e) for k_, v in r["var"].items() if v["valid"] and k_ != "freeR") / tp["eps_mc"]
            r["SF_model"] = vmax
            print(f"[data] {run} {bn}: N {f['N']:.0f} +- {f.get('N_err', 0):.0f} eps {e:.4f} +- {s:.4f} (boot {r['boot_rms']}, minos {f.get('eps_minos')}) SF {r['SF']:.3f} +- {r['SF_hesse']:.3f} model {vmax:.3f} "
                  f"NR {[round(x) for x in f['NR']]} NP {[round(x) for x in f['NP']]} NC {[round(x) for x in f['NC']]} sc {f['sc']:.3f} dsc {f['dsc']:.3f} p {f['p']:.3f} lim {f['at_limit']} | "
                  + " ".join(f"{k}:{v['eps'] / tp['eps_mc']:.3f}" for k, v in r["var"].items()), flush=True)
            R[f"{run}|{bn}"] = r
    json.dump(R, open(out, "w"), indent=1)


# ------------------------------------------------------------------ toys (G_B2)
def _toy_task(a):
    H, tp = a
    f = fit2d(H, tp)
    return (f["eps"], f.get("eps_err"), f["N"], f["valid"])


def cmd_toys(out, datajson, ntoy=30, only=None):
    D = json.load(open(datajson))
    fd, mb = load("fd_nd"), load("mb_nd")
    R = {}
    mc_ctr = 0.5 * (ME[1:] + ME[:-1])
    for run, (smp, cf) in RUNS.items():
        T = load(smp)
        okc = cf(T["cls"])
        for bn in PTB:
            key = f"{run}|{bn}"
            if only and only not in key:
                continue
            d = D[key]["fit"]
            tp = templates(fd, bn)
            w = sel(T, bn, okc)
            H = hist2(T, w, bn)
            nS = d["N"]
            nB = [H[c].sum() - nS * (d["eps"] if c == 0 else 1 - d["eps"]) for c in range(2)]
            # tlo z pasm masy danych: rozklad (dm | kategoria) z |m - m0| > 45 MeV; masa: exp(liniowy) dopasowany do pasm (per kategoria)
            sb = np.abs(mc_ctr - M0) > 0.045
            B = []
            for c in range(2):
                hd = H[c][sb].sum(0); hm = H[c][:, :].sum(1)
                x = (mc_ctr - M0) / 0.075
                co = np.polyfit(x[sb], np.log(np.maximum(hm[sb], 1.0)), 1)
                sm = np.exp(np.polyval(co, x)); sm /= sm.sum()
                B.append(np.outer(sm, hd / hd.sum()))
            # R z mb_nd (prawda): ksztalt 2D, licznosc = NR_MC/NS_MC x nS, podzial pass/fail z prawdy
            wr = sel(mb, bn, is_R(mb)); ws = sel(mb, bn, is_S(mb))
            HR = hist2(mb, wr, bn); fR = wr.sum() / ws.sum()
            # sygnal: fd_nd prawda, losowanie zdarzen
            wsf = np.nonzero(sel(fd, bn, is_S(fd)))[0]
            pf = fd["score"][wsf] >= THR[bn]
            res = {}
            for lab, sf_inj in (("eps_MC", 1.0), ("0.85eps_MC", 0.85)):
                e_true = tp["eps_mc"] * sf_inj
                rng = np.random.default_rng(100 + int(sf_inj * 100))
                tasks = []
                for i in range(ntoy):
                    n = rng.poisson(nS)
                    npass = rng.binomial(n, e_true)
                    ip = rng.choice(wsf[pf], npass); iF = rng.choice(wsf[~pf], n - npass)
                    Hs = np.array([np.histogram2d(fd["m"][ix], fd["ds_dm"][ix], [ME, DE])[0] for ix in (ip, iF)])
                    Ht = Hs.copy()
                    for c in range(2):
                        nr = fR * nS * HR[c].sum() / max(HR.sum(), 1)
                        Ht[c] += rng.poisson(HR[c] / max(HR[c].sum(), 1) * nr)
                        Ht[c] += rng.poisson(B[c] * max(nB[c] - nr, 0))
                    tasks.append((Ht, tp))
                with Pool(6) as p:
                    o = p.map(_toy_task, tasks)
                e = np.array([x[0] for x in o if x[3] and x[1]]); s = np.array([x[1] for x in o if x[3] and x[1]])
                pull = (e - e_true) / s
                res[lab] = {"eps_true": e_true, "n": len(e), "mean": float(e.mean()), "rms": float(e.std(ddof=1)), "mean_hesse": float(s.mean()), "pull_mean": float(pull.mean()), "pull_rms": float(pull.std(ddof=1)),
                            "bias": float(e.mean() - e_true), "bias_err": float(e.std(ddof=1) / math.sqrt(len(e))), "N_mean": float(np.mean([x[2] for x in o if x[3] and x[1]]))}
                q = res[lab]
                print(f"[toys] {key} {lab}: eps_true {e_true:.3f} <fit> {q['mean']:.3f} rms {q['rms']:.3f} <hesse> {q['mean_hesse']:.3f} bias {q['bias']:+.4f} +- {q['bias_err']:.4f} pull {q['pull_mean']:+.2f} / {q['pull_rms']:.2f} "
                      f"N {q['N_mean']:.0f} (wstrz. {nS:.0f}) n {q['n']}", flush=True)
            res["nS"] = nS; res["nB"] = [float(x) for x in nB]
            R[key] = res
    json.dump(R, open(out, "w"), indent=1)


if __name__ == "__main__":
    cmd, CACHE, out = sys.argv[1:4]
    if cmd == "mc":
        cmd_mc(out)
    elif cmd == "data":
        cmd_data(out)
    elif cmd == "toys":
        cmd_toys(out, sys.argv[4], only=(sys.argv[5] if len(sys.argv) > 5 else None))
