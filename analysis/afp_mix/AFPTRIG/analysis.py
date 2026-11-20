#!/usr/bin/env python3
"""AFPTRIG/analysis.py — metoda kształtu na danych: estymator nominalny, warianty, domknięcia w CR, kategorie D.
Wynik: results_raw_<run>.json. Żadnych fitów masy / uzysków sygnału — tylko liczności w oknach i frakcje tła przypadkowego.

Estymator (strona s, komórka c = blok LB × klasa nSelBL):  N_acc(SR) = Σ_c N_CR,c(dane) · R_c,  R_c = T_c(SR)/T_c(CR),
T_c = wszystkie pary (zdarzenie centralne z puli bez protonu po stronie s, proton „niemożliwy” z tego samego bloku LB)."""
import json, sys, time
import numpy as np
import lib, core
from core import DELTA, M_NOM

INF = np.inf
NCLS = len(lib.NSEL_EDGES) - 1
POOL_THR = 0.14


class Mix:
    """Mieszanie pula centralna × protony w blokach: ułamek par z Δ ≤ e dla każdego zdarzenia puli (dokładnie, searchsorted)."""
    def __init__(self, xc_sorted, g_c, xp, g_p):
        self.x = xc_sorted; self.pos = {}; self.prot = {}
        for g in np.unique(g_c):
            self.pos[int(g)] = np.flatnonzero(g_c == g)
            p = np.sort(xp[g_p == g]); self.prot[int(g)] = p
        self.n = len(xc_sorted); self.n_prot = len(xp)

    def cum(self, e):
        """F_i(e) = ułamek protonów bloku z x_i − ξ ≤ e  (⇔ ξ ≥ x_i − e)."""
        if e == INF:
            F = np.ones(self.n, np.float32)
        elif e == -INF:
            F = np.zeros(self.n, np.float32)
        else:
            F = np.empty(self.n, np.float32)
        for g, pos in self.pos.items():
            p = self.prot[g]
            if len(p) == 0:
                F[pos] = np.nan; continue
            if np.isfinite(e):
                F[pos] = 1.0 - np.searchsorted(p, self.x[pos] - e, side="left") / len(p)
        return F

    def frac(self, lo, hi):
        return self.cum(hi) - self.cum(lo)

    def relvar(self, sr, cr):
        """względna wariancja R od liczby protonów per blok (wagi puli = 1)."""
        out = {}
        for g, pos in self.pos.items():
            p = self.prot[g]; xs = self.x[pos]
            if len(p) < 2 or len(xs) == 0:
                out[g] = np.nan; continue
            b = []
            for lo, hi in (sr, cr):
                b.append((np.searchsorted(xs, p + hi, side="right") - np.searchsorted(xs, p + lo, side="right")).astype(float))
            S0, S1 = b[0].sum(), b[1].sum()
            out[g] = float(((b[0] - S0 / S1 * b[1]) ** 2).sum() / S1 ** 2 / (S0 / S1) ** 2) if S0 > 0 and S1 > 0 else np.nan
        return out


class Run:
    def __init__(self, run):
        self.run = run; d = self.d = lib.load(run)
        self.has = d["hasD"]; self.cls = np.clip(np.searchsorted(lib.NSEL_EDGES, d["nsel"], side="right") - 1, 0, NCLS - 1)
        mb = d["T_HLT_noalg_mb_L1MBTS_1"] | d["T_HLT_noalg_mb_L1MBTS_2"]
        self.S = {"afp_any": d["afp"], "pt4_L1AFP": d["T_HLT_mb_sptrk_pt4_L1AFP_A_OR_C"], "noalg_L1AFP": d["T_HLT_noalg_L1AFP_A_OR_C"], "nonAFP": d["non"],
                  "twin": d["T_HLT_mb_sptrk_pt4_L1MBTS_2"], "mb": mb | d["T_HLT_mb_sptrk_pt2_L1MBTS_2"] | d["T_HLT_mb_mbts_L1MBTS_1"] | d["T_HLT_mb_sptrk_L1MBTS_1"],
                  "j12": d["T_HLT_noalg_L1J12"]}
        self.S["mbts"] = self.S["twin"] | self.S["mb"]
        self.delta = {s: d[f"xical_{s}"] - d[f"xi1_{s}"] for s in lib.SIDES}
        self.cache = {}

    def groups(self, s, stream, nlb, var):
        k = ("g", s, stream, nlb, var)
        if k not in self.cache:
            ref = self.S[stream] & self.d[f"tag_{s}"] & self.has
            x = self.d["lb"] if var == "lb" else self.d["mu"]
            if nlb <= 1:
                self.cache[k] = np.zeros(len(x), np.int32)
            else:
                e = np.unique(np.quantile(x[ref], np.linspace(0, 1, nlb + 1)[1:-1]))
                self.cache[k] = np.searchsorted(e, x, side="right").astype(np.int32)
        return self.cache[k]

    def protons(self, s, stream, kind):
        d = self.d; t = d[f"tag_{s}"]
        if kind == "impossible":
            return self.S[stream] & t & (d[f"xical_{s}"] > POOL_THR)
        if kind == "impossible_hasD":
            return self.S[stream] & t & (d[f"xical_{s}"] > POOL_THR) & self.has
        if kind == "alltags":
            return self.S[stream] & t
        if kind == "nonafp_veto":
            return d["non"] & t & (d[f"xical_{s}"] > POOL_THR)
        raise KeyError(kind)

    def pool_idx(self, s, stream, pool):
        k = ("i", s, stream, pool)
        if k not in self.cache:
            d = self.d
            pm = (self.S[stream] & d[f"tag_{s}"] & self.has) if pool == "self" else (self.S[pool] & (d[f"np_{s}"] == 0) & self.has)
            i = np.flatnonzero(pm)
            self.cache[k] = i[np.argsort(d[f"xical_{s}"][i], kind="stable")]
        return self.cache[k]

    def mix(self, s, stream, pool, prot, nlb, var):
        k = ("m", s, stream, pool, prot, nlb, var)
        if k not in self.cache:
            d = self.d; g = self.groups(s, stream, nlb, var); idx = self.pool_idx(s, stream, pool); qm = self.protons(s, stream, prot)
            self.cache[k] = Mix(d[f"xical_{s}"][idx], g[idx], d[f"xi1_{s}"][qm], g[qm])
        return self.cache[k]

    def est(self, s, stream, pool, prot="impossible", nlb=8, var="lb", strat=True, m=M_NOM, sr=(-INF, DELTA), cr=None, w=None):
        """Jeden estymator. w = wagi per zdarzenie (None → 1 dla zdarzeń z D). → dict."""
        d = self.d
        cr = (DELTA + m, INF) if cr is None else cr
        M = self.mix(s, stream, pool, prot, nlb, var); idx = self.pool_idx(s, stream, pool); g = self.groups(s, stream, nlb, var)
        ka = ("a", s, stream, pool, prot, nlb, var, sr, cr)
        if ka not in self.cache:
            self.cache[ka] = (M.frac(*sr), M.frac(*cr), M.relvar(sr, cr))
        a_sr, a_cr, prv = self.cache[ka]
        W = self.has.astype(float) if w is None else w
        ncell = (int(g.max()) + 1) * NCLS
        cell2g = np.arange(ncell) // NCLS
        di = np.flatnonzero(self.S[stream] & d[f"tag_{s}"] & (W > 0)); dl = self.delta[s][di]
        dcell = g[di] * NCLS + (self.cls[di] if strat else 0)
        in_cr = (dl > cr[0]) & (dl <= cr[1]); in_sr = (dl > sr[0]) & (dl <= sr[1])
        wc = W[idx]; kc = wc > 0; ic = idx[kc]
        ccell = g[ic] * NCLS + (self.cls[ic] if strat else 0)
        r = core.shape_estimate(dcell, W[di], in_cr, ccell, wc[kc], a_sr[kc].astype(float), a_cr[kc].astype(float), cell2g=cell2g, prv=prv, ncell=ncell)
        r.update(n_sr=float(W[di][in_sr].sum()), n_sr_events=int(in_sr.sum()), n_tag=float(W[di].sum()), n_pool=int(kc.sum()), n_protons=M.n_prot)
        r["f_acc"] = r["n_acc"] / r["n_sr"] if r["n_sr"] > 0 else None
        r["f_acc_err_stat"] = r["err"] / r["n_sr"] if r["n_sr"] > 0 else None
        return r

    def drop(self, prefix):
        for k in [k for k in self.cache if k[0] == prefix]:
            del self.cache[k]

    def hist_pred(self, s, stream, pool, edges, m=M_NOM, nlb=8, strat=True):
        """Przewidywany histogram Δ tła przypadkowego (zdarzenia z D), normalizacja w CR nominalnym per komórka."""
        d = self.d; M = self.mix(s, stream, pool, "impossible", nlb, "lb"); idx = self.pool_idx(s, stream, pool); g = self.groups(s, stream, nlb, "lb")
        ncell = (int(g.max()) + 1) * NCLS
        ccell = g[idx] * NCLS + (self.cls[idx] if strat else 0)
        di = np.flatnonzero(self.S[stream] & d[f"tag_{s}"] & self.has); dl = self.delta[s][di]
        dcell = g[di] * NCLS + (self.cls[di] if strat else 0)
        N = np.bincount(dcell[dl > DELTA + m], minlength=ncell).astype(float)
        Scr = np.bincount(ccell, 1.0 - M.cum(DELTA + m), ncell)
        k = np.where(Scr > 0, N / np.where(Scr > 0, Scr, 1), 0.0)
        C = np.array([(np.bincount(ccell, M.cum(e), ncell) * k).sum() for e in edges])
        return np.diff(C)


def nominal_pool(run, stream):
    if stream == "nonAFP":
        return "nonAFP"
    if run == 427929:
        return "twin" if stream == "pt4_L1AFP" else "nonAFP"
    return {"pt4_L1AFP": "twin", "noalg_L1AFP": "mb", "afp_any": "mbts"}[stream]


KEYS = ("n_sr", "n_cr", "n_acc", "err", "f_acc", "f_acc_err_stat")


def main():
    t0 = time.time(); OUT = {}
    runs = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else lib.RUNS
    for run in runs:
        R = Run(run); O = OUT[str(run)] = {}
        for stream in ("afp_any", "pt4_L1AFP", "noalg_L1AFP", "nonAFP"):
            np_ = nominal_pool(run, stream); Os = O[stream] = {"nominal_pool": np_}
            pools = ["nonAFP", "twin", "self"] if run == 427929 else ["nonAFP", "mbts", "twin", "mb", "j12", "self"]
            for s in lib.SIDES:
                o = Os[s] = {}
                o["nominal"] = R.est(s, stream, np_)
                o["pool"] = {p: R.est(s, stream, p) for p in pools}
                o["pool_nostrat"] = {p: R.est(s, stream, p, strat=False) for p in pools}
                o["protons"] = {k: R.est(s, stream, np_, prot=k) for k in ("impossible_hasD", "alltags", "nonafp_veto")}
                o["strat"] = {"nlb1": R.est(s, stream, np_, nlb=1), "nlb4": R.est(s, stream, np_, nlb=4), "nlb16": R.est(s, stream, np_, nlb=16),
                              "mu8": R.est(s, stream, np_, var="mu"), "no_nsel": R.est(s, stream, np_, strat=False)}
                o["margin"] = {f"{m:g}": R.est(s, stream, np_, m=m) for m in (0.0, 0.01, 0.025, 0.04, 0.06, 0.08, 0.12)}
                cl = {}
                for nm, cr, sr in (("norm>0.10_pred_0.045-0.10", (0.10, INF), (0.045, 0.10)), ("norm>0.10_pred_0.045-0.07", (0.10, INF), (0.045, 0.07)),
                                   ("norm>0.10_pred_0.07-0.10", (0.10, INF), (0.07, 0.10)), ("norm>0.20_pred_0.10-0.20", (0.20, INF), (0.10, 0.20)),
                                   ("norm>0.20_pred_0.045-0.10", (0.20, INF), (0.045, 0.10)), ("norm_0.045-0.10_pred_>0.10", (0.045, 0.10), (0.10, INF)),
                                   ("norm>0.10_pred_0.03-0.045", (0.10, INF), (0.03, 0.045)), ("norm>0.10_pred_0.02-0.03", (0.10, INF), (0.02, 0.03)),
                                   ("norm>0.045_pred_0.03-0.045", (0.045, INF), (0.03, 0.045)), ("norm>0.045_pred_0.02-0.03", (0.045, INF), (0.02, 0.03))):
                    cl[nm] = {p: R.est(s, stream, p, sr=sr, cr=cr) for p in dict.fromkeys([np_, "nonAFP"])}
                    cl[nm]["nostrat"] = R.est(s, stream, np_, sr=sr, cr=cr, strat=False)
                o["closure"] = cl
                dep = {}
                for nm, x, q in (("lb", R.d["lb"], 6), ("mu", R.d["mu"], 6)):
                    ref = R.S[stream] & R.d[f"tag_{s}"] & R.has
                    e = np.unique(np.quantile(x[ref], np.linspace(0, 1, q + 1))).astype(float); e[-1] += 1e-6
                    rows = []
                    for i in range(len(e) - 1):
                        w = (R.has & (x >= e[i]) & (x < e[i + 1])).astype(float)
                        rr = R.est(s, stream, np_, w=w); rows.append({"lo": float(e[i]), "hi": float(e[i + 1]), **{k: rr[k] for k in KEYS}})
                    dep[nm] = rows
                rows = []
                for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 50), (50, 10000)):
                    w = (R.has & (R.d["nsel"] >= lo) & (R.d["nsel"] < hi)).astype(float)
                    rr = R.est(s, stream, np_, w=w); rows.append({"lo": lo, "hi": hi, **{k: rr[k] for k in KEYS}})
                dep["nsel"] = rows
                o["dep"] = dep
                if stream in ("afp_any", "pt4_L1AFP"):
                    cats = {}
                    for i, c in enumerate(lib.CATS):
                        w = R.d["cnt"][:, i].astype(float)
                        cats[c] = {"nominal": R.est(s, stream, np_, w=w), "alt_nonAFP": R.est(s, stream, "nonAFP", w=w), "nostrat": R.est(s, stream, np_, strat=False, w=w)}
                    o["cats"] = cats
                ed = np.round(np.concatenate([np.arange(-0.10, 0.3001, 0.01), [0.35, 0.4, 0.5, 0.7, 1.0, 3.0]]), 6)
                ed = np.sort(np.append(ed, 0.045))
                dm = R.S[stream] & R.d[f"tag_{s}"] & R.has
                o["hist_delta"] = {"edges": ed.tolist(), "data": np.histogram(R.delta[s][dm], ed)[0].tolist(),
                                   "pred": {p: R.hist_pred(s, stream, p, ed).tolist() for p in dict.fromkeys([np_, "nonAFP"])},
                                   "pred_nostrat": R.hist_pred(s, stream, np_, ed, strat=False).tolist()}
                print(run, stream, s, "nominal", {k: round(v, 4) if isinstance(v, float) else v for k, v in o["nominal"].items()}, f"[{time.time() - t0:.0f}s]", flush=True)
                R.drop("a"); R.drop("m")
            R.cache.clear()
            json.dump(OUT, open("results_raw_%s.json" % "_".join(str(r) for r in runs), "w"), indent=1)


if __name__ == "__main__":
    main()
