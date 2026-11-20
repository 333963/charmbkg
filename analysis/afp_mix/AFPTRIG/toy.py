#!/usr/bin/env python3
"""AFPTRIG/toy.py — test metody kształtu na toy (generator charmbkg.afp_mix_toy.generate, bez zmian) z dodanym tu
trygerem wymagającym protonu: L1 AFP odpala od każdego widzianego protonu z wydajnością ε_A = 0,9, ε_C = 0,45 (jak w danych),
opcjonalnie zależną od ξ; strumień AFP = zdarzenia z odpalonym L1; strumień „MB” = wszystkie zdarzenia (jak 427929)
albo podpróbka obciążona aktywnością (klasa trygera 1 ≈ J12). Prawda: tag przypadkowy = dokładnie jeden proton po stronie,
w oknie, przechodzi cięcie, i NIE jest to własny proton zdarzenia. → toy_results.json"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.expanduser("~/Projects/charm/charm-bkg"))
from charmbkg import afp_mix_toy as T
import lib, core
from analysis import Run, NCLS

SC = {
    "nominal": {},
    "halo_xi_tilt_low (inne widmo ξ halo przy ξ_cal < 0,14)": {"halo_xi_low": (0.14, 0.9)},
    "xical_step_0.12_x1.5 (częstość ×1,5 przy ξ_cal < 0,12)": {"halo_xical_step": (0.12, 1.5)},
    "xical_step_0.07_x2 (częstość ×2 przy ξ_cal < 0,07)": {"halo_xical_step": (0.07, 2.0)},
    "rate_vs_nsel (częstość zależna od nSelBL)": {"halo_nsel": 3.0},
    "xical_smear_0.5 (przeciek własnych protonów)": {"xical_smear": 0.5},
    "trig_eff_vs_xi (ε_L1 rośnie z ξ: ×(0,6 + 4ξ))": {"_trig_slope": True},
    "halo_x3 (większa częstość przypadkowych)": {"halo0": 0.06},
}


def make(seed, kw, n_lb=40, n_per_lb=75000):
    kw = dict(kw); slope = kw.pop("_trig_slope", False)
    tc = T.ToyConfig(n_runs=1, n_lb=n_lb, n_per_lb=n_per_lb, seed=seed, **kw)
    _, ev, tr = T.generate(tc)
    rng = np.random.default_rng(1000 + seed)
    n = len(ev); d = {}
    d["lb"] = ev["lumiBlock"].to_numpy().astype(np.int32); d["mu"] = ev["actual_mu"].to_numpy().astype(np.float32)
    d["nsel"] = ev["nSelBL"].to_numpy().astype(np.int32); d["hasD"] = ev["has_D"].to_numpy().astype(bool)
    fire = np.zeros(n, bool); acc = {}
    for s, eff in (("A", 0.9), ("C", 0.45)):
        npn = ev[f"afp_nProtons_{s}"].to_numpy(); x1 = ev[f"afp_xi1_{s}"].to_numpy()
        d[f"np_{s}"] = npn.astype(np.int16); d[f"xi1_{s}"] = x1.astype(float); d[f"xical_{s}"] = ev[f"xi_cal_{s}"].to_numpy().astype(float)
        d[f"tag_{s}"] = (npn == 1) & (x1 >= lib.XI_LO) & (x1 <= lib.XI_HI)
        e = eff * (np.clip(0.6 + 4 * x1, 0, 1 / eff) if slope else 1.0)
        # każdy proton po stronie odpala niezależnie; przy n ≥ 2 przybliżenie: 1 − (1 − ε)^n z ε najmniejszego ξ
        fire |= (npn > 0) & (rng.random(n) < 1 - (1 - np.clip(e, 0, 1)) ** np.maximum(npn, 1))
        acc[s] = d[f"tag_{s}"] & (d[f"xical_{s}"] <= x1 + lib.DELTA) & ~tr[f"own_tag_{s}"].to_numpy()
    d["afp"] = fire; d["non"] = np.ones(n, bool)
    R = Run.__new__(Run); R.run = 0; R.d = d; R.has = d["hasD"]
    R.cls = np.clip(np.searchsorted(lib.NSEL_EDGES, d["nsel"], side="right") - 1, 0, NCLS - 1)
    j12 = ev["trig_class"].to_numpy() == 1
    R.S = {"afp_any": fire, "nonAFP": d["non"], "j12": j12}
    R.delta = {s: d[f"xical_{s}"] - d[f"xi1_{s}"] for s in lib.SIDES}; R.cache = {}
    return R, acc


def one(seed, kw):
    R, acc = make(seed, kw); out = {}
    for s in lib.SIDES:
        true = float((R.S["afp_any"] & R.has & acc[s]).sum())
        r = {"true": true}
        for nm, a in (("nominal", dict(pool="nonAFP")), ("no_nsel", dict(pool="nonAFP", strat=False)), ("pool_J12_strat", dict(pool="j12")),
                      ("pool_J12_nostrat", dict(pool="j12", strat=False)), ("pool_self", dict(pool="self")), ("prot_alltags", dict(pool="nonAFP", prot="alltags")),
                      ("m0.01", dict(pool="nonAFP", m=0.01)), ("m0.06", dict(pool="nonAFP", m=0.06))):
            e = R.est(s, "afp_any", **a); r[nm] = [e["n_acc"], e["err"]]
            if nm == "nominal":
                r["n_sr"] = e["n_sr"]; r["n_cr"] = e["n_cr"]
        c = R.est(s, "afp_any", "nonAFP", sr=(0.045, 0.10), cr=(0.10, np.inf)); r["closure"] = [c["n_sr"], c["n_acc"], c["err"]]
        out[s] = r
    return out


if __name__ == "__main__":
    seeds = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "1,2,3,4,5").split(",")]
    OUT = {}; t0 = time.time()
    for sc, kw in SC.items():
        OUT[sc] = [one(sd, kw) for sd in seeds]
        A = OUT[sc]
        tr = np.array([[o[s]["true"] for s in "AC"] for o in A])
        line = f"{sc:60s} true A/C {tr[:,0].mean():.0f}/{tr[:,1].mean():.0f} f_acc {np.mean([o['A']['true']/o['A']['n_sr'] for o in A]):.3f}"
        for nm in ("nominal", "no_nsel", "pool_J12_strat", "pool_J12_nostrat", "pool_self", "prot_alltags", "m0.01", "m0.06"):
            rr = np.array([[o[s][nm][0] / o[s]["true"] for s in "AC"] for o in A])
            line += f" | {nm} {rr.mean():.3f}±{rr.std(ddof=1)/np.sqrt(rr.size):.3f}"
        cl = np.array([[o[s]["closure"][0] / o[s]["closure"][1] for s in "AC"] for o in A])
        print(line, f"| closure {cl.mean():.3f}±{cl.std(ddof=1)/np.sqrt(cl.size):.3f} [{time.time()-t0:.0f}s]", flush=True)
        json.dump({"seeds": seeds, "scenarios": OUT}, open("toy_results.json", "w"), indent=1)
