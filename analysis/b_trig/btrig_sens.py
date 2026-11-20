#!/usr/bin/env python3
"""B_trig — czułość: uzysk D* per klasa triggera / run i oczekiwany błąd SF95 (skalowanie z nominału B).

Wejście: json z histogramami Δm D* (|m(Kπ) − m_D0| < 55 MeV, D0Index ≥ 0) per „<łańcuch|all|grl>|<bin>” i próbka
(`lx/dshist_root.py` na ROOT v2.0, tylko odczyt; albo dowolne źródło w tym formacie). Klasy rozłączne 428770 nie są w
histogramach (tylko łańcuchy) → przybliżenie: J12only ≈ J12 − (j20MBTS ∪ sptrk4) (górne oszacowanie nakładania: suma).
Metoda: fit_single (dstar_refit) → S (sygnał), B = tło w oknie ±2,5 MeV wokół μ; σ(SF)_klasa = σ_nom · [√(S+B)/S]_klasa /
[√(S+B)/S]_nom, gdzie nom = 428770 „all” (σ_nom: 5–8 → 0,10, 8–30 → 0,04, nominał B). Z(SF95 = 0,83 vs 1) = 0,17/σ.
Użycie: btrig_sens.py --hist dshist.json --ev-json EV_dstar/lxplus/ev_analysis_cm.json --out DIR
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import btrig_common as C  # noqa: E402

C.add_paths()
SIG_NOM = {"5-8": 0.10, "8-30": 0.04, "3-5": float("nan")}
CLS_SHOW = {"data_428770": ["all", "HLT_noalg_L1J12", "HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_mb_L1MBTS_2"],
            "data_427929": ["all", "grl", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_j20_L1MBTS_2", "HLT_noalg_L1J12"]}


def sb_from_hist(h, shape):
    from charmtag import dstar_refit as DR
    from scipy.stats import norm
    f = DR.fit_single(np.asarray(h, float), C.E_DM, shape)
    S = float(f["N"])
    e = C.E_DM
    c = 0.5 * (e[1:] + e[:-1])
    win = np.abs(c - f.get("mu", 145.45)) < 2.5
    sc = shape["sc"] * f.get("scale", 1.0)
    fs = shape["f"] * (norm.cdf(2.5 / sc) - norm.cdf(-2.5 / sc)) + (1 - shape["f"]) * (
        norm.cdf(2.5 / (sc * shape["k"])) - norm.cdf(-2.5 / (sc * shape["k"])))
    Sw = S * fs
    B = max(float(np.asarray(h)[win].sum()) - Sw, 1.0)
    return {"S": S, "S_err": float(f["N_err"]), "S_win": Sw, "B_win": B, "valid": bool(f["valid"])}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist", required=True)
    ap.add_argument("--ev-json", default=os.path.join(C.CW, "EV_dstar", "lxplus", "ev_analysis_cm.json"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    H = json.load(open(a.hist))
    shapes = {b: R["shape_mc"] for b, R in json.load(open(a.ev_json))["bins"].items()}
    res = {"samples": {}, "rows": []}
    for samp, v in H.items():
        res["samples"][samp] = {k: v.get(k) for k in ("files", "nfiles_total", "entries", "nev")}
        for key, h in v["ds"].items():
            cl, b = key.split("|")
            if cl not in CLS_SHOW.get(samp, []) or b not in shapes:
                continue
            r = sb_from_hist(h, shapes[b])
            r.update({"sample": samp, "cls": cl, "bin": b, "n_ev": v["nev"].get(cl, 0), "n_ev_all": v["entries"],
                      "S_per_kev": 1e3 * r["S"] / max(v["entries"], 1), "frac_files": v["files"] / max(v["nfiles_total"], 1)})
            res["rows"].append(r)
    ref = {r["bin"]: r for r in res["rows"] if r["sample"] == "data_428770" and r["cls"] == "all"}
    for r in res["rows"]:
        R = ref.get(r["bin"])
        sig = SIG_NOM.get(r["bin"], np.nan)
        # ekstrapolacja do pełnego runu, gdy histogram z części plików
        fS = 1.0 / r["frac_files"]
        fR = 1.0 / R["frac_files"] if R else 1.0
        if R and np.isfinite(sig) and r["S_win"] > 0 and R["S"] > 3 * R["S_err"] and r["S"] > 2 * r["S_err"]:
            r["sigma_SF95"] = C.sigma_scale(R["S_win"] * fR, R["B_win"] * fR, sig, r["S_win"] * fS, r["B_win"] * fS)
            r["Z_083_vs_1"] = 0.17 / r["sigma_SF95"]
        r["S_full"] = r["S"] * fS
        r["S_ratio_vs_428770all"] = r["S_full"] / (R["S"] * fR) if R and R["S"] > 0 else np.nan
    L = ["# B_trig — uzysk D* per trigger/run i czułość na SF95", "",
         "S = sygnał D* z fitu Δm (okno m(Kπ) ±55 MeV, bez cięcia H), ekstrapolowany do pełnego runu; σ(SF95) skalowane z nominału 428770 "
         "(5–8: 0,10; 8–30: 0,04) przez √(S+B)/S w oknie ±2,5 MeV.", "",
         "| próbka | klasa | bin | S (pełny run) | S/S(428770 all) | S/B okno | σ(SF95) | Z(0,83 vs 1) |", "|---|---|---|---|---|---|---|---|"]
    for r in sorted(res["rows"], key=lambda r: (r["bin"], r["sample"], r["cls"])):
        L.append(f"| {r['sample']} | {r['cls']} | {r['bin']} | {r['S_full']:.0f} ± {r['S_err'] / r['frac_files']:.0f} | "
                 f"{r.get('S_ratio_vs_428770all', np.nan):.2f} | {r['S_win'] / r['B_win']:.3f} | {r.get('sigma_SF95', np.nan):.2f} | "
                 f"{r.get('Z_083_vs_1', np.nan):.1f} |")
    L += ["", "Zdarzenia (pliki użyte / wszystkie):"] + [f"- {k}: {v['entries']} zd. z {v['files']}/{v['nfiles_total']} plików" for k, v in res["samples"].items()]
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(os.path.join(a.out, "btrig_sens.json"), "w"), indent=1, default=float)
    open(os.path.join(a.out, "btrig_sens.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())
