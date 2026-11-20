#!/usr/bin/env python3
"""FITBIAS_Dplus — tabele danych i zapis uzysków: out/data.json + out/select.json + out/mc.json → tabele/data.md, out/dplus_yields.json.
Wybór (KRYTERIA §8): A2; bramka K3 (p ≥ 0,01) per komórka w kolejności A2 → A8 → A1 → A3 → A4 → A0."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fd_lib as D  # noqa: E402

O = os.path.join(D.BASE, "out")
R = json.load(open(os.path.join(O, "data.json")))
MC = {(r["bin"], r["model"]): r for r in json.load(open(os.path.join(O, "mc.json")))}
SEL = json.load(open(os.path.join(O, "select.json")))
OOS = {r["bin"]: r for r in json.load(open(os.path.join(O, "mc_oos.json"))) if r["comp"] == "pk"}
ORDER = ["A2", "A8", "A1", "A3", "A4", "A0"]
CH = "A2"
Y = {"_meta": {"zadanie": "FITBIAS_Dplus_2026-11-20", "selekcja": "D+ -> K- pi+ pi+, drzewo H, WP70, pik inkluzywny, bez AFP",
               "model_nominalny": "A2 = szablon sygnalu fd_nd z ograniczeniem Gaussa (dmu 0,5 MeV; sc 0,10) + szablon DMM+CM z mb_nd wprost z drzewa H przy WP70 (kappa_MC x 1, zwiazany z N_S) + exp(Czebyszew 3), 1,70-2,20 GeV, biny 5 MeV",
               "N": "uzysk (normalizacja szablonu sygnalu na pelnym zakresie 1,70-2,20)", "sigma_stat": "max(1,3 x HESSE, bootstrap 30 replik po zdarzeniach)",
               "sigma_model": "max(|N - N_A0|, |N - N_A4|) (+ |N - N_A2| gdy model komorki != A2)", "p": "p(chi2 Pearsona) modelu komorki",
               "fit_bias_rel": "|b| (+) sigma_b z domkniecia mb_nd (pelne tlo, WP70) modelem komorki; NIE jest w sigma_model; uzysk nie jest poprawiany o b",
               "template_stat_bias_rel": "oszacowanie obciazenia od szumu szablonu DMM+CM przy pelnej statystyce mb_nd = 0,5 x wynik poza probka z polowek (mc_oos.json); informacyjnie, miesci sie w sigma_model",
               "uwaga": "zadne ustawienie nie spelnilo K1 (5-8 GeV: b = +9 +- 6 %); wybor wg Q; post hoc"}, "cells": []}
T = ["# FITBIAS_Dplus — dane: uzyski D⁺ WP70 (κ × 1), wszystkie ustawienia", "",
     "N ± HESSE [bootstrap]; sc; dmu [MeV]; p(χ²). Wybrane: **A2**; bramka K3 p ≥ 0,01 per komórka w kolejności A2 → A8 → A1 → A3 → A4 → A0.", "",
     "| run | bin | " + " | ".join(D.SET) + " |", "|---|---|" + "---|" * len(D.SET)]
for r in sorted(R, key=lambda r: (r["run"], D.BINS.index(r["bin"]))):
    row = []
    for m in D.SET:
        o = r["models"][m]
        row.append(f"{o['N']:.0f} ± {o['err']:.0f} [{o['boot_rms']:.0f}]; {o['sc']:.2f}; {1e3 * o['dmu']:+.1f}; {o['p']:.2g}")
    T.append(f"| {r['run']} | {r['bin']} | " + " | ".join(row) + " |")
T += ["", "## Nominał → wybrany", "", "| run | bin | N A0 (obecny) ± σ | p A0 | model komórki | N ± σ_stat ± σ_model | ΔN/N | σ_stat/N | σ_model/N | p | sc ± σ | fit_bias (MC) | N kształt stały (opisowo) |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in sorted(R, key=lambda r: (r["run"], D.BINS.index(r["bin"]))):
    M = r["models"]
    cm = next(m for m in ORDER if M[m]["p"] >= 0.01)
    o, a0 = M[cm], M["A0"]
    st = max(1.3 * o["err"], o["boot_rms"])
    alts = [a0["N"], M["A4"]["N"]] + ([M[CH]["N"]] if cm != CH else [])
    sm = max(abs(o["N"] - x) for x in alts)
    f = MC[(r["bin"], cm)]["full"]
    fb = float(np.hypot(f["b"], f["b_err"]))
    st0 = max(a0["err"], a0["boot_rms"])
    T.append(f"| {r['run']} | {r['bin']} | {a0['N']:.0f} ± {st0:.0f} | {a0['p']:.3f} | {cm}{'' if cm == CH else ' (A2: p = %.2g → odrzucony)' % M[CH]['p']} | **{o['N']:.0f} ± {st:.0f} ± {sm:.0f}** | {100 * (o['N'] / a0['N'] - 1):+.1f} % | "
             f"{100 * st / o['N']:.1f} % | {100 * sm / o['N']:.1f} % | {o['p']:.3f} | {o['sc']:.2f} ± {o['sc_err']:.2f} | {100 * f['b']:+.1f} ± {100 * f['b_err']:.1f} % → {100 * fb:.1f} % | {o['fix']['N']:.0f} ({100 * (o['fix']['N'] / o['N'] - 1):+.1f} %; p = {o['fix']['p']:.2g}) |")
    Y["cells"].append({"run": r["run"], "bin": r["bin"], "model": cm, "N": round(o["N"], 1), "sigma_stat": round(st, 1), "sigma_model": round(sm, 1), "p": round(o["p"], 4),
                       "sigma_hesse": round(o["err"], 1), "sigma_boot": round(o["boot_rms"], 1), "sc": round(o["sc"], 3), "sc_err": round(o["sc_err"], 3), "dmu_MeV": round(1e3 * o["dmu"], 2),
                       "chi2": round(o["chi2"], 1), "ndf": o["ndf"], "N_A0": round(a0["N"], 1), "sigma_A0": round(st0, 1), "p_A0": round(a0["p"], 4), "dN_rel_vs_A0": round(o["N"] / a0["N"] - 1, 4),
                       "N_A4": round(M["A4"]["N"], 1), "N_A2": round(M[CH]["N"], 1), "p_A2": round(M[CH]["p"], 5), "N_fixshape": round(o["fix"]["N"], 1),
                       "template_stat_bias_rel": round(0.5 * OOS[r["bin"]]["mean"], 4), "closure_b": round(f["b"], 4), "closure_b_err": round(f["b_err"], 4), "fit_bias_rel": round(fb, 4), "K3_ok": bool(o["p"] >= 0.01)})
T += ["", "## Gdzie w masie siedzi χ² (A0 i A2; χ²/liczba binów 5 MeV, w nawiasie Σ(dane − model))", "",
      "| run | bin | ustawienie | " + " | ".join(n for n, _, _ in __import__("fd_data").REG) + " | razem χ²/ndf (p) |", "|---|---|---|" + "---|" * 7]
for r in sorted(R, key=lambda r: (r["run"], D.BINS.index(r["bin"]))):
    for m in ("A0", CH):
        o = r["models"][m]
        T.append(f"| {r['run']} | {r['bin']} | {m} | " + " | ".join(f"{g['chi2']:.1f}/{g['nb']} ({g['sum_res']:+.0f})" for g in o["reg"]) + f" | {o['chi2']:.0f}/{o['ndf']} ({o['p']:.3f}) |")
T += ["", "## Reszty znormalizowane w binach 25 MeV (A2; 1,70 → 2,20 GeV)", "", "| run | bin | " + " | ".join(f"{1.70 + 0.025 * i:.3f}" for i in range(20)) + " |", "|---|---|" + "---|" * 20]
for r in sorted(R, key=lambda r: (r["run"], D.BINS.index(r["bin"]))):
    T.append(f"| {r['run']} | {r['bin']} | " + " | ".join(f"{x:+.1f}" for x in r["models"][CH]["pull25"]) + " |")
T += ["", "## Dodatkowy Gauss w 2,0103 GeV (σ = 8 MeV; D*⁺→D⁰π⁺, D⁰→K⁻π⁺), normalizacja κ_g × N_S — skan (A2; opisowo)", "", "| run | bin | κ_g najlepsze | N_g = κ_g × N_S | −2Δln L względem κ_g = 0 | zmiana N_S |", "|---|---|---|---|---|---|"]
for r in sorted(R, key=lambda r: (r["run"], D.BINS.index(r["bin"]))):
    s = np.array(r["models"][CH]["dstar_scan"])
    i, i0 = s[:, 1].argmin(), np.abs(s[:, 0]).argmin()
    T.append(f"| {r['run']} | {r['bin']} | {s[i, 0]:+.2f} | {s[i, 0] * s[i, 2]:+.0f} | {s[i, 1] - s[i0, 1]:+.1f} | {100 * (s[i, 2] / s[i0, 2] - 1):+.1f} % |")
open(os.path.join(D.BASE, "tabele", "data.md"), "w").write("\n".join(T) + "\n")
json.dump(Y, open(os.path.join(O, "dplus_yields.json"), "w"), indent=1, ensure_ascii=False)
print("\n".join(T[16:27])); print("\n".join(T[-8:]))
