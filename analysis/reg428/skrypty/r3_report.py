#!/usr/bin/env python3
"""REG428 r3 (20.11.2026, projekt B; lokalnie) — tabele: uzyski per region, SF z tagu D*, wydajność L1_J12 dla zdarzeń z D (odniesienie MBTS 427929),
przekrój z R3 i test spójności z 427929. POST HOC. Wejście: ../out/{yields,dstar,sbeff}.json, runs/lxplus/bits_summary.json, XSEC (tylko odczyt).
Wyjście: ../tabele/tabele.md, ../out/summary.json"""
import json
import math
import os

import numpy as np
from scipy.stats import chi2 as C2

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT, TAB = os.path.join(HERE, "..", "out"), os.path.join(HERE, "..", "tabele")
BINS = ["3-5", "5-8", "8-30"]
L = {"R1": (0.4687, 0.0104), "R2": (27.65, 0.61), "R3": (29.41, 0.68), "R3hi": (26.83, 0.62), "j20AFP": (28.12, 0.62), "MB427": (0.4648, 0.0103)}   # nb^-1 (LUMI, XSEC §8.1)
BR = {"D0": (0.03951, 0.0003), "Dplus": (0.0931, 0.0014)}
N_J12_EV_427 = 61082          # zdarzenia MB ∧ GRL ∧ HLT_noalg_L1J12 w 427929 (J12EFF.md, odniesienie mbts1)

Y = {(r["sp"], r["bin"], r["run"], r["cls"]): r for r in json.load(open(os.path.join(OUT, "yields.json")))}
X = json.load(open(os.path.join(CW, "XSEC_2026-11-19", "out", "xsec_results.json")))
BS = json.load(open(os.path.join(HERE, "..", "runs", "lxplus", "bits_summary.json")))
SB = json.load(open(os.path.join(OUT, "sbeff.json")))
try:
    DSJ = json.load(open(os.path.join(OUT, "dstar.json")))
except FileNotFoundError:
    DSJ = None


def pm(x, e, nd=0):
    if x is None or e is None or not (math.isfinite(x) and math.isfinite(e)):
        return "–"
    return f"{x:.{nd}f} ± {e:.{nd}f}"


def yv(sp, bn, run, cls, model=None):
    r = Y[(sp, bn, run, cls)]
    m = model or r["chosen"]
    x = r[m]
    return x["N"], x["stat"], m, x.get("p"), x.get("at_limit")


def model_for(sp, bn):
    """Jeden model na bin dla całego runu 428770: wybór z klasy allGRL (FITBIAS: 8–30 → M0)."""
    return Y[(sp, bn, "428770", "allGRL")]["chosen"]


def chi2c(v, e):
    v, e = np.array(v, float), np.array(e, float)
    w = 1 / e ** 2
    mu = (w * v).sum() / w.sum()
    c = float((w * (v - mu) ** 2).sum())
    return mu, 1 / math.sqrt(w.sum()), c, len(v) - 1, float(C2.sf(c, len(v) - 1))


S = {}
T = ["# REG428 — tabele (generowane przez `skrypty/r3_report.py`; POST HOC; tylko bity triggera)", ""]

# ---------------------------------------------------------------- T0 zdarzenia
ce = BS["classes_events"]
T += ["## T0. Zdarzenia drzewa 428770 (wszystkie w GRL) per region bitowy", "", "| klasa | zdarzeń | udział | L [nb⁻¹] | zdarzeń / µb⁻¹ |", "|---|---|---|---|---|"]
for c, lab, lk in (("allGRL", "wszystkie", None), ("R1", "R1 `noalg_L1AFP_A_OR_C`", "R1"), ("R2", "R2 `j20_L1AFP_A_OR_C` ∧ ¬R1", "R2"), ("R3x", "R3x `noalg_L1J12` ∧ ¬R1 ∧ ¬R2", None),
                  ("rest", "reszta (żaden z trzech)", None), ("J12", "łańcuch `noalg_L1J12` (cały)", "R3"), ("j20AFP", "łańcuch `j20_L1AFP_A_OR_C` (cały)", "j20AFP"),
                  ("R1_J12", "R1 ∧ J12", None), ("R1_j20", "R1 ∧ j20_L1AFP", None), ("J12_AFPtbp", "J12 ∧ L1_AFP_A_OR_C (przed preskalą)", None), ("J12_j20AFP", "J12 ∧ j20_L1AFP", None)):
    T.append(f"| {lab} | {ce[c]} | {100 * ce[c] / ce['allGRL']:.2f} % | {L[lk][0] if lk else '–'} | {ce[c] / L[lk][0] / 1e3:.1f} |" if lk else f"| {lab} | {ce[c]} | {100 * ce[c] / ce['allGRL']:.2f} % | – | – |")
rate428, rate427 = ce["J12"] / (L["R3"][0] * 1e3), N_J12_EV_427 / (L["MB427"][0] * 1e3)
S["j12_event_rate_per_ub"] = {"428770": rate428, "427929": rate427, "ratio_428_over_427": rate428 / rate427}
T += ["", f"Zdarzenia J12 na µb⁻¹: 428770 {rate428:.1f}, 427929 (MB ∧ J12, {N_J12_EV_427} zdarzeń) {rate427:.1f} → stosunek **{rate428 / rate427:.3f}** (odpowiedź L1_J12 / pociągi paczek; XSEC §8.1: 1/0,869).",
      f"P(J12 | R1) = {100 * ce['R1_J12'] / ce['R1']:.3f} %, P(j20 | R1) = {100 * ce['R1_j20'] / ce['R1']:.2f} %, P(L1_AFP | J12) = {100 * ce['J12_AFPtbp'] / ce['J12']:.2f} %, P(j20_L1AFP | J12 ∧ L1_AFP) = {100 * ce['J12_j20AFP'] / ce['J12_AFPtbp']:.1f} %.", ""]

# ---------------------------------------------------------------- T1 uzyski
T += ["## T1. Uzyski WP70 drzewa H per region (σ = max(HESSE, rms 30 replik bootstrapu po zdarzeniach))", "",
      "D⁰: model per bin wybrany na klasie „wszystkie” (M8; 8–30 GeV M0, bo p(χ²) M8 < 0,01); w nawiasie drugi model i p(χ²) M8 w danej klasie. D⁺: M0.", "",
      "| gat. | pT | wszystkie | R1 | R2 | R3x | reszta | Σ rozłącznych / wszystkie | łańcuch J12 | łańcuch j20_L1AFP |", "|---|" + "---|" * 9]
S["yields"] = {}
for sp in ("D0", "Dplus"):
    for bn in BINS:
        m = model_for(sp, bn)
        row = [sp.replace("D0", "D⁰").replace("Dplus", "D⁺"), bn.replace("-", "–")]
        tot = 0.0
        for c in ("allGRL", "R1", "R2", "R3x", "rest"):
            n, e, _, p, lim = yv(sp, bn, "428770", c, m)
            s = pm(n, e)
            if sp == "D0":
                o = "M0" if m == "M8" else "M8"
                s += f" ({o}: {Y[(sp, bn, '428770', c)][o]['N']:.0f}; p₈ {Y[(sp, bn, '428770', c)]['M8']['p']:.2g})"
            if lim:
                s += " ᴸ"
            row.append(s)
            S["yields"][f"{sp}|{bn}|{c}"] = {"N": n, "stat": e, "model": m, "p": p}
            if c != "allGRL":
                tot += n
        row.append(f"{tot / yv(sp, bn, '428770', 'allGRL', m)[0]:.3f}")
        for c in ("J12", "j20AFP"):
            n, e, _, p, lim = yv(sp, bn, "428770", c, m)
            row.append(pm(n, e) + (" ᴸ" if lim else ""))
            S["yields"][f"{sp}|{bn}|{c}"] = {"N": n, "stat": e, "model": m, "p": p}
        T.append("| " + " | ".join(row) + " |")
T += ["", "ᴸ = skala szerokości sygnału na granicy (fit mało określony). Udziały uzysku: R3x " +
      ", ".join(f"{100 * S['yields'][f'{sp}|{bn}|R3x']['N'] / S['yields'][f'{sp}|{bn}|allGRL']['N']:.0f} %" for sp in ("D0", "Dplus") for bn in BINS) +
      "; reszta " + ", ".join(f"{100 * S['yields'][f'{sp}|{bn}|rest']['N'] / S['yields'][f'{sp}|{bn}|allGRL']['N']:.0f} %" for sp in ("D0", "Dplus") for bn in BINS) +
      " (D⁰ 3–5 / 5–8 / 8–30, D⁺ 3–5 / 5–8 / 8–30).", ""]

# ---------------------------------------------------------------- T2 SF D*
if DSJ:
    D = {(r["run"], r["cls"], r["bin"], r["wp"]): r for r in DSJ["fits"]}
    T += ["## T2. SF drzewa H z tagu D* (D⁰ z D*; `fit_passfail`, tło progowe) per region, WP70 [WP50]", "",
          "| pT | wszystkie | R3x | R2 | reszta | R1 | łańcuch J12 | ¬J12 | 427929 MB | 427929 MB ∧ J12 | zgodność R3x / R2 / reszta: χ²/ndf (p) |", "|---|" + "---|" * 10]
    S["sf"] = {}
    for bn in BINS:
        row = [bn.replace("-", "–")]
        for run, c in (("428770", "all"), ("428770", "R3x"), ("428770", "R2"), ("428770", "rest"), ("428770", "R1"), ("428770", "J12"), ("428770", "notJ12"), ("427929", "MB"), ("427929", "MB_J12")):
            r, r5 = D[(run, c, bn, 0.7)], D[(run, c, bn, 0.5)]
            ok = r["determined"] and r["valid"] and r["SF"] is not None
            s = (pm(r["SF"], r["SF_err"], 3) if ok else "nieokreślony") + f" (N {r['N']:.0f})"
            if ok and r5["determined"] and r5["valid"]:
                s += f" [{r5['SF']:.3f} ± {r5['SF_err']:.3f}]"
            row.append(s)
            S["sf"][f"{run}|{c}|{bn}"] = {"SF": r["SF"], "err": r["SF_err"], "ok": bool(ok), "N": r["N"], "n_rows": r["n_rows"]}
        v = [(D[("428770", c, bn, 0.7)]["SF"], D[("428770", c, bn, 0.7)]["SF_err"]) for c in ("R3x", "R2", "rest") if S["sf"][f"428770|{c}|{bn}"]["ok"]]
        if len(v) >= 2:
            mu, er, c2, nd, p = chi2c([a for a, _ in v], [b for _, b in v])
            row.append(f"{c2:.2f}/{nd} ({p:.2f}); średnia {mu:.3f} ± {er:.3f}")
            S["sf"][f"compat|{bn}"] = {"chi2": c2, "ndf": nd, "p": p, "mean": mu, "err": er, "n": len(v)}
        else:
            row.append("–")
        a, b = S["sf"][f"428770|J12|{bn}"], S["sf"][f"428770|notJ12|{bn}"]
        if a["ok"] and b["ok"]:
            S["sf"][f"J12_vs_not|{bn}"] = {"diff": a["SF"] - b["SF"], "err": math.hypot(a["err"], b["err"]), "z": (a["SF"] - b["SF"]) / math.hypot(a["err"], b["err"])}
        T.append("| " + " | ".join(row) + " |")
    T += ["", "J12 − ¬J12: " + "; ".join(f"{bn.replace('-', '–')}: {v['diff']:+.3f} ± {v['err']:.3f} ({v['z']:+.1f}σ)" for bn in BINS for v in [S["sf"].get(f"J12_vs_not|{bn}")] if v) + ".", ""]

# ---------------------------------------------------------------- T3 wydajność J12 dla zdarzeń z D (427929 MB)
T += ["## T3. Wydajność `HLT_noalg_L1J12` dla zdarzeń z D przy WP70 — odniesienie MBTS, run 427929 (ε = N(MB ∧ J12)/N(MB), fit piku pass i fail osobno; błąd HESSE)", "",
      "| gat. | pT | model | N(MB ∧ J12) | N(MB ∧ ¬J12) | ε_J12(D) [%] | błąd wzgl. | z odjęcia pasm [%] | kandydaci (cała masa) [%] |", "|---|" + "---|" * 8]
S["eps"] = {}
for sp in ("D0", "Dplus"):
    for bi, bn in enumerate(BINS):
        m = model_for(sp, bn)
        npv, npe, _, _, lim = yv(sp, bn, "427929", "MB_J12", m)
        nf, nfe, _, _, _ = yv(sp, bn, "427929", "MB_notJ12", m)
        e = npv / (npv + nf)
        ee = math.sqrt((nf * npe) ** 2 + (npv * nfe) ** 2) / (npv + nf) ** 2
        sb = SB["res"][f"{sp}|427929|MB_J12/MB|bins"]
        ca = SB["res"][f"{sp}|427929|MB_J12/MB|cand_all_mass"]
        S["eps"][f"{sp}|{bn}"] = {"eps": e, "err": ee, "Np": npv, "Np_err": npe, "Nf": nf, "Nf_err": nfe, "model": m, "at_limit": bool(lim), "sb": sb["eps"][bi], "sb_err": sb["eps_err"][bi],
                                  "cand": ca["pass"][bi] / ca["all"][bi]}
        T.append(f"| {sp.replace('D0', 'D⁰').replace('Dplus', 'D⁺')} | {bn.replace('-', '–')} | {m} | {pm(npv, npe, 1)}{' ᴸ' if lim else ''} | {pm(nf, nfe)} | {pm(100 * e, 100 * ee, 1)} | {100 * ee / e:.0f} % | "
                 f"{pm(100 * sb['eps'][bi], 100 * sb['eps_err'][bi], 1)} | {100 * ca['pass'][bi] / ca['all'][bi]:.1f} |")
T += ["", "| pT | ε_J12 średnia D⁰ ⊕ D⁺ [%] | błąd wzgl. | χ² (1 st. sw.) |", "|---|---|---|---|"]
for bn in BINS:
    a, b = S["eps"][f"D0|{bn}"], S["eps"][f"Dplus|{bn}"]
    mu, er, c2, nd, p = chi2c([a["eps"], b["eps"]], [a["err"], b["err"]])
    S["eps"][f"comb|{bn}"] = {"eps": mu, "err": er, "chi2": c2, "p": p}
    T.append(f"| {bn.replace('-', '–')} | {pm(100 * mu, 100 * er, 1)} | {100 * er / mu:.0f} % | {c2:.2f} |")
fe = SB["fine_edges"]
T += ["", "Zależność od pT kandydata (odjęcie pasm, 427929 MB, WP70; D⁰ przy pT < 5 GeV nieużyteczne — pasma nie opisują tła z odbiciem):", "",
      "| pT [GeV] | " + " | ".join(f"{fe[i]:g}–{fe[i + 1]:g}" for i in range(len(fe) - 1)) + " |", "|---|" + "---|" * (len(fe) - 1)]
for sp in ("D0", "Dplus"):
    r = SB["res"][f"{sp}|427929|MB_J12/MB|fine"]
    T.append(f"| {sp.replace('D0', 'D⁰').replace('Dplus', 'D⁺')} ε [%] | " + " | ".join(pm(100 * a, 100 * b, 0) for a, b in zip(r["eps"], r["eps_err"])) + " |")
T += ["", "### R1 runu 428770 jako odniesienie (metoda z zadania) — nieużyteczne", "",
      "| gat. | pT | kandydaci WP70 w R1 (cała masa) | z bitem J12 | z bitem j20_L1AFP | sygnał R1 z odjęcia pasm | sygnał ∧ J12 | sygnał ∧ j20 |", "|---|" + "---|" * 7]
for sp in ("D0", "Dplus"):
    a, b = SB["res"][f"{sp}|428770|R1_J12/R1|bins"], SB["res"][f"{sp}|428770|R1_j20/R1|bins"]
    ca, cb = SB["res"][f"{sp}|428770|R1_J12/R1|cand_all_mass"], SB["res"][f"{sp}|428770|R1_j20/R1|cand_all_mass"]
    for bi, bn in enumerate(BINS):
        T.append(f"| {sp.replace('D0', 'D⁰').replace('Dplus', 'D⁺')} | {bn.replace('-', '–')} | {ca['all'][bi]:.0f} | {ca['pass'][bi]:.0f} ({100 * ca['pass'][bi] / ca['all'][bi]:.1f} %) | {cb['pass'][bi]:.0f} ({100 * cb['pass'][bi] / cb['all'][bi]:.0f} %) | "
                 f"{pm(a['S_all'][bi], a['S_all_err'][bi])} | {pm(a['S_pass'][bi], a['S_pass_err'][bi])} | {pm(b['S_pass'][bi], b['S_pass_err'][bi])} |")
        S.setdefault("r1", {})[f"{sp}|{bn}"] = {"cand": ca["all"][bi], "cand_J12": ca["pass"][bi], "cand_j20": cb["pass"][bi], "S": a["S_all"][bi], "S_err": a["S_all_err"][bi]}
T.append("")

# ---------------------------------------------------------------- T4 przekrój z R3 i test spójności
T += ["## T4. Łańcuch J12 (R3, L = 29,41 nb⁻¹): stosunek tempa ρ = [N_J12(428770)/L] / [N(MB ∧ J12; 427929)/L] i przekrój σ_R3 = N_J12 / (C · SF · B · L · ε_J12)", "",
      "C, SF, B i σ(427929) jak w XSEC (C z mb_nd, SF z tagu D* runu 427929; D⁺ i D⁰ 3–5: SF = 1, scenariusz). ε_J12 z T3 (ten sam model fitu po obu stronach). Błąd σ_R3: stat. N ⊕ ε ⊕ L (2,3 % ⊕ 2,2 %); "
      "bez systematyki przeniesienia ε między runami (T0: tempo zdarzeń J12 różni się o 16 %) i bez systematyk wspólnych z 427929.", "",
      "| gat. | pT | N_J12 (428770) | N_J12/L [µb] | N(MB ∧ J12)/L (427929) [µb] | ρ | σ_R3 [µb] | błąd wzgl. | σ(427929) [µb], ten sam model fitu (tylko stat.; nominał XSEC D⁰ 8–30: 132 z M8) | σ_R3/σ_427929 | R3hi (bez LB 4–16): N/L [µb] |", "|---|" + "---|" * 10]
S["xsec"] = {}
for sp in ("D0", "Dplus"):
    for bn in BINS:
        m = model_for(sp, bn)
        c = X["cells"][f"{sp}|427929|{bn}"]
        n, ne, _, _, _ = yv(sp, bn, "428770", "J12", m)
        nh, nhe, _, _, _ = yv(sp, bn, "428770", "J12hi", m)
        ep = S["eps"][f"{sp}|{bn}"]
        r428, r427 = n / (L["R3"][0] * 1e3), ep["Np"] / (L["MB427"][0] * 1e3)
        rho = r428 / r427
        lrel = math.hypot(L["R3"][1] / L["R3"][0], L["MB427"][1] / L["MB427"][0])
        rho_rel = math.sqrt((ne / n) ** 2 + (ep["Np_err"] / ep["Np"]) ** 2 + lrel ** 2)
        k = c["C"] * c["SF"] * BR[sp][0]
        sig = n / (k * ep["eps"] * L["R3"][0] * 1e3)
        sig_rel = math.sqrt((ne / n) ** 2 + (ep["err"] / ep["eps"]) ** 2 + (L["R3"][1] / L["R3"][0]) ** 2)
        nmb, nmbe, _, _, _ = yv(sp, bn, "427929", "MB", m)
        s427 = nmb / (k * L["MB427"][0] * 1e3)
        S["xsec"][f"{sp}|{bn}"] = {"N_J12": n, "N_err": ne, "rho": rho, "rho_err": rho * rho_rel, "sigma_R3": sig, "sigma_R3_rel": sig_rel, "sigma_427": s427, "sigma_427_stat": s427 * nmbe / nmb,
                                   "model": m, "ratio": sig / s427, "rate_hi": nh / (L["R3hi"][0] * 1e3), "rate": r428, "C": c["C"], "SF": c["SF"]}
        T.append(f"| {sp.replace('D0', 'D⁰').replace('Dplus', 'D⁺')} | {bn.replace('-', '–')} | {pm(n, ne)} ({m}) | {pm(r428, r428 * ne / n, 3)} | {pm(r427, r427 * ep['Np_err'] / ep['Np'], 3)} | {pm(rho, rho * rho_rel, 2)} | "
                 f"{sig:.0f} ± {sig * sig_rel:.0f} | {100 * sig_rel:.0f} % | {pm(s427, s427 * nmbe / nmb)} | {pm(sig / s427, rho * rho_rel, 2)} | {pm(nh / (L['R3hi'][0] * 1e3), nhe / (L['R3hi'][0] * 1e3), 3)} |")
v = [(S["xsec"][k]["rho"], S["xsec"][k]["rho_err"]) for k in S["xsec"]]
mu, er, c2, nd, p = chi2c([a for a, _ in v], [b for _, b in v])
S["rho_all"] = {"mean": mu, "err": er, "chi2": c2, "ndf": nd, "p": p}
v2 = [(S["xsec"][k]["rho"], S["xsec"][k]["rho_err"]) for k in S["xsec"] if not k.endswith("3-5")]
mu2, er2, c22, nd2, p2 = chi2c([a for a, _ in v2], [b for _, b in v2])
S["rho_ge5"] = {"mean": mu2, "err": er2, "chi2": c22, "ndf": nd2, "p": p2}
T += ["", f"ρ średnie (6 komórek): **{mu:.2f} ± {er:.2f}**, χ²/ndf = {c2:.1f}/{nd} (p = {p:.2f}); pT ≥ 5 GeV (4 komórki): **{mu2:.2f} ± {er2:.2f}**, χ²/ndf = {c22:.1f}/{nd2} (p = {p2:.2f}). "
      f"Oczekiwane z tempa wszystkich zdarzeń J12: {S['j12_event_rate_per_ub']['ratio_428_over_427']:.2f}.",
      "σ_R3/σ_427929 ≡ ρ (te same C, SF, B; N(MB) się skraca) — to test stałości tempa D ∧ J12 na jednostkę świetlności, nie niezależny pomiar przekroju.", ""]

# ---------------------------------------------------------------- T5 R2
T += ["## T5. R2 / łańcuch `j20_L1AFP_A_OR_C` (L = 28,12 nb⁻¹) względem R1 (0,4687 nb⁻¹): ε_j20(D | L1_AFP) = [N_j20AFP/L] / [N_R1/L]", "",
      "| gat. | pT | N łańcuch j20_L1AFP | N R1 (ten sam model) | ε_j20(D \\| AFP) | błąd wzgl. |", "|---|---|---|---|---|---|"]
S["r2"] = {}
for sp in ("D0", "Dplus"):
    for bn in BINS:
        m = model_for(sp, bn)
        n, ne, _, _, _ = yv(sp, bn, "428770", "j20AFP", m)
        n1, n1e, _, _, _ = yv(sp, bn, "428770", "R1", m)
        e = (n / L["j20AFP"][0]) / (n1 / L["R1"][0]) if n1 > 0 else float("nan")
        rel = math.sqrt((ne / n) ** 2 + (n1e / n1) ** 2) if n1 > 0 else float("nan")
        S["r2"][f"{sp}|{bn}"] = {"N_j20": n, "N_j20_err": ne, "N_R1": n1, "N_R1_err": n1e, "eps": e, "rel": rel}
        T.append(f"| {sp.replace('D0', 'D⁰').replace('Dplus', 'D⁺')} | {bn.replace('-', '–')} | {pm(n, ne)} | {pm(n1, n1e)} | {e:.3f} | {100 * rel:.0f} % |")
T.append("")
os.makedirs(TAB, exist_ok=True)
open(os.path.join(TAB, "tabele.md"), "w").write("\n".join(T) + "\n")
json.dump(S, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("\n".join(T))
