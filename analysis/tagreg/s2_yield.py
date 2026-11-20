#!/usr/bin/env python3
"""s2 (część B): PRZEWIDYWANIE uzysku tagowanych D w 427929. Wejście: uzyski INKLUZYWNE z istniejących raportów (wpisane niżej ze źródłem),
ułamki tagu kandydatów z PASM BOCZNYCH (s1.json), f_acc z mieszania. Nic tu nie jest fitem regionu tagowanego."""
import json, os, numpy as np, common as C
S1 = json.load(open(os.path.join(C.OUT, "s1.json")))
# uzyski inkluzywne 427929 (GRL, wszystkie trygery); MC_templates/runs/lxplus/ana/WYNIKI.md §3 (fit C0), D*: §2a (N_S(D*) z tagu)
INC = {"WP80": {"D0": [(6773, 388), (2010, 120), (813, 51)], "Dplus": [(4999, 156), (1833, 62), (767, 36)], "Dstar": [(513, 38), (265, 18), (107, 11)]},
       "WP50": {"D0": [(2966, 118), (1160, 55), (493, 25)], "Dplus": [(2329, 72), (1123, 37), (450, 25)], "Dstar": [(296, 20), (157, 13), (74, 9)]}}
INC428 = {"WP80": {"D0": [(12570, 628), (5762, 275), (9029, 204)], "Dplus": [(10392, 269), (4890, 145), (8724, 145)], "Dstar": [(787, 60), (608, 32), (1094, 36)]}}
SB_DP70 = [0.157, 0.698, 1.336]            # S/B D+ WP70 427929 (DLZ_dplus/DATA_2026-11-17/RAPORT.md §2)
BINS = ["3-5", "5-8", "8-30"]
FACC = {"427929": (0.0087, 0.0044), "428770": (0.048, 0.014)}
# scenariusz M (model): R_okno × A·ε_AFP × ε_ξcal
R_NOM, R_LO, R_HI = 0.0075, 0.0020, 0.020
AE_NOM, AE_LO, AE_HI = 0.55, 0.35, 0.80
EPS_XC = 0.99
PM = R_NOM * AE_NOM * EPS_XC; PM_LO = R_LO * AE_LO * EPS_XC; PM_HI = R_HI * AE_HI * EPS_XC
out = {"P_model": [PM, PM_LO, PM_HI], "rows": []}
r = "427929"; t = S1[r]["tot"]; fa, dfa = FACC[r]
lines = ["| WP | gatunek | pT | N_incl (fit inkl.) | P_tag,SB [%] | S: prawdziwe, kotwica SB (×0,5–2) | M: prawdziwe, model (zakres) | przypadkowe | Z_max = √S (SB) | Z przy S/B inkl. |", "|---|---|---|---|---|---|---|---|---|---|"]
for wp in INC:
    tot = {}
    for sp in INC[wp]:
        sS = sM = sA = 0.0; z2 = 0.0
        for i, b in enumerate(BINS):
            n, dn = INC[wp][sp][i]; key = f"{sp}|{b}"
            g = t["GRL"]["sb_cands"][key]; k = t["T2"]["sb_cands"][key]; p = k / g; dp = np.sqrt(k) / g
            S = n * p * (1 - fa); dS = S * np.hypot(dn / n, dp / p); A = n * p * fa; dA = A * dfa / fa
            M = n * PM
            zmax = np.sqrt(S); zsb = ""
            if sp == "Dplus":
                sb = SB_DP70[i]; zz = S / np.sqrt(S + S / sb); zsb = f"{zz:.1f}"; z2 += zz ** 2
            lines.append(f"| {wp} | {sp} | {b} | {n} ± {dn} | {100*p:.3f} ± {100*dp:.3f} | {S:.1f} ± {dS:.1f} ({S/2:.1f}–{2*S:.1f}) | {M:.1f} ({n*PM_LO:.1f}–{n*PM_HI:.1f}) | {A:.2f} ± {dA:.2f} | {zmax:.1f} | {zsb} |")
            out["rows"].append(dict(wp=wp, sp=sp, pt=b, N=n, P_sb=p, S=S, dS=dS, M=M, M_lo=n * PM_LO, M_hi=n * PM_HI, A=A, dA=dA))
            sS += S; sM += n; sA += A
        lines.append(f"| {wp} | {sp} | **3–30** | {sM:.0f} | | **{sS:.0f}** ({sS/2:.0f}–{2*sS:.0f}) | **{sM*PM:.0f}** ({sM*PM_LO:.1f}–{sM*PM_HI:.0f}) | {sA:.2f} | {np.sqrt(sS):.1f} | {np.sqrt(z2):.1f} |" if sp == "Dplus" else
                     f"| {wp} | {sp} | **3–30** | {sM:.0f} | | **{sS:.0f}** ({sS/2:.0f}–{2*sS:.0f}) | **{sM*PM:.0f}** ({sM*PM_LO:.1f}–{sM*PM_HI:.0f}) | {sA:.2f} | {np.sqrt(sS):.1f} | |")
# porównanie z 428770: tylko liczby zdarzeń/kandydatów T2 (bez uzysków)
cmp_ = ["| wielkość (T2 = region sygnału, liczby bez fitu) | 427929 | 428770 | stosunek 428770/427929 |", "|---|---|---|---|"]
for k in ("zdarzenia", "z kand. D (slim)", "z kand. D0 pT>=3", "z D0 score>0.6 (bt all_)"):
    a = S1["427929"]["tot"]["T2"][k]; b = S1["428770"]["tot"]["T2"][k]; cmp_.append(f"| {k} | {a} | {b} | {b/a:.2f} |")
for key in S1["427929"]["tot"]["T2"]["sb_cands"]:
    a = S1["427929"]["tot"]["T2"]["sb_cands"][key]; b = S1["428770"]["tot"]["T2"]["sb_cands"][key]
    if a: cmp_.append(f"| kandydaci SB {key} | {a} | {b} | {b/a:.1f} |")
# oczekiwanie analogiczne (kotwica SB) dla 428770 — wyłącznie do porównania rzędu wielkości
t8 = S1["428770"]["tot"]; fa8 = FACC["428770"][0]; l8 = ["| gatunek | pT | N_incl 428770 WP80 | P_tag,SB [%] | S (kotwica SB) | przypadkowe | S(428770)/S(427929) |", "|---|---|---|---|---|---|---|"]
for sp in INC428["WP80"]:
    for i, b in enumerate(BINS):
        n, dn = INC428["WP80"][sp][i]; key = f"{sp}|{b}"; p = t8["T2"]["sb_cands"][key] / t8["GRL"]["sb_cands"][key]
        S8 = n * p * (1 - fa8); s9 = [x["S"] for x in out["rows"] if x["wp"] == "WP80" and x["sp"] == sp and x["pt"] == b][0]
        l8.append(f"| {sp} | {b} | {n} ± {dn} | {100*p:.3f} | {S8:.0f} | {n*p*fa8:.1f} | {S8/s9:.1f} |")
open(os.path.join(C.OUT, "s2_tables.md"), "w").write("\n".join(lines) + "\n\n" + "\n".join(cmp_) + "\n\n" + "\n".join(l8) + "\n")
C.jdump(out, os.path.join(C.OUT, "s2.json")); print(open(os.path.join(C.OUT, "s2_tables.md")).read())
