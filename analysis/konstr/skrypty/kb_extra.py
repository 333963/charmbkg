#!/usr/bin/env python3
"""KONSTR B — dodatek po fitach okna (POST HOC, poza kryteriami): ułamek f = N_tag/N_incl z mianownikiem tego samego modelu oraz wariant „wszystkie szerokie składowe w tle”
(model nietagowany M7 z f_refl = 0: sygnał = sam wąski pik; odpowiednik rachunku okno − pasma). → tabele/B_extra.md"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kb_tag as B
T = B.T
R = B.load_incl_rows(); TG = B.tagged(blind=False)
sp = "D0"; r = R["sp"] == 0; t = TG[sp]; Rm, Rn, Ri = R["m"][r], R["nsel"][r], R["ipt"][r]
w = np.ones(len(Rm))
L = ["# B — dodatek (post hoc, poza kryteriami): D⁰ 427929 MBTS, tło z konstrukcji, wagi: brak", "", "| model rozkładu nietagowanego (definicja N) | zakres | N_incl (ten sam model) | N_tag | −/+ | Z | f_raw [%] |", "|---|---|---|---|---|---|---|"]
orig = B.V.wcat


def nore(*a, **k):
    c = orig(*a, **k); c["frefl"] = 0.0; c["refT"] = None
    return c
for lab, mname, rg, patch in (("M8 (sygnał + odbicie + DMM/CM związane)", "M8", None, False), ("M7 (DMM/CM w tle)", "M7", None, False), ("M0 (DMM/CM w tle)", "M0", (1.75, 2.00), False),
                              ("M7 bez odbicia (wszystkie szerokie składowe w tle; N = wąski pik)", "M7", None, True), ("M0 bez odbicia", "M0", (1.75, 2.00), True)):
    B.V.wcat = nore if patch else orig
    Bh, Sm, Nw = B.incl_model(Rm, Ri, w, sp, mname, rg)
    q = B.CFit(Sm, Bh).full(T.hist5(t["m"], sp, rg), cl95=False)
    L.append(f"| {lab} | {'1,70–2,05' if rg is None else '1,75–2,00'} | {Nw:.0f} | {q['N']:.1f} | −{q['N'] - q['lo68']:.1f} +{q['hi68'] - q['N']:.1f} | {q['Z']:.2f} | {100 * q['N'] / Nw:.2f} −{100 * (q['N'] - q['lo68']) / Nw:.2f} +{100 * (q['hi68'] - q['N']) / Nw:.2f} |")
open(os.path.join(B.T_, "tabele", "B_extra.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
