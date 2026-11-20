#!/usr/bin/env python3
"""s1b: diagnostyka niezgodności preskali efektywnej HLT_noalg_L1AFP_A_OR_C (428770): zależność od łańcucha odniesienia, LB, krotności, tagu."""
import numpy as np, os, common as C
d = C.load(428770); cl, b = C.classes(d); T0, T1, T2 = C.regions(d)
X = d["T_" + C.AN]; tbp = d["TBP_L1_AFP_A_OR_C"]; out = []
def row(lab, m):
    n = int(m.sum()); k = int((m & X).sum()); s = f"{lab:55s} n_ref {n:8d} n_pass {k:6d} p_eff " + (f"{n/k:7.2f} ± {n/k/np.sqrt(k):.2f}" if k else "—")
    print(s); out.append(s)
refs = {"J12": d["T_" + C.J12], "sptrk_pt4_L1AFP": d["T_" + C.AS], "j20_L1MBTS_2": d["T_HLT_j20_L1MBTS_2"], "sptrk_pt4_L1MBTS_2": d["T_HLT_mb_sptrk_pt4_L1MBTS_2"],
        "noalg_mb_L1MBTS_2": d["T_HLT_noalg_mb_L1MBTS_2"], "J12 & !sptrk4AFP": d["T_" + C.J12] & ~d["T_" + C.AS], "J12 & sptrk4AFP": d["T_" + C.J12] & d["T_" + C.AS],
        "sptrk4AFP & !J12 & !MBTS": d["T_" + C.AS] & ~b["J"] & ~b["M"]}
for nm, y in refs.items():
    row(f"Y={nm}, TBP_AFP", y & tbp)
    for lab, m in (("LB<=16", d["lb"] <= 16), ("LB>16", d["lb"] > 16)):
        row(f"   {lab}", y & tbp & m)
y = d["T_" + C.AS] & tbp
for lab, m in (("tag A|C", d["tag_A"] | d["tag_C"]), ("bez tagu", ~(d["tag_A"] | d["tag_C"])), ("nsel<20", d["nsel"] < 20), ("nsel 20-40", (d["nsel"] >= 20) & (d["nsel"] < 40)), ("nsel>=40", d["nsel"] >= 40),
               ("npu==0", d["npu"] == 0), ("npu>0", d["npu"] > 0), ("hasD", d["hasD"]), ("bez D", ~d["hasD"]), ("TBP_J12", d["TBP_L1_J12"]), ("!TBP_J12", ~d["TBP_L1_J12"]),
               ("TBP_MBTS_2", d["TBP_L1_MBTS_2"]), ("!TBP_MBTS_2", ~d["TBP_L1_MBTS_2"])):
    row(f"Y=sptrk_pt4_L1AFP, {lab}", y & m)
# odwrotnie: P(sptrk4AFP | noalg) i P(J12 pass | TBP_J12)
m = d["T_" + C.AN]; s = f"P(sptrk_pt4_L1AFP | noalg_L1AFP) = {(m & d['T_' + C.AS]).sum() / m.sum():.4f}; P(TBP_AFP | J12) = {(d['T_'+C.J12] & tbp).sum() / d['T_'+C.J12].sum():.4f}; P(TBP_AFP | wszystkie) = {tbp.mean():.4f}"
print(s); out.append(s)
s = f"zdarzenia z TBP_AFP bez żadnego bitu AFP: {(tbp & ~b['A']).sum()} z {tbp.sum()}; klasa none z TBP_AFP: {(cl['none'] & tbp).sum()} z {cl['none'].sum()}; none z TBP_J12: {(cl['none'] & d['TBP_L1_J12']).sum()}; none z TBP_MBTS_1: {(cl['none'] & d['TBP_L1_MBTS_1']).sum()}"
print(s); out.append(s)
open(os.path.join(C.OUT, "s1b.log"), "w").write("\n".join(out) + "\n")
