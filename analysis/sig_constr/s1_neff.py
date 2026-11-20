#!/usr/bin/env python3
"""SIG_constr A1: ile efektywnych zdarzeń sygnału daje tag D* — N_sig i n_eff (SBS, sWeights) per klasa × obszar (PRE / H06) × bin pT."""
from sc_common import *

M = load("fd_nd"); shapes = mc_shapes(M)
rows = []
for name, run in (("data_427929", 427929), ("data_428770", 428770)):
    D = load(name); cls = data_classes(D, run)
    for cn, cm in cls.items():
        for reg, sthr in (("PRE", None), ("H06", 0.6)):
            for b, (lo, hi) in list(PTB.items()) + [("3-30", (3.0, 30.0))]:
                sel = cm & base_sel(D, sthr, lo, hi)
                n = int(sel.sum())
                if n < 30:
                    continue
                wb = sbs_weights(D, sel)
                dm = D["dm"][sel] * 1e3
                nwin = int(((dm > IC.SIG[0]) & (dm < IC.SIG[1])).sum())
                if b != "3-30" and wb.sum() / max(nwin, 1) >= 0.25 and wb.sum() > 80:     # sWeights tylko przy czystości okna ≥ 0,25 (niżej fit niestabilny)
                    ws, inf1 = sweights(D["dm"][sel], shapes[b]); nserr = inf1["Ns_err"]
                else:
                    ws, nserr = np.zeros(1), float("nan")
                rows.append({"cls": cn, "reg": reg, "bin": b, "n_pairs": n, "n_win": nwin, "Ns_sbs": float(wb.sum()), "Ns_sbs_err": float(np.sqrt((wb ** 2).sum())),
                             "neff_sbs": neff(wb), "pur_win": float(wb.sum() / max(nwin, 1)), "Ns_sw": float(ws.sum()), "Ns_sw_err": nserr, "neff_sw": neff(ws)})
                r = rows[-1]
                print(f"{cn:16s} {reg} {b:5s} par={n:7d} okno={nwin:6d} SBS {r['Ns_sbs']:7.0f}±{r['Ns_sbs_err']:4.0f} n_eff {r['neff_sbs']:6.0f} czyst {r['pur_win']:.3f} | sW {r['Ns_sw']:7.0f}±{r['Ns_sw_err']:4.0f} n_eff {r['neff_sw']:6.0f}", flush=True)
# MC odniesienie: liczba prawdziwych par sygnałowych
mc = {}
for reg, sthr in (("PRE", None), ("H06", 0.6)):
    for b, (lo, hi) in list(PTB.items()) + [("3-30", (3.0, 30.0))]:
        mc[f"{reg}|{b}"] = int(mc_sig(M, sthr, lo, hi).sum())
print(mc)
jdump({"rows": rows, "mc_sig": mc, "shapes": shapes}, "s1_neff.json")
