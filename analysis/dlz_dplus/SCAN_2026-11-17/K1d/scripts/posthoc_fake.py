"""K1d POST HOC: fałszywy sygnał na samym tle MC (fit z kształtem ustalonym, pol3) vs WP, w % B w oknie i w σ bootstrapu; Δ = b_DLZ − b_H vs WP. → posthoc.md"""
import json, os
import numpy as np
K1D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = json.load(open(os.path.join(K1D, "out", "k1d.json")))["cells"]
L = ["# K1d POST HOC — fałszywy sygnał na samym tle MC i Δ(DLZ − H) w funkcji WP (pol3)\n",
     "N_fake = N_fit na histogramie bez prawdziwego sygnału (kształt ustalony); σ_boot = odch. std. (N_fit − N_true) w replikach; Δ = b_DLZ − b_H ze sparowanym błędem.\n",
     "| wariant | pT | klas. | WP50: N_fake (% B_okna; /σ_boot) | WP70 | WP80 | Δ WP50 | Δ WP70 | Δ WP80 |", "|---|---|---|---|---|---|---|---|---|"]
for var in ["a", "b1_427929", "b1_428770"]:
    for bn in ["3-5", "5-8", "8-30"]:
        D = []
        for wp in [0.5, 0.7, 0.8]:
            x = {}
            for c in ["H", "DLZ"]:
                q = C[f"{var}|{bn}|{c}|{wp}"]["fit"]["3"]
                x[c] = (q["b"], np.array(q["boot"]["N"]) / np.array(q["boot"]["nt"]), np.array(q["boot"]["ok"]) > 0)
            ok = x["H"][2] & x["DLZ"][2]
            d, sd = x["DLZ"][0] - x["H"][0], (x["DLZ"][1][ok] - x["H"][1][ok]).std(ddof=1)
            D.append(f"{d:+.3f} ± {sd:.3f} ({d/sd:+.1f}σ)")
        for c in ["H", "DLZ"]:
            s = []
            for wp in [0.5, 0.7, 0.8]:
                cell = C[f"{var}|{bn}|{c}|{wp}"]
                q = cell["fit"]["3"]
                l = q["loo"]["bkg_only_fixshape"]
                N, nt, ok = np.array(q["boot"]["N"]), np.array(q["boot"]["nt"]), np.array(q["boot"]["ok"]) > 0
                sb = (N[ok] - nt[ok]).std(ddof=1)
                s.append(f"{l['dN']:+.0f} ({100*l['dN']/cell['B_win']:+.1f} %; {l['dN']/sb:+.1f})")
            L.append(f"| {var} | {bn} | {c} | " + " | ".join(s) + " | " + (" | ".join(D) if c == "DLZ" else " | | ") + " |")
open(os.path.join(K1D, "posthoc.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
