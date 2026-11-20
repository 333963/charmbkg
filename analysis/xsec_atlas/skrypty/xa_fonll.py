#!/usr/bin/env python3
"""XSEC_atlas A — calki FONLL (formularz sieciowy, pliki fonll/*.txt; metoda calkowania jak XSEC_2026-11-19/xsec2/fonll/integ.py) i skalowanie energii
7 TeV (|eta| < 2,1) -> 13,6 TeV (|eta| < 2,2) dla D*: K = [f sigma_c + sigma_B->D*](13,6) / [...](7); f(c -> D*+) = 0,2429 (1509.01061), wariant 0,156 (ALICE)."""
import glob, json, math, os
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def load(fn):
    return [[float(x) for x in l.split()] for l in open(fn) if l.strip() and not l.startswith("#") and l.split()[0][0].isdigit()]
def integ(rows, col):
    s = 0
    for a, b in zip(rows[:-1], rows[1:]):
        x1, x2, y1, y2 = a[0], b[0], a[col], b[col]
        k = math.log(y2 / y1) / math.log(x2 / x1)
        s += y1 * x1 * math.log(x2 / x1) if abs(k + 1) < 1e-9 else y1 / (k + 1) * (x2 * (x2 / x1) ** k - x1)
    return s / 1e6                                   # pb -> ub
I = {}
for f in sorted(glob.glob(os.path.join(H, "fonll", "*_pt*.txt"))):
    r = load(f)
    if len(r) < 3:
        print("PUSTY", os.path.basename(f)); continue
    I[os.path.basename(f)[:-4]] = [integ(r, 1), integ(r, 2), integ(r, 3)]
out = {"integrals_ub": I, "K": {}}
L = ["| bin pT | FONLL 7 TeV, \\|η\\| < 2,1: c→D* (f = 1) / B→D* [µb] | FONLL 13,6 TeV, \\|η\\| < 2,2 [µb] | σ_FONLL(D*±) 13,6 TeV [µb] (min–max skali c) | K (f = 0,2429) | K (f = 0,156) | K: min/min, max/max |", "|---|---|---|---|---|---|---|"]
for bn in ("3.5-5", "5-8", "8-30"):
    try:
        c7, b7 = I[f"prompt_Dst_c7_eta2.1_pt{bn}"], I[f"BtoDst_c7_eta2.1_pt{bn}"]; c19, b19 = I[f"prompt_Dst_c19_eta2.2_pt{bn}"], I[f"BtoDst_c19_eta2.2_pt{bn}"]
    except KeyError as e:
        print("brak", e); continue
    K = lambda f, i=0: (f * c19[i] + b19[i]) / (f * c7[i] + b7[i])  # noqa: E731
    f = 0.2429
    out["K"][bn] = {"K": K(f), "K_alice": K(0.156), "K_min": K(f, 1), "K_max": K(f, 2), "sigma13p6": [f * c19[i] + b19[i] for i in range(3)], "sigma7": [f * c7[i] + b7[i] for i in range(3)],
                    "fb13p6": b19[0] / (f * c19[0] + b19[0])}
    L.append(f"| {bn} | {c7[0]:.1f} / {b7[0]:.2f} | {c19[0]:.1f} / {b19[0]:.2f} | {f * c19[0] + b19[0]:.1f} ({f * c19[1] + b19[1]:.1f}–{f * c19[2] + b19[2]:.1f}) | {K(f):.3f} | {K(0.156):.3f} | {K(f, 1):.3f} / {K(f, 2):.3f} |")
print("\n".join(L))
json.dump(out, open(os.path.join(H, "out/fonll.json"), "w"), indent=1)
open(os.path.join(H, "tabele/fonll_skalowanie.md"), "w").write("# FONLL: skalowanie energii 7 → 13,6 TeV dla D*± (formularz sieciowy FONLL v1.3.4, NNPDF30 NLO, pasmo skal ⊕ mas ⊕ PDF)\n\n" + "\n".join(L) + "\n")
