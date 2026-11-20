#!/usr/bin/env python3
"""FITBIAS — skutek dla danych (KRYTERIA §2.4). Wejście: out/var_data.json, out/var_mc.json, out/select.json → tabele/data.md, out/data_summary.json"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
O = lambda n: os.path.join(HERE, "..", "out", n)  # noqa: E731
f = lambda x, n=3: "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:+.{n}f}".replace(".", ",")  # noqa: E731
u = lambda x, n=3: "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{n}f}".replace(".", ",")  # noqa: E731
MS = [f"M{i}" for i in range(10)]
BINS = ["3-5", "5-8", "8-30"]


def wmean(v, e):
    v, e = np.array(v, float), np.array(e, float)
    w = 1 / e ** 2
    return float((v * w).sum() / w.sum()), float(1 / math.sqrt(w.sum()))


def rms(x):
    x = [v for v in x if v is not None]
    return float(np.std(x)) if len(x) > 5 else float("nan")


D = json.load(open(O("var_data.json")))["res"]
MC = {(r["sp"], r["bin"]): r for r in json.load(open(O("var_mc.json")))["res"]}
SEL = json.load(open(O("select.json")))
from scipy.stats import chi2 as C2  # noqa: E402
L = ["# FITBIAS — skutek dla danych: stary model (M0) vs wybrany; wszystko post hoc", ""]
summ = {}
# ---- R1
L += ["## R1. Uzyski przy WP70: M0 vs wybrany (σ = max(HESSE, rms 30 replik bootstrapu po zdarzeniach))", "",
      "| gat. | run | pT | wybrany | N (M0) ± HESSE [boot] | sc; dmu MeV; p(χ²) | N (nowy) ± HESSE [boot] | sc; dmu MeV; p(χ²) | ΔN/N | σ_N/N: M0 → nowy | σ_nowy/σ_M0 |", "|---|" + "---|" * 10]
for r in D:
    ch = SEL[r["sp"]]["choice"]
    o = r["S"]["WP70"]
    a, b = o["M0"], o[ch]
    ba, bb = rms(a.get("boot", [])), rms(b.get("boot", []))
    sa, sb = max(a["err"] or 0, ba if math.isfinite(ba) else 0), max(b["err"] or 0, bb if math.isfinite(bb) else 0)
    g = lambda x: f"{u(x['sc'], 2)}; {f(1000 * x['dmu'], 1)}; {C2.sf(x['chi2'], x['ndf']):.2g}"  # noqa: E731
    summ[f"{r['sp']}|{r['run']}|{r['bin']}"] = {"N0": a["N"], "s0": sa, "N1": b["N"], "s1": sb, "choice": ch}
    L.append(f"| {r['sp']} | {r['run']} | {r['bin']} | {ch} | {a['N']:.0f} ± {a['err']:.0f} [{ba:.0f}] | {g(a)} | {b['N']:.0f} ± {b['err']:.0f} [{bb:.0f}] | {g(b)} | {f(b['N'] / a['N'] - 1)} | {u(sa / a['N'])} → {u(sb / b['N'])} | {u((sb / b['N']) / (sa / a['N']), 2)} |")
L += ["", "## R1b. Uzyski WP70 wszystkimi ustawieniami (opisowo; N ± HESSE)", "", "| gat. | run | pT | " + " | ".join(MS) + " |", "|---|---|---|" + "---|" * 10]
for r in D:
    o = r["S"]["WP70"]
    L.append(f"| {r['sp']} | {r['run']} | {r['bin']} | " + " | ".join(f"{o[m]['N']:.0f} ± {(o[m]['err'] or float('nan')):.0f}" for m in MS) + " |")
# ---- R2 SF_rel
L += ["", "## R2. Pass/fail: SF_rel = f_dane/f_MC(fd) (σ = max(HESSE, boot)); to samo na mb_nd; stosunek podwójny dane/MC-fit", ""]
dd = {}
for mode in ("M0", "new"):
    L += [f"### {'stary model M0' if mode == 'M0' else 'wybrany model'}", "", "| gat. | run | ε_t | 3–5 | 5–8 | 8–30 | średnia dane | średnia MC-fit (prawda MC) | podwójny stosunek dane/MC |", "|---|" + "---|" * 8]
    for sp in ("D0", "Dplus"):
        m = "M0" if mode == "M0" else SEL[sp]["choice"]
        for run in ("427929", "428770"):
            for et in ("0.5", "0.4", "0.35"):
                v, e, mv, me, tv = [], [], [], [], []
                for bn in BINS:
                    r = next(x for x in D if x["sp"] == sp and x["run"] == run and x["bin"] == bn)
                    o = r["J"][et]; x = o[m]
                    b = rms(x.get("boot_f", []))
                    s = max(x.get("f_err") or 0, b if math.isfinite(b) else 0)
                    v.append(x["f"] / o["f_mc"]); e.append(s / o["f_mc"])
                    q = MC[(sp, bn)]["J"][et]; y = q[m]
                    mv.append(y["f"] / q["f_mc"]); me.append(max(y.get("f_err") or 0, y.get("f_rms") or 0) / q["f_mc"]); tv.append(q["f_true"] / q["f_mc"])
                a, ae = wmean(v, e); c, ce = wmean(mv, me)
                dr, dre = a / c, a / c * math.hypot(ae / a, ce / c)
                dd[(mode, sp, run, et)] = {"data": (a, ae), "mc": (c, ce), "truth": float(np.mean(tv)), "dr": (dr, dre), "cells": list(zip(v, e)), "mc_cells": list(zip(mv, me))}
                L.append(f"| {sp} | {run} | {et.replace('.', ',')} | " + " | ".join(f"{u(a_)} ± {u(b_)}" for a_, b_ in zip(v, e)) + f" | {u(a)} ± {u(ae)} | {u(c)} ± {u(ce)} ({u(float(np.mean(tv)))}) | {u(dr, 2)} ± {u(dre, 2)} |")
    L.append("")
# ---- R3 drift
L += ["## R3. Dryf d(ε) = [N/ε_MC](ε)/[N/ε_MC](0,70) − 1 z fitów jednej kategorii: dane (błąd: rms sparowanych replik) i mb_nd-fit (błąd przybliżony √|σ_a² − σ_b²|)", "",
      "| gat. | próbka | model | d(0,50): 3–5 | 5–8 | 8–30 | średnia | d(0,35): 3–5 | 5–8 | 8–30 | średnia |", "|---|" + "---|" * 10]
drift = {}
for sp in ("D0", "Dplus"):
    for mode in ("M0", "new"):
        m = "M0" if mode == "M0" else SEL[sp]["choice"]
        for src in ("427929", "428770", "mb_nd"):
            cells = {}
            for wp in ("WP50", "WP35"):
                v, e = [], []
                for bn in BINS:
                    if src == "mb_nd":
                        S = MC[(sp, bn)]["S"]
                        a, b = S[wp][m], S["WP70"][m]
                        d = (a["N"] / S[wp]["eps"]) / (b["N"] / S["WP70"]["eps"]) - 1
                        ra, rb = (a["err"] or float("nan")) / a["N"], (b["err"] or float("nan")) / b["N"]
                        v.append(d); e.append((1 + d) * math.sqrt(abs(ra ** 2 - rb ** 2)))
                    else:
                        r = next(x for x in D if x["sp"] == sp and x["run"] == src and x["bin"] == bn)
                        a, b = r["S"][wp][m], r["S"]["WP70"][m]
                        ea, eb = r["S"][wp]["eps"], r["S"]["WP70"]["eps"]
                        d = (a["N"] / ea) / (b["N"] / eb) - 1
                        reps = [(x / ea) / (y / eb) - 1 for x, y in zip(a.get("boot", []), b.get("boot", [])) if x is not None and y is not None and y > 0]
                        v.append(d); e.append(float(np.std(reps)) if len(reps) > 5 else float("nan"))
                cells[wp] = (v, e, wmean(v, e) if all(np.isfinite(e)) else (float("nan"), float("nan")))
            drift[(sp, mode, src)] = {k: v[2] for k, v in cells.items()}
            L.append(f"| {sp} | {src} | {m} | " + " | ".join(" | ".join(f"{f(x)} ± {u(y)}" for x, y in zip(*cells[wp][:2])) + f" | {f(cells[wp][2][0])} ± {u(cells[wp][2][1])}" for wp in ("WP50", "WP35")) + " |")
L.append("")
# ---- R4 skan pełny (dane), M0 vs nowy: N/ε względem WP70
L += ["## R4. Skan: [N/ε_MC](ε)/[N/ε_MC](0,70) − 1 w wybranych punktach (bez błędów; HESSE N w R1)", "", "| gat. | próbka | pT | model | 0,90 | 0,80 | 0,60 | 0,50 | 0,40 | 0,35 | 0,30 | 0,20 |", "|---|" + "---|" * 11]
for sp in ("D0", "Dplus"):
    for src in ("427929", "428770", "mb_nd"):
        for bn in BINS:
            for mode in ("M0", "new"):
                m = "M0" if mode == "M0" else SEL[sp]["choice"]
                S = MC[(sp, bn)]["S"] if src == "mb_nd" else next(x for x in D if x["sp"] == sp and x["run"] == src and x["bin"] == bn)["S"]
                ref = S["WP70"][m]["N"] / S["WP70"]["eps"]
                L.append(f"| {sp} | {src} | {bn} | {m} | " + " | ".join(f((S[w][m]["N"] / S[w]["eps"]) / ref - 1, 2) for w in ("WP90", "WP80", "WP60", "WP50", "WP40", "WP35", "WP30", "WP20")) + " |")
open(os.path.join(HERE, "..", "tabele", "data.md"), "w").write("\n".join(L) + "\n")
json.dump({"wp70": summ, "sf": {"|".join(k): v for k, v in dd.items()}, "drift": {"|".join(k): v for k, v in drift.items()}}, open(O("data_summary.json"), "w"))
print("OK tabele/data.md")
