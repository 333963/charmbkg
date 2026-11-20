#!/usr/bin/env python3
"""WPSCAN — diagnostyki POST HOC (poza KRYTERIA.md): (a) stosunki przy ustalonych ε (bez wyboru optimum), (b) dryf N_prawda/ε na MC mb_nd
(czy ε z fd_nd przenosi się na mb_nd wzdłuż skanu), (c) średni ważony dryf po komórkach.  wps_extra.py wps_fit.json"""
import json, math, sys
import numpy as np
J = json.load(open(sys.argv[1])); R = J["res"]; info = J["info"]
idx = {(r["sp"], r["tag"], r["bin"], r["src"], r["eps_nom"]): r for r in R}
kind = lambda t: "drzewo" if ".xgb." in t else "sieć"
f = lambda x, n=3: "—" if x is None or not math.isfinite(x) else f"{x:+.{n}f}".replace(".", ",")
g = lambda x, n=2: "—" if x is None or not math.isfinite(x) else f"{x:.{n}f}".replace(".", ",")
FIX = [0.5, 0.35]
def rel(r):
    n = r["nom"]; return n["err"] / n["N"] if n["err"] and n["N"] > 0 and n["valid"] else float("nan")
def dd(r, r0):
    bN = np.array([x if x is not None else np.nan for x in r["boot_N"]], float); b0 = np.array([x if x is not None else np.nan for x in r0["boot_N"]], float)
    d = (r["nom"]["N"] / r["eps"]) / (r0["nom"]["N"] / r0["eps"]) - 1
    return d, float(np.nanstd((bN / r["eps"]) / (b0 / r0["eps"]) - 1))
print("### T5 (post hoc). Ustalone punkty ε = 0,50 i 0,35 względem WP70 — bez wyboru optimum\n")
print("| gatunek | próbka | pT | klasyf. | σ_HESSE/N: WP70 | stosunek (0,50) | stosunek (0,35) | d dane (0,50) | d dane (0,35) | d prawda MC mb_nd (0,50) | (0,35) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
acc = {}
for sp, tags in info["tags"].items():
    for src in ("427929|MB427929", "428770|all"):
        for b in info["bins"]:
            bn = f"{b[0]}-{b[1]}"
            for t in tags:
                r0 = idx[(sp, t, bn, src, 0.7)]; m0 = idx[(sp, t, bn, "mb|all", 0.7)]
                row = [sp, src.split("|")[0], bn, kind(t), g(rel(r0), 3)]
                rr, ds, dm = [], [], []
                for e in FIX:
                    r = idx[(sp, t, bn, src, e)]; m = idx[(sp, t, bn, "mb|all", e)]
                    rr.append(g(rel(r) / rel(r0)))
                    d, sd = dd(r, r0); ds.append(f"{f(d)} ± {g(sd, 3)}")
                    acc.setdefault((sp, kind(t), src.split("|")[0], e), []).append((d, sd))
                    dt = (m["N_truth"] / m["eps"]) / (m0["N_truth"] / m0["eps"]) - 1
                    sdt = math.sqrt(abs(1 / m["N_truth"] - 1 / m0["N_truth"])) * (1 + dt)
                    dm.append(f"{f(dt)} ± {g(sdt, 3)}")
                print("| " + " | ".join(row + rr + ds + dm) + " |")
print("\n### T6 (post hoc). Średni ważony dryf d po 3 binach pT\n")
print("| gatunek | klasyf. | próbka | d(0,50) | d(0,35) |")
print("|---|---|---|---|---|")
for sp in info["tags"]:
    for k in ("drzewo", "sieć"):
        for s in ("427929", "428770"):
            o = []
            for e in FIX:
                v = acc[(sp, k, s, e)]; w = np.array([1 / x[1] ** 2 for x in v]); d = np.array([x[0] for x in v])
                o.append(f"{f(float((w * d).sum() / w.sum()))} ± {g(float(1 / math.sqrt(w.sum())), 3)}")
            print(f"| {sp} | {k} | {s} | {o[0]} | {o[1]} |")
