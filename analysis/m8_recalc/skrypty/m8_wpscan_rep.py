#!/usr/bin/env python3
"""Tabele pkt 1 z out/wpscan_m8.json (+ domknięcie mb_nd z FITBIAS/out/var_mc.json) → tabele/wpscan_m8.md, out/wpscan_m8_summary.json"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "..", "out", "wpscan_m8.json")))["res"]
MC = {r["bin"]: r for r in json.load(open(os.path.join(HERE, "..", "..", "FITBIAS_2026-11-19", "out", "var_mc.json")))["res"] if r["sp"] == "D0"}
PMIN = 1e-3
pl = lambda x, n=3: f"{x:.{n}f}".replace(".", ",")  # noqa: E731


def arr(x):
    return np.array([np.nan if v is None else v for v in x], float)


def rel(o):
    """σ/N: HESSE, boot rms, max."""
    b = arr(o["boot_N"])
    rms = float(np.nanstd(b))
    return o["err"] / o["N"], rms / o["N"], max(o["err"], rms) / o["N"]


def ratio(c, m, e):
    """stosunek [σ_H/N](ε)/[σ_H/N](0,70): wartość, rms po replikach; to samo dla max(H, boot)."""
    a, b = c["S"][str(e)][m], c["S"]["0.7"][m]
    r0 = (a["err"] / a["N"]) / (b["err"] / b["N"])
    rr = (arr(a["boot_err"]) / arr(a["boot_N"])) / (arr(b["boot_err"]) / arr(b["boot_N"]))
    return r0, float(np.nanstd(rr)), rel(a)[2] / rel(b)[2], rr


def drift(c, m, e):
    a, b = c["S"][str(e)][m], c["S"]["0.7"][m]
    ea, eb = c["S"][str(e)]["eps"], c["S"]["0.7"]["eps"]
    d0 = (a["N"] / ea) / (b["N"] / eb) - 1
    dd = (arr(a["boot_N"]) / ea) / (arr(b["boot_N"]) / eb) - 1
    return d0, float(np.nanstd(dd)), dd


L = ["# M8_przeliczenie — pkt 1: skan punktu pracy drzewa H dla D⁰, M0 (stary) vs M8 (κ × 1,364); wszystko post hoc", ""]
S = {}
# kontrola adekwatności M8 wzdłuż skanu
L += ["## W0. Adekwatność M8 wzdłuż skanu (dane): punkty z p(χ²) < 10⁻³, na granicy dmu/sc albo nieważne (ε ∈ [0,10; 0,95], 18 punktów)", "",
      "| run | pT | p(χ²) M8: 0,70 / 0,50 / 0,35 | p(χ²) M0: 0,70 / 0,50 / 0,35 | punkty M8 odrzucone (ε: p) | granica / nieważny M8 | werdykt |", "|---|---|---|---|---|---|---|"]
REJ = {}
for c in R:
    bad = [(e, o["M8"]["p"]) for e, o in c["S"].items() if o["M8"]["p"] < PMIN]
    lim = [e for e, o in c["S"].items() if o["M8"]["at_limit"] or not o["M8"]["valid"]]
    rej = any(float(e) in (0.7, 0.5, 0.35) for e, _ in bad)
    REJ[(c["run"], c["bin"])] = rej
    L.append(f"| {c['run']} | {c['bin']} | " + " / ".join(f"{c['S'][e]['M8']['p']:.2g}" for e in ("0.7", "0.5", "0.35")) + " | " + " / ".join(f"{c['S'][e]['M0']['p']:.2g}" for e in ("0.7", "0.5", "0.35")) +
             f" | {'; '.join(f'{e}: {p:.1g}' for e, p in bad) or '—'} | {', '.join(lim) or '—'} | {'**dane odrzucają szablon → zostaje M0**' if rej else 'M8'} |")
L += ["", "## W1. Błąd statystyczny uzysku: σ_N/N przy WP70 i stosunek do WP70 przy ε = 0,50 / 0,35 (H = sam HESSE, jak T5 WPSCAN; max = max(HESSE, bootstrap 30 replik))", "",
      "| run | pT | model | σ/N WP70: H / max | stosunek H (0,50) | stosunek H (0,35) | stosunek max (0,50) | (0,35) | σ/N max przy 0,50 / 0,35 | Δ stosunku H (M8 − M0): 0,50; 0,35 [σ sparowane] |", "|---|---|---|---|---|---|---|---|---|---|"]
for c in R:
    key = f"{c['run']}|{c['bin']}"
    S[key] = {"m8_rejected": REJ[(c["run"], c["bin"])]}
    rr = {}
    for m in ("M0", "M8"):
        h70, _, x70 = rel(c["S"]["0.7"][m])
        q = {e: ratio(c, m, e) for e in (0.5, 0.35)}
        rr[m] = q
        S[key][m] = {"relH70": h70, "relmax70": x70, **{f"rH_{e}": q[e][0] for e in q}, **{f"rH_err_{e}": q[e][1] for e in q}, **{f"rmax_{e}": q[e][2] for e in q},
                     **{f"relmax_{e}": rel(c["S"][str(e)][m])[2] for e in q}, **{f"relH_{e}": rel(c["S"][str(e)][m])[0] for e in q}}
    dl = []
    for e in (0.5, 0.35):
        d = rr["M8"][e][0] - rr["M0"][e][0]
        s = float(np.nanstd(rr["M8"][e][3] - rr["M0"][e][3]))
        dl.append(f"{d:+.2f} ({d / s:+.1f}σ)")
        S[key][f"dr_{e}"] = [d, s]
    for m in ("M0", "M8"):
        v, q = S[key][m], rr[m]
        tag = " (odrzucony)" if m == "M8" and REJ[(c["run"], c["bin"])] else ""
        L.append(f"| {c['run']} | {c['bin']} | {m}{tag} | {pl(v['relH70'])} / {pl(v['relmax70'])} | {pl(q[0.5][0], 2)} ± {pl(q[0.5][1], 2)} | {pl(q[0.35][0], 2)} ± {pl(q[0.35][1], 2)} | {pl(q[0.5][2], 2)} | {pl(q[0.35][2], 2)} | "
                 f"{pl(v['relmax_0.5'])} / {pl(v['relmax_0.35'])} | {'; '.join(dl) if m == 'M8' else ''} |")
L += ["", "Błąd stosunku H = rms po 30 replikach (repliki wspólne dla wszystkich ε i obu modeli). Stosunek „max” nie ma błędu (rms z 30 replik ma sam błąd ±13 %).", "",
      "## W2. N, N/ε_MC i dryf d(ε) = [N/ε](ε)/[N/ε](0,70) − 1 (błąd: rms sparowanych replik)", "",
      "| run | pT | model | N(0,70) | N(0,50) | N(0,35) | d(0,50) | d(0,35) |", "|---|---|---|---|---|---|---|---|"]
for c in R:
    key = f"{c['run']}|{c['bin']}"
    for m in ("M0", "M8"):
        d = {e: drift(c, m, e) for e in (0.5, 0.35)}
        S[key][m].update({f"d_{e}": d[e][0] for e in d}); S[key][m].update({f"d_err_{e}": d[e][1] for e in d})
        S[key][m].update({f"N_{e}": c["S"][str(e)][m]["N"] for e in (0.7, 0.5, 0.35)})
        S[key][m].update({f"Nerr_{e}": rel(c["S"][str(e)][m])[2] * c["S"][str(e)][m]["N"] for e in (0.7, 0.5, 0.35)})
        L.append(f"| {c['run']} | {c['bin']} | {m} | " + " | ".join(f"{S[key][m][f'N_{e}']:.0f} ± {S[key][m][f'Nerr_{e}']:.0f}" for e in (0.7, 0.5, 0.35)) +
                 f" | {d[0.5][0]:+.3f} ± {d[0.5][1]:.3f} | {d[0.35][0]:+.3f} ± {d[0.35][1]:.3f} |")
# średnie ważone dryfu po binach
L += ["", "Średni ważony dryf po 3 binach (nominał: M8, w 428770 8–30 GeV M0):", "", "| run | wariant | d(0,50) | d(0,35) |", "|---|---|---|---|"]
for run in ("427929", "428770"):
    for lab, pick in (("M0", lambda k: "M0"), ("M8 wszędzie", lambda k: "M8"), ("nominał (M8; M0 gdzie odrzucony)", lambda k: "M0" if S[k]["m8_rejected"] else "M8")):
        row = []
        for e in (0.5, 0.35):
            ks = [k for k in S if k.startswith(run)]
            v = np.array([S[k][pick(k)][f"d_{e}"] for k in ks]); s = np.array([S[k][pick(k)][f"d_err_{e}"] for k in ks])
            w = 1 / s ** 2
            row.append(f"{(v * w).sum() / w.sum():+.3f} ± {1 / math.sqrt(w.sum()):.3f}")
            S.setdefault(f"mean|{run}|{lab}", {})[f"d_{e}"] = [float((v * w).sum() / w.sum()), float(1 / math.sqrt(w.sum()))]
        L.append(f"| {run} | {lab} | {row[0]} | {row[1]} |")
L += ["", "## W3. Pass/fail wewnątrz WP70 (fit łączny A = ε ≥ ε_t, B = WP70 ∧ ¬A): stosunek błędów [E = N_A + N_B]/[N_70], HESSE (jak WPTIGHT) i max(HESSE, boot)", "",
      "| run | pT | model | H: ε_t = 0,50 / 0,40 / 0,35 | max: 0,50 / 0,40 / 0,35 | E/N_70 − 1 (0,50) | p(χ²) (0,50) |", "|---|---|---|---|---|---|---|"]
for c in R:
    key = f"{c['run']}|{c['bin']}"
    for m in ("M0", "M8"):
        s70 = c["S"]["0.7"][m]
        h, x = [], []
        for et in ("0.5", "0.4", "0.35"):
            j = c["J"][et][m]
            h.append(j["E_err"] / s70["err"] if j.get("E_err") else float("nan"))
            x.append(max(j.get("E_err") or 0, float(np.nanstd(arr(j["boot_E"])))) / (rel(s70)[2] * s70["N"]))
        S[key][m]["EH"] = h; S[key][m]["Emax"] = x
        j = c["J"]["0.5"][m]
        L.append(f"| {c['run']} | {c['bin']} | {m} | " + " / ".join(pl(v, 2) for v in h) + " | " + " / ".join(pl(v, 2) for v in x) + f" | {j['E'] / s70['N'] - 1:+.3f} | {j['p']:.2g} |")
L += ["", "## W4. Domknięcie mb_nd wzdłuż skanu (b = N_fit/N_prawda − 1 ± HESSE; κ_MC; z FITBIAS/out/var_mc.json)", "",
      "| pT | model | 0,90 | 0,70 | 0,50 | 0,35 | 0,20 |", "|---|---|---|---|---|---|---|"]
for bn in ("3-5", "5-8", "8-30"):
    for m in ("M0", "M8"):
        row = []
        for wp in ("WP90", "WP70", "WP50", "WP35", "WP20"):
            o = MC[bn]["S"][wp]
            row.append(f"{o[m]['N'] / o['N_true'] - 1:+.2f} ± {(o[m]['err'] or float('nan')) / o['N_true']:.2f}")
            S.setdefault(f"mc|{bn}|{m}", {})[wp] = [o[m]["N"] / o["N_true"] - 1, (o[m]["err"] or float("nan")) / o["N_true"]]
        L.append(f"| {bn} | {m} | " + " | ".join(row) + " |")
open(os.path.join(HERE, "..", "tabele", "wpscan_m8.md"), "w").write("\n".join(L) + "\n")
json.dump(S, open(os.path.join(HERE, "..", "out", "wpscan_m8_summary.json"), "w"), indent=1)
print("\n".join(L))
