#!/usr/bin/env python3
"""FITBIAS — rozkład obciążenia na przyczyny (z out/diag_T.json): przyrosty b przy dokładaniu składowych (kształt stały), wpływ wolnego kształtu, obciążenie własne (zabawki).
→ tabele/decomp.md"""
import json, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "..", "out", "diag_T.json")))["res"]
f = lambda x, n=3: "—" if x is None or not math.isfinite(x) else f"{x:+.{n}f}".replace(".", ",")
u = lambda x, n=3: "—" if x is None or not math.isfinite(x) else f"{x:.{n}f}".replace(".", ",")
L = ["# FITBIAS — rozkład obciążenia na przyczyny (MC mb_nd z kategoriami prawdy, zestaw T). Wszystko w jednostkach N_prawda.", "",
     "Przyrost Δ = b(po dołożeniu składowej) − b(przed), kolejność I, model nominalny z kształtem STAŁYM (reżim addytywny); błąd przyrostu √|σ_k² − σ_(k−1)²|/N.",
     "„wolny − stały” = b(wszystko, wolny) − b(wszystko, stały) na tej samej realizacji; „własne” = średnie obciążenie w zabawkach z modelu (wolny kształt; stały ma 0,00 ± 0,01).", "",
     "| gat. | pT | kat. | S/B | szablon sygnału+odbicie (Sig+Refl) | +DMM | +CM | +PartComb | +Comb | suma (stały) ± σ_fit | wolny − stały | własne (zabawki, wolny) | DMM przy wolnym kształcie |", "|---|" + "---|" * 12]
agg = {}
for r in R:
    nt = r["N_true"]
    rows = {x["add"]: x for x in r["seq"]["I|fix"]}
    rf = {x["add"]: x for x in r["seq"]["I|free"]}
    order = [g for g in ("Sig", "Refl", "DMM", "CM", "PC", "Comb") if g in rows]
    base = "Refl" if "Refl" in rows else "Sig"
    inc, prev = {}, None
    for g in order:
        if prev is not None:
            e = math.sqrt(abs((rows[g]["err"] or 0) ** 2 - (rows[prev]["err"] or 0) ** 2)) / nt
            inc[g] = (rows[g]["b"] - rows[prev]["b"], e)
        prev = g
    b0 = rows[base]["b"]
    tot, tote = rows["Comb"]["b"], (rows["Comb"]["err"] or float("nan")) / nt
    dfree = rf["Comb"]["b"] - rows["Comb"]["b"]
    own, owne = r["toys"]["free"]["bias"], r["toys"]["free"]["bias_err"]
    dmm_free = rf["DMM"]["b"] - rf[base]["b"]
    L.append(f"| {r['sp']} | {r['bin']} | {r['cat']} | {u(r['SB_win'])} | {f(b0)} | " + " | ".join(f"{f(inc[g][0])} ± {u(inc[g][1])}" for g in ("DMM", "CM", "PC", "Comb")) + f" | {f(tot)} ± {u(tote)} | {f(dfree)} | {f(own)} ± {u(owne)} | {f(dmm_free)} |")
    agg.setdefault((r["sp"], r["cat"]), []).append((b0, inc, tot, tote, dfree, own, owne, dmm_free))
L += ["", "## Średnie po trzech binach pT (ważone 1/σ² dla przyrostów z błędem; arytmetyczne dla pozostałych)", "",
      "| gat. | kat. | Sig+Refl | +DMM | +CM | +PartComb | +Comb | suma (stały) | wolny − stały (średnia; zakres) | własne (wolny) | DMM przy wolnym |", "|---|" + "---|" * 10]
def wm(v, e):
    v, e = np.array(v), np.maximum(np.array(e), 1e-3); w = 1 / e ** 2
    return (v * w).sum() / w.sum(), 1 / math.sqrt(w.sum())
for (sp, cat), xs in agg.items():
    cells = []
    for g in ("DMM", "CM", "PC", "Comb"):
        m, e = wm([x[1][g][0] for x in xs], [x[1][g][1] for x in xs]); cells.append(f"{f(m)} ± {u(e)}")
    m, e = wm([x[2] for x in xs], [x[3] for x in xs])
    df = [x[4] for x in xs]
    om, oe = wm([x[5] for x in xs], [x[6] for x in xs])
    L.append(f"| {sp} | {cat} | {f(float(np.mean([x[0] for x in xs])))} | " + " | ".join(cells) + f" | {f(m)} ± {u(e)} | {f(float(np.mean(df)))} ({f(min(df))}…{f(max(df))}) | {f(om)} ± {u(oe)} | {f(float(np.mean([x[7] for x in xs])))} |")
open(os.path.join(HERE, "..", "tabele", "decomp.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L[-14:]))
