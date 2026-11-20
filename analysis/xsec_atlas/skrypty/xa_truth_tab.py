#!/usr/bin/env python3
"""XSEC_atlas — N_gen po PRAWDZIWYM trybie rozpadu (lista corek) vs po decayType; tabela per bin pT, |eta| < 2,2. Tylko MC (mb_nd)."""
import glob, json, os
import numpy as np
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = {}
for f in sorted(glob.glob(os.path.join(H, "in/truth/xa_truth_*.npz"))):
    z = np.load(f)
    for k in z.files:
        T.setdefault(k, []).append(z[k])
nev = int(sum(x[0] for x in T.pop("nev"))); T.pop("files"); T = {k: np.concatenate(v) for k, v in T.items()}
ap = np.abs(T["pdg"]); fid = np.abs(T["eta"]) < 2.2
BINS = [("3-5", 3, 5), ("3.5-5", 3.5, 5), ("5-8", 5, 8), ("8-30", 8, 30), ("3.5-30", 3.5, 30)]
# D*: lista corek zawiera D0 ORAZ jego corki -> D0 K- pi+ pi+ = 100201
# kody topologii: K- + 10 K+ + 100 pi+ + 1000 pi- + 10000 pi0 + 100000 D0 (znaki wzgledem rodzica)
EX = {"D0": (ap == 421) & (T["code"] == 101) & (T["mc"] <= 1), "Dplus": (ap == 411) & (T["code"] == 201) & (T["mc"] <= 1),
      "Dstar": (ap == 413) & (T["code"] == 100201) & (T["mc"] <= 1) & (T["d0code"] == 101) & (T["d0mc"] <= 1)}
DT = {"D0": (ap == 421) & (T["dty"] == 10), "Dplus": (ap == 411) & (T["dty"] == 4), "Dstar": (ap == 413) & (T["dty"] == 1)}
ALL = {"D0": ap == 421, "Dplus": ap == 411, "Dstar": ap == 413}
out = {"nev": nev, "bins": {}}
print("zdarzen", nev, "wierszy", len(ap))
print("decayType D*:", np.unique(T["dty"][ap == 413], return_counts=True))
print("D* -> D0 pi+ bez dowiazania corki D0:", int(((ap == 413) & (T["dty"] == 1) & (T["d0code"] < 0) & fid & (T["pt"] > 3)).sum()), "z", int(((ap == 413) & (T["dty"] == 1) & fid & (T["pt"] > 3)).sum()))
L = ["| gatunek | bin pT | N wszystkie | N decayType | N ekskl. | ekskl./decayType | BR_MC ekskl. [%] | f_b ekskl. |", "|---|---|---|---|---|---|---|---|"]
for sp in ("D0", "Dplus", "Dstar"):
    for bn, lo, hi in BINS:
        b = fid & (T["pt"] >= lo) & (T["pt"] < hi)
        na, nd, ne = int((ALL[sp] & b).sum()), int((DT[sp] & b).sum()), int((EX[sp] & b).sum())
        nrad = int((EX[sp] & b & (T["mc"] == 1)).sum()) if sp != "Dstar" else int((EX[sp] & b & (T["d0mc"] == 1)).sum())
        fb = float(T["fromB"][EX[sp] & b].mean())
        r = ne / nd; er = np.sqrt(r * (1 - r) / nd)
        out["bins"][f"{sp}|{bn}"] = {"N_all": na, "N_dty": nd, "N_excl": ne, "N_excl_rad": nrad, "f_b": fb, "ratio": r, "ratio_err": er}
        L.append(f"| {sp} | {bn} | {na} | {nd} | {ne} | {r:.4f} ± {er:.4f} | {100 * ne / na:.3f} ± {100 * np.sqrt(ne) / na:.3f} | {fb:.3f} |")
print("\n".join(L))
json.dump(out, open(os.path.join(H, "out/truth_ngen.json"), "w"), indent=1)
open(os.path.join(H, "tabele/ngen_prawda.md"), "w").write("# N_gen w mb_nd (|η| < 2,2): po prawdziwym trybie rozpadu (lista córek) vs po `decayType`\n\n"
    f"mb_nd: {nev} zdarzeń. Ekskl. = dokładnie K⁻π⁺ / K⁻π⁺π⁺ / D⁰(→K⁻π⁺)π⁺ (znaki względem rodzica), brak brakującej masy (m_D − m_wid < 30 MeV) lub radiacyjny (m_brak < 100 MeV). D*: lista córek D* zawiera D⁰ i jego córki (kod D⁰ K⁻ π⁺ π⁺), `decayType` 1.\n\n" + "\n".join(L) + "\n")
