#!/usr/bin/env python3
"""XSEC_atlas C — model generatora: zmiana C_ekskl przy przewazeniu widma pT w binie (nachylenie +-10 % na szerokosc binu; dla 8-30: w log pT) oraz
rozkladu |eta| (+-10 % dla |eta| > 1,1). Licznik wazony pT/eta reco, mianownik pT/eta prawdy. Tylko MC (mb_nd)."""
import glob, json, os, sys
import numpy as np
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(H, "skrypty"))
import xa_lib as X
import xa_fit as XF
T = {}
for f in sorted(glob.glob(os.path.join(H, "in/truth/xa_truth_*.npz"))):
    z = np.load(f)
    for k in ("pdg", "pt", "eta", "code", "mc", "d0code", "d0mc"):
        T.setdefault(k, []).append(z[k])
T = {k: np.concatenate(v) for k, v in T.items()}
ap = np.abs(T["pdg"])
EX = {"D0": (ap == 421) & (T["code"] == 101) & (T["mc"] <= 1), "Dplus": (ap == 411) & (T["code"] == 201) & (T["mc"] <= 1),
      "Dstar": (ap == 413) & (T["code"] == 100201) & (T["mc"] <= 1) & (T["d0code"] == 101) & (T["d0mc"] <= 1)}


def wts(pt, eta, lo, hi, kind, s):
    if kind == "pT":
        u = (np.log(pt / lo) / np.log(hi / lo) - 0.5)
        return 1 + s * 0.1 * 2 * u            # +-10 % na krancach binu
    return np.where(np.abs(eta) > 1.1, 1 + s * 0.1, 1.0)


out = {}; L = ["| gatunek | selekcja | bin | ΔC/C widmo pT ± [%] | ΔC/C η ± [%] | razem [%] |", "|---|---|---|---|---|---|"]
P = X.derived(X.load("mb_nd", "D0__pair_"), "D0")
for sp in ("D0", "Dplus", "Dstar"):
    MB = XF.load("mb_nd", sp) if sp != "Dstar" else None
    for bn, lo, hi in X.BINS:
        g = EX[sp] & (np.abs(T["eta"]) < 2.2) & (T["pt"] >= lo) & (T["pt"] < hi)
        if sp == "Dstar":
            m = (P["Lxy"] >= 0.25) & (np.abs(P["m"] - X.M0["D0"]) < 0.040) & (P["ptmin"] > 1.0) & (np.abs(P["eta"]) < 2.2) & (P["ds_pt"] >= lo) & (P["ds_pt"] < hi) & np.isin(P["ds_truth"], (0, 1))
            sels = {"Δm": (P["ds_pt"][m], P["eta"][m])}
        else:
            b = (MB["pt"] >= lo) & (MB["pt"] < hi) & (np.abs(MB["eta"]) < 2.2) & np.isin(MB["truth"], (0, 1)) & ((MB["swap"] == 0) if sp == "D0" else True)
            sels = {}
            for sn, s in (("WP70", XF.sel_wp(MB, sp, bn)), ("cięcia", XF.sel_cut(MB, sp))):
                sels[sn] = (MB["pt"][b & s], MB["eta"][b & s])
        for sn, (rpt, reta) in sels.items():
            d = {}
            for kind in ("pT", "eta"):
                r = [np.sum(wts(rpt, reta, lo, hi, kind, s)) / np.sum(wts(T["pt"][g], T["eta"][g], lo, hi, kind, s)) / (len(rpt) / g.sum()) - 1 for s in (+1, -1)]
                d[kind] = float(max(abs(r[0]), abs(r[1])))
            d["tot"] = float(np.hypot(d["pT"], d["eta"])); out[f"{sp}|{sn}|{bn}"] = d
            L.append(f"| {sp} | {sn} | {bn} | {100 * d['pT']:.1f} | {100 * d['eta']:.1f} | {100 * d['tot']:.1f} |")
print("\n".join(L))
json.dump(out, open(os.path.join(H, "out/gen_model.json"), "w"), indent=1)
open(os.path.join(H, "tabele/model_generatora.md"), "w").write("# Model generatora: czułość C_ekskl na widmo pT i η (mb_nd; przeważenie ±10 % na krańcach binu w log pT; ±10 % dla |η| > 1,1)\n\nZależność od krotności: nie do policzenia z tabel (brak krotności w prawdzie mianownika) — wiersz z literatury (δ5 = 2,0–2,5 %, 1512.02913).\n\n" + "\n".join(L) + "\n")
