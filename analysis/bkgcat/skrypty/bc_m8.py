#!/usr/bin/env python3
"""BKGCAT (2)+(3): szablon składowych pikujących (DMM + CM) wprost z drzewa H (bkgcat.npz) w formacie modelu M8 z FITBIAS, domknięcie na mb_nd i uzyski w danych.
Format (out/peak_templates_H.npz): '<gatunek>|<bin>|<kategoria>|h5' = histogram 5 MeV na pełnym zakresie (D0 1,70–2,05; D+ 1,70–2,20 — krawędzie fb_var.full_edges),
'…|kap' = N_(DMM+CM)/N_Sig (pełny zakres), '…|h5_pdg', '…|kap_pdg' = to samo po przeważeniu trybów do PDG (out/br.json), '…|nsig', '…|npk'.
Kategorie: WP95…WP10 (skumulowane) oraz B50, B40, B35 ([70, x)), L90-70 — nazwy jak fb_var.t_mc. Wczytanie do M8: k = fb_var.wcat(...); k['peakT'] = fb_fit.TmplE(full_edges, h5);
k['kap'] = kap (funkcja patch() niżej) — FITBIAS bez zmian.
Warianty: ref = FITBIAS (szablon T, κ interpolowane), dir = wprost z H, pdg = wprost + przeważenie BR; skala κ r = 1,364 (FITBIAS) albo 1 (BKGCAT: BR_MC = PDG).
POST HOC. Dane: tylko inkluzywne piki z WPSCAN (bez AFP). Wyjście: out/m8.json, tabele/m8.md.   bc_m8.py [--in bkgcat.npz] [--nproc 6]"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bc_lib as L  # noqa: E402
F = L.F
import fb_var as V  # noqa: E402

B = W = None
CATDEF = {f"WP{int(round(e * 100))}": (e, None) for e in [round(0.95 - 0.05 * i, 2) for i in range(18)]}
CATDEF.update({"B50": (0.7, 0.5), "B40": (0.7, 0.4), "B35": (0.7, 0.35), "L90-70": (0.9, 0.7)})


def _init(fn):
    global B, W
    V._init()
    B = L.BK(fn)
    W = json.load(open(os.path.join(L.BASE, "out", "br.json")))["weights"]


def mode_w(sp, key):
    """Waga BR_PDG/BR_MC dla klucza kategorii (DMM: tryb prawdziwego D gatunku kandydata; CM|411, CM|421: tryb tamtego hadronu; reszta 1)."""
    p = key.split("|")
    g = L.group(key)
    own = "421" if sp == "D0" else "411"
    if g in ("Sig", "Refl", "DMM"):
        return W.get(f"{own}|{p[1]}", 1.0)
    if g == "CM" and p[1] in ("411", "421"):
        return W.get(f"{p[1]}|{p[2]}", 1.0)
    return 1.0


def peak(sp, bn, elo, ehi=None):
    c = B.cat(sp, bn, elo, ehi)
    z = np.zeros(int(round((F.SP[sp]["mr"][1] - F.SP[sp]["mr"][0]) / 0.001)))
    pk, pkw, ns, nsw = z.copy(), z.copy(), 0.0, 0.0
    for k, h in c.items():
        g = L.group(k)
        if g in ("DMM", "CM"):
            pk += h; pkw += mode_w(sp, k) * h
        elif g == "Sig":
            ns += h.sum(); nsw += mode_w(sp, k) * h.sum()
    return {"h5": L.h5(pk), "h5_pdg": L.h5(pkw), "kap": pk.sum() / ns, "kap_pdg": pkw.sum() / nsw, "nsig": ns, "npk": pk.sum()}


def patch(k, sp, bn, elo, ehi, mode, r):
    """k z fb_var.wcat → wariant szablonu pikującego. mode: ref | dir | pdg; r = skala κ."""
    k = dict(k)
    if mode != "ref":
        p = peak(sp, bn, elo, ehi)
        k["peakT"] = F.TmplE(V.full_edges(sp), p["h5" if mode == "dir" else "h5_pdg"])
        k["kap"] = float(p["kap" if mode == "dir" else "kap_pdg"])
    k["kap"] = k["kap"] * r
    return k


RFB = {"D0": 1.364, "Dplus": 1.068}                 # skala κ użyta w FITBIAS (var_data.json)
VARS = [("ref", "fb"), ("ref", 1.0), ("dir", "fb"), ("dir", 1.0), ("pdg", 1.0)]


def task(t):
    sp, bn = t
    res = {"sp": sp, "bin": bn, "kap": {}, "mc": {}, "data": {}}
    k0 = V.wcat(sp, bn, 0.7)
    p = peak(sp, bn, 0.7)
    res["kap"] = {"ref": k0["kap"], "dir": float(p["kap"]), "pdg": float(p["kap_pdg"]), "nsig": float(p["nsig"]), "npk": float(p["npk"]),
                  "kap_err_dir": float(np.sqrt(p["npk"]) / p["nsig"])}
    # kształt: χ² szablonu T vs H (5 MeV → 25 MeV, zakres alt), znormalizowane
    rng = F.SP[sp]["alt"]
    e = np.arange(rng[0], rng[1] + F.BW / 2, F.BW)
    hT = k0["peakT"].bins(e, F.SP[sp]["m0"])
    hH = L.cutr(p["h5"], sp, rng)
    nb = (len(hH) // 5) * 5
    a, b = hH[:nb].reshape(-1, 5).sum(1), hT[:nb].reshape(-1, 5).sum(1)
    b = b / b.sum() * a.sum()
    res["shape_chi2"] = [float(np.sum((a - b) ** 2 / np.maximum(b, 1e-9))), int(len(a) - 1)]
    res["frac_win"] = {"dir": float(hH[np.abs(0.5 * (e[1:] + e[:-1]) - F.SP[sp]["m0"]) < F.SP[sp]["win"]].sum() / max(hH.sum(), 1)),
                       "ref": float(hT[np.abs(0.5 * (e[1:] + e[:-1]) - F.SP[sp]["m0"]) < F.SP[sp]["win"]].sum() / max(hT.sum(), 1e-9))}
    # domknięcie mb_nd (WP70 i kategorie luźne) — κ bez skali (MC = MC)
    for cn in ("WP70", "WP50", "WP90", "B50"):
        elo, ehi = CATDEF[cn]
        kk = V.wcat(sp, bn, elo, ehi)
        H = V.wdata(kk, "mb|all")
        nt = float(V.wdata(kk, "mb|sig").sum())
        o = {"N_true": nt}
        for mname in ("M0", "M4", "M8"):
            for mode in ("ref", "dir"):
                if mname == "M0" and mode == "dir":
                    continue
                f = V.sres(V.fit_set(H[None], [patch(kk, sp, bn, elo, ehi, mode, 1.0)], sp, mname))
                o[f"{mname}|{mode}"] = {"b": f["N"] / nt - 1, "b_err": (f["err"] or float("nan")) / nt, "valid": f["valid"], "sc": f["sc"], "chi2": f["chi2"], "ndf": f["ndf"]}
        res["mc"][cn] = o
    # dane WP70
    for run, src in V.SRC.items():
        D = V.wdata(k0, src)
        o = {}
        for mname in ("M0", "M4", "M8"):
            for mode, r in (VARS if mname != "M0" else [("ref", 1.0)]):
                f = V.sres(V.fit_set(D[0][None], [patch(k0, sp, bn, 0.7, None, mode, RFB[sp] if r == "fb" else r)], sp, mname))
                o[f"{mname}|{mode}|{r}"] = {"N": f["N"], "err": f["err"], "valid": f["valid"], "sc": f["sc"], "dmu": f["dmu"], "chi2": f["chi2"], "ndf": f["ndf"], "p": F.pval(f["chi2"], f["ndf"])}
        res["data"][run] = o
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=L.IN); ap.add_argument("--nproc", type=int, default=6)
    a = ap.parse_args()
    _init(a.inp)
    # 1) zapis szablonów w formacie M8
    out = {}
    for sp in L.SPS:
        out[f"edges|{sp}"] = V.full_edges(sp)
        for bn in L.BINS:
            for cn, (elo, ehi) in CATDEF.items():
                p = peak(sp, bn, elo, ehi)
                for kk, v in p.items():
                    out[f"{sp}|{bn}|{cn}|{kk}"] = np.asarray(v)
    np.savez_compressed(os.path.join(L.BASE, "out", "peak_templates_H.npz"), **out)
    with Pool(a.nproc, initializer=_init, initargs=(a.inp,)) as pool:
        R = pool.map(task, [(sp, bn) for sp in L.SPS for bn in L.BINS])
    json.dump(R, open(os.path.join(L.BASE, "out", "m8.json"), "w"))
    D = {(r["sp"], r["bin"]): r for r in R}
    fb = {(r["sp"], r["bin"], r["run"]): r for r in json.load(open(os.path.join(L.CW, "FITBIAS_2026-11-19", "out", "var_data.json")))["res"]}
    T = ["# BKGCAT — szablon składowych pikujących wprost z drzewa H w modelu M8 (POST HOC względem FITBIAS)", "",
         "ref = FITBIAS (szablon z MC_templates, pokrewny klasyfikator, WP80; κ interpolowane w ε); dir = DMM + CM z drzewa H przy tym samym WP (kształt i κ);",
         "pdg = dir + przeważenie trybów do PDG (`br.md`). r = skala κ w danych (FITBIAS: r_FB = 1,364 dla D⁰ i 1,068 dla D⁺; BKGCAT: 1 — BR trybu sygnału w MC = PDG). Błędy HESSE (FITBIAS: rzeczywisty rozrzut ~1,3×).", "",
         "## κ = N(DMM + CM)/N(Sig) przy WP70 i kształt", "", "| gat. | bin | κ ref (interp.) | κ dir ± stat | κ pdg | N szablonu H | udział w oknie: ref | dir | χ²/ndf kształtu T vs H (25 MeV) |", "|---|---|---|---|---|---|---|---|---|"]
    for (sp, bn), r in D.items():
        k = r["kap"]
        T.append(f"| {sp} | {bn} | {k['ref']:.3f} | {k['dir']:.3f} ± {k['kap_err_dir']:.3f} | {k['pdg']:.3f} | {k['npk']:.0f} | {r['frac_win']['ref']:.3f} | {r['frac_win']['dir']:.3f} | {r['shape_chi2'][0]:.1f}/{r['shape_chi2'][1]} |")
    T += ["", "## Domknięcie na mb_nd: b = N_fit/N_prawda − 1 [%]", "", "| gat. | bin | kategoria | M0 | M4 ref | M4 dir | M8 ref | M8 dir | sc (M8 dir) | p(χ²) M8 dir |", "|---|---|---|---|---|---|---|---|---|---|"]
    for (sp, bn), r in D.items():
        for cn, o in r["mc"].items():
            g = lambda k: L.pm(o[k]["b"], o[k]["b_err"], 1, True) + ("" if o[k]["valid"] else " (!)")  # noqa: E731
            T.append(f"| {sp} | {bn} | {cn} | {g('M0|ref')} | {g('M4|ref')} | {g('M4|dir')} | {g('M8|ref')} | {g('M8|dir')} | {o['M8|dir']['sc']:.2f} | {F.pval(o['M8|dir']['chi2'], o['M8|dir']['ndf']):.2f} |")
    T += ["", "## Dane, WP70: uzysk N (M8) per wariant szablonu", "",
          "| gat. | run | bin | FITBIAS M8 (zapis) | M8 ref r=r_FB | M8 ref r=1 | M8 dir r=r_FB | **M8 dir r=1** | M8 pdg r=1 | Δ(dir r=1 − ref r=r_FB) | Δ(dir − ref, r=1) | Δ(pdg − dir) | p(χ²): ref r_FB / dir 1 / pdg | M0 |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (sp, bn), r in D.items():
        for run, o in r["data"].items():
            g = lambda k: f"{o[k]['N']:.0f} ± {o[k]['err']:.0f}" if o[k]["err"] else f"{o[k]['N']:.0f} (!)"  # noqa: E731
            n = lambda k: o[k]["N"]  # noqa: E731
            ref = fb[(sp, bn, run)]["S"]["WP70"]["M8"]["N"]
            T.append(f"| {sp} | {run} | {bn} | {ref:.0f} | {g('M8|ref|fb')} | {g('M8|ref|1.0')} | {g('M8|dir|fb')} | **{g('M8|dir|1.0')}** | {g('M8|pdg|1.0')} | "
                     f"{100 * (n('M8|dir|1.0') / n('M8|ref|fb') - 1):+.1f} % | {100 * (n('M8|dir|1.0') / n('M8|ref|1.0') - 1):+.1f} % | {100 * (n('M8|pdg|1.0') / n('M8|dir|1.0') - 1):+.1f} % | "
                     f"{o['M8|ref|fb']['p']:.2g} / {o['M8|dir|1.0']['p']:.2g} / {o['M8|pdg|1.0']['p']:.2g} | {g('M0|ref|1.0')} |")
    T += ["", "## Dane, WP70: to samo dla M4 (zakres nominalny)", "", "| gat. | run | bin | M4 ref r=r_FB | M4 ref r=1 | M4 dir r=1 | M4 pdg r=1 | Δ(dir − ref, r=1) |", "|---|---|---|---|---|---|---|---|"]
    for (sp, bn), r in D.items():
        for run, o in r["data"].items():
            g = lambda k: f"{o[k]['N']:.0f} ± {o[k]['err']:.0f}" if o[k]["err"] else f"{o[k]['N']:.0f} (!)"  # noqa: E731
            T.append(f"| {sp} | {run} | {bn} | {g('M4|ref|fb')} | {g('M4|ref|1.0')} | {g('M4|dir|1.0')} | {g('M4|pdg|1.0')} | {100 * (o['M4|dir|1.0']['N'] / o['M4|ref|1.0']['N'] - 1):+.1f} % |")
    open(os.path.join(L.BASE, "tabele", "m8.md"), "w").write("\n".join(T) + "\n")
    print("zapisano out/peak_templates_H.npz, out/m8.json, tabele/m8.md")


if __name__ == "__main__":
    main()
