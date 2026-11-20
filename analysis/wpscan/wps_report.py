#!/usr/bin/env python3
"""WPSCAN — odczyt wg KRYTERIA.md z wps_fit.json → tabele markdown (stdout) + wps_summary.json.
  wps_report.py --fit wps_fit.json --out wps_summary.json [--scan]"""
import argparse
import json
import math
import sys

import numpy as np

SRCS = ["427929|MB427929", "428770|all"]
LO, HI, REF = 0.2, 0.95, 0.7


def f3(x, nd=3):
    return "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}".replace(".", ",")


def point(r, clos):
    """Wielkości pochodne jednego punktu."""
    o = {"eps_nom": r["eps_nom"], "eps": r["eps"], "ok": False, "why": []}
    n = r.get("nom")
    if n is None:
        o["why"].append("brak MC"); return o
    N, eh = n["N"], n["err"]
    bN = np.array([x if x is not None else np.nan for x in r["boot_N"]], float)
    bE = np.array([x if x is not None else np.nan for x in r["boot_err"]], float)
    o["bN"], o["bE"] = bN, bE
    nb = int(np.isfinite(bN).sum())
    sb = float(np.nanstd(bN)) if nb >= 10 else None
    o.update({"N": N, "err_hesse": eh, "err_boot": sb, "nboot_ok": nb, "chi2": n["chi2"], "ndf": n["ndf"], "p": n["p"], "maxpull": n["maxpull"],
              "sc": n["sc"], "dmu": n["dmu"], "frefl": r.get("frefl", 0.0), "sb_ratio": N / max(n["bkg_win"], 1e-9), "bkg_win": n["bkg_win"]})
    if not n["valid"] or eh is None:
        o["why"].append("fit nieważny")
    if n["at_limit"]:
        o["why"].append("granica dmu/sc")
    if N <= 0:
        o["why"].append("N ≤ 0")
    lowmu = n["n_lowmu"] > 0.2 * (n["ndf"] + 6)
    o["lowmu"] = bool(lowmu)
    if not lowmu and n["p"] < 1e-3:
        o["why"].append(f"p(χ²) = {n['p']:.1e}")
    if n["maxpull"] >= 4.5:
        o["why"].append(f"pull {n['maxpull']:.1f} @ {n['m_maxpull']:.3f}")
    ords = [v for k, v in r["var"].items() if k.startswith("exp") and v["valid"]]
    if not ords:
        o["why"].append("brak wariantów rzędu")
    if eh is None or N <= 0:
        return o
    o["syst_bkg"] = max(abs(v["N"] - N) for v in ords) if ords else None
    ext = [v for k, v in r["var"].items() if v["valid"]]
    o["syst_ext"] = max(abs(v["N"] - N) for v in ext) if ext else None
    o["stat"] = max(eh, sb or 0.0)
    o["F_stat"] = o["stat"] / N
    o["F"] = math.sqrt(o["stat"] ** 2 + (o["syst_bkg"] or 0.0) ** 2) / N
    o["F_ext"] = math.sqrt(o["stat"] ** 2 + (o["syst_ext"] or 0.0) ** 2) / N
    o["mcstat"] = math.sqrt(max(r["eps"] * (1 - r["eps"]), 0) / max(r["n_mc_tot"], 1)) / max(r["eps"], 1e-9)
    o["y"] = N / r["eps"]
    if clos is not None and clos.get("nom") and clos["nom"]["err"]:
        cn = clos["nom"]
        b = cn["N"] / max(clos["N_truth"], 1) - 1
        se = cn["err"] / max(clos["N_truth"], 1)
        o["clos_bias"], o["clos_err"] = b, se
        if cn["valid"] and abs(b) > 2.5 * se and abs(b) > 0.05:
            o["why"].append(f"domknięcie MC {b:+.2f} ± {se:.2f}")
    o["ok"] = not o["why"]
    return o


def analyse(P):
    """P: lista punktów (malejące ε). → słownik podsumowania klasyfikatora w komórce."""
    by = {p["eps_nom"]: p for p in P}
    grid = sorted(by)
    good = [e for e in grid if by[e]["ok"] and LO - 1e-9 <= e <= HI + 1e-9]
    S = {"n_good": len(good), "bad": {f"{e:.2f}": by[e]["why"] for e in grid if not by[e]["ok"] and LO - 1e-9 <= e <= HI + 1e-9}}
    if not by.get(REF, {}).get("ok") or len(good) < 8:
        S["undetermined"] = "WP70 nieadekwatny" if not by.get(REF, {}).get("ok") else f"adekwatnych punktów {len(good)} < 8"
        return S, by
    for e in good:                                   # wygładzanie: średnia z adekwatnych sąsiadów ±0,05
        nb = [by[x]["F"] for x in (round(e - 0.05, 2), e, round(e + 0.05, 2)) if x in by and by[x]["ok"] and LO - 1e-9 <= x <= HI + 1e-9]
        by[e]["Fs"] = float(np.mean(nb))
    eo = min(good, key=lambda e: by[e]["Fs"])
    po, pr = by[eo], by[REF]
    rel = lambda p: p["bE"] / p["bN"]  # noqa: E731
    gb = rel(po) / rel(pr)
    S.update({"eps_opt": eo, "F70": pr["F"], "Fs70": pr["Fs"], "Fopt": po["F"], "Fsopt": po["Fs"], "g": po["Fs"] / pr["Fs"], "g_raw": po["F"] / pr["F"],
              "g_err": float(np.nanstd(gb)) if eo != REF else 0.0, "g_stat_only": po["F_stat"] / pr["F_stat"],
              "F70_stat": pr["F_stat"], "Fopt_stat": po["F_stat"], "F70_bkg": (pr["syst_bkg"] or 0) / pr["N"], "Fopt_bkg": (po["syst_bkg"] or 0) / po["N"],
              "N70": pr["N"], "Nopt": po["N"], "relb_opt": rel(po).tolist(), "relb_70": rel(pr).tolist()})
    nbr = [round(eo - 0.05, 2), round(eo + 0.05, 2)]
    S["nbr"] = {f"{x:.2f}": (by[x]["F"] / pr["F"] if x in by and by[x].get("F") and by[x]["ok"] else None) for x in nbr}
    S["nbr_ok"] = all(v is not None and v < 0.95 for v in S["nbr"].values())
    # dryf N/ε względem WP70
    d, sd = {}, {}
    for e in good:
        p = by[e]
        d[e] = p["y"] / pr["y"] - 1
        rb = (p["bN"] / p["eps"]) / (pr["bN"] / pr["eps"]) - 1
        sd[e] = float(np.nanstd(rb)) if e != REF else 0.0
        p["d"], p["d_err"] = d[e], sd[e]
    oth = [e for e in good if e != REF]
    S["drift_rms"] = float(np.sqrt(np.mean([d[e] ** 2 for e in oth])))
    S["drift_chi2"] = float(sum((d[e] / max(sd[e], 1e-9)) ** 2 for e in oth)); S["drift_ndf"] = len(oth)
    S["drift_max"] = max(((abs(d[e]) / max(sd[e], 1e-9)), e, d[e], sd[e]) for e in oth)
    S["d_opt"], S["d_opt_err"] = d[eo], sd[eo]
    S["drift_sig_opt"] = bool(eo != REF and abs(d[eo]) > 2 * sd[eo])
    Feff = math.sqrt(po["Fs"] ** 2 + d[eo] ** 2) if S["drift_sig_opt"] else po["Fs"]
    S["g_eff"] = Feff / pr["Fs"]
    if S["g"] > 0.97:
        S["P1"] = "WP70 w optimum (płasko)"
    elif S["g_eff"] <= 0.90 and S["g_eff"] + 2 * S["g_err"] < 1 and S["nbr_ok"]:
        S["P1"] = "OPŁACA SIĘ"
    elif S["g"] <= 0.90:
        S["P1"] = "zysk niepotwierdzony (" + ("dryf zjada zysk" if S["g_eff"] > 0.90 else "sąsiedzi" if not S["nbr_ok"] else "2σ") + ")"
    else:
        S["P1"] = "mały zysk, bez rekomendacji"
    return S, by


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--out", default=""); ap.add_argument("--scan", action="store_true")
    a = ap.parse_args(argv)
    J = json.load(open(a.fit))
    info, R = J["info"], J["res"]
    idx = {(r["sp"], r["tag"], r["bin"], r["src"], r["j"]): r for r in R}
    bins = [f"{b[0]}-{b[1]}" for b in info["bins"]]
    nE = len(info["eps_desc"])
    SUM, PTS = {}, {}
    for sp, tags in info["tags"].items():
        for src in SRCS:
            for bn in bins:
                for tag in tags:
                    P = []
                    for j in range(nE):
                        r = idx.get((sp, tag, bn, src, j))
                        if r is None:
                            continue
                        P.append(point(r, idx.get((sp, tag, bn, "mb|all", j))))
                    if not P:
                        continue
                    S, by = analyse(P)
                    SUM[(sp, src, bn, tag)] = S; PTS[(sp, src, bn, tag)] = by
    kind = lambda t: "drzewo" if ".xgb." in t else "sieć"  # noqa: E731
    print("### T1. Optimum punktu pracy i zysk względem WP70 (F = √(stat² + tło²)/N; wygładzone 3 punkty)\n")
    print("| gatunek | próbka | pT | klasyf. | ε_opt | F(WP70) [stat / tło] | F(ε_opt) [stat / tło] | g = F_opt/F_70 | g surowe | g sam stat | sąsiedzi F/F70 | werdykt P1 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (sp, src, bn, tag), S in SUM.items():
        if "undetermined" in S:
            print(f"| {sp} | {src.split('|')[0]} | {bn} | {kind(tag)} | — | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE: {S['undetermined']} |")
            continue
        print(f"| {sp} | {src.split('|')[0]} | {bn} | {kind(tag)} | {f3(S['eps_opt'], 2)} | {f3(S['Fs70'])} [{f3(S['F70_stat'])} / {f3(S['F70_bkg'])}] | "
              f"{f3(S['Fsopt'])} [{f3(S['Fopt_stat'])} / {f3(S['Fopt_bkg'])}] | {f3(S['g'])} ± {f3(S['g_err'])} | {f3(S['g_raw'])} | {f3(S['g_stat_only'])} | "
              f"{' / '.join(f3(v, 2) for v in S['nbr'].values())} | {S['P1']} |")
    print("\n### T2. Stabilność N/ε_MC wzdłuż skanu (d = [N/ε]/[N/ε]_70 − 1; punkty adekwatne, ε ∈ [0,2; 0,95])\n")
    print("| gatunek | próbka | pT | klasyf. | rms d | χ²/ndf (d = 0) | największe odchylenie (ε: d ± σ) | d(ε_opt) | istotny w optimum | g_eff (z dryfem) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for (sp, src, bn, tag), S in SUM.items():
        if "undetermined" in S:
            continue
        mx = S["drift_max"]
        print(f"| {sp} | {src.split('|')[0]} | {bn} | {kind(tag)} | {f3(S['drift_rms'])} | {f3(S['drift_chi2'], 1)}/{S['drift_ndf']} | "
              f"{f3(mx[1], 2)}: {mx[2]:+.3f} ± {mx[3]:.3f} | {S['d_opt']:+.3f} ± {S['d_opt_err']:.3f} | {'TAK' if S['drift_sig_opt'] else 'nie'} | {f3(S['g_eff'])} |")
    print("\n### T3. Sieć w swoim optimum vs drzewo w swoim optimum (ρ = F_min(sieć)/F_min(drzewo))\n")
    print("| gatunek | próbka | pT | ε_opt drzewo / sieć | F_min drzewo | F_min sieć | ρ (optimum) | ρ przy WP70 | ρ sam stat (optimum) | werdykt P3 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    P3 = {}
    for sp, tags in info["tags"].items():
        for src in SRCS:
            for bn in bins:
                St, Sn = SUM.get((sp, src, bn, tags[0])), SUM.get((sp, src, bn, tags[1]))
                if St is None or Sn is None:
                    continue
                if "undetermined" in St or "undetermined" in Sn:
                    print(f"| {sp} | {src.split('|')[0]} | {bn} | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE |")
                    continue
                rho = Sn["Fsopt"] / St["Fsopt"]
                rb = np.array(Sn["relb_opt"]) / np.array(St["relb_opt"])
                se = float(np.nanstd(rb))
                rho70 = Sn["Fs70"] / St["Fs70"]
                v = "sieć realnie lepsza" if rho <= 0.90 and rho + 2 * se < 1 else ("zysk mały" if rho <= 1.0 else "sieć nie lepsza")
                P3[(sp, src, bn)] = {"rho": rho, "rho_err": se, "rho70": rho70, "verdict": v}
                print(f"| {sp} | {src.split('|')[0]} | {bn} | {f3(St['eps_opt'], 2)} / {f3(Sn['eps_opt'], 2)} | {f3(St['Fsopt'])} | {f3(Sn['Fsopt'])} | "
                      f"{f3(rho)} ± {f3(se)} | {f3(rho70)} | {f3(Sn['Fopt_stat'] / St['Fopt_stat'])} | {v} |")
    print("\n### T4. Adekwatność fitu — punkty odrzucone w ε ∈ [0,2; 0,95]\n")
    print("| gatunek | próbka | pT | klasyf. | adekwatnych / 16 | odrzucone (ε: powód) |")
    print("|---|---|---|---|---|---|")
    for (sp, src, bn, tag), S in SUM.items():
        bad = "; ".join(f"{e}: {', '.join(w)}" for e, w in sorted(S["bad"].items())) or "—"
        print(f"| {sp} | {src.split('|')[0]} | {bn} | {kind(tag)} | {S['n_good']} | {bad} |")
    if a.scan:
        for (sp, src, bn, tag), by in PTS.items():
            print(f"\n#### Skan: {sp} {src} {bn} GeV, {kind(tag)} ({tag})\n")
            print("| ε | ε_MC | N ± stat (HESSE / boot) | tło ±1 | tło rozsz. | F | F wygł. | N/ε | d vs WP70 | S/B okno | χ²/ndf | p | max pull | sc | f_refl | domkn. MC | adekw. |")
            print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for e in sorted(by, reverse=True):
                p = by[e]
                if "N" not in p:
                    print(f"| {e:.2f} | {f3(p['eps'])} | — |" + " — |" * 13 + f" {'; '.join(p['why'])} |"); continue
                cb = f"{p['clos_bias']:+.3f} ± {p['clos_err']:.3f}" if "clos_bias" in p else "—"
                dd = f"{p['d']:+.3f} ± {p['d_err']:.3f}" if "d" in p else "—"
                print(f"| {e:.2f} | {f3(p['eps'])} | {p['N']:.0f} ± {p.get('stat', float('nan')):.0f} ({(p['err_hesse'] or float('nan')):.0f} / {(p['err_boot'] or float('nan')):.0f}) | "
                      f"{(p.get('syst_bkg') or float('nan')):.0f} | {(p.get('syst_ext') or float('nan')):.0f} | {f3(p.get('F'))} | {f3(p.get('Fs'))} | {p.get('y', float('nan')):.0f} | {dd} | "
                      f"{f3(p['sb_ratio'])} | {p['chi2']:.0f}/{p['ndf']} | {p['p']:.3f} | {p['maxpull']:.1f} | {f3(p['sc'], 2)} | {f3(p['frefl'], 2)} | {cb} | {'tak' if p['ok'] else '; '.join(p['why'])} |")
    if a.out:
        clean = lambda S: {k: v for k, v in S.items() if k not in ("relb_opt", "relb_70")}  # noqa: E731
        json.dump({"cells": {"|".join(k): clean(S) for k, S in SUM.items()}, "P3": {"|".join(k): v for k, v in P3.items()}}, open(a.out, "w"), indent=1,
                  default=lambda o: o.item() if hasattr(o, "item") else str(o))


if __name__ == "__main__":
    sys.exit(main())
