#!/usr/bin/env python3
"""WPTIGHT — odczyt wg KRYTERIA.md: wpt_mc.json + wpt_data.json → tabele markdown (stdout) + wpt_summary.json.
  wpt_report.py --mc wpt_mc.json --data wpt_data.json --out wpt_summary.json > wpt_tables.md"""
import argparse
import json
import math
import sys

import numpy as np

G2_MAX, G1_NS = 0.05, 2.5
SF70_SCEN = [0.0, 0.05, 0.10, 0.15]
MAIN = ["427929|MB427929", "428770|all"]


def F(x, nd=3):
    return "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}".replace(".", ",")


def pm(x, e, nd=3):
    return f"{F(x, nd)} ± {F(e, nd)}"


def arr(bl, k):
    return np.array([b[k] if b is not None and b.get(k) is not None else np.nan for b in bl], float)


def rms(x):
    x = x[np.isfinite(x)]
    return float(np.std(x)) if len(x) >= 10 else None


def jk_ratio_rms(a, b):
    """Jackknife błędu stosunku rms(a)/rms(b) po replikach (te same repliki)."""
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = len(a)
    if n < 10:
        return None
    th = np.array([np.std(np.delete(a, i)) / np.std(np.delete(b, i)) for i in range(n)])
    return float(math.sqrt((n - 1) / n * np.sum((th - th.mean()) ** 2)))


def ordvars(var):
    return [v for k, v in var.items() if k.startswith("exp") and v.get("valid")]


def analyse(r):
    """Wielkości pochodne jednego zadania (dane lub MC)."""
    j, s70 = r["joint"], r["wp70"]
    o = {"sp": r["sp"], "bin": r["bin"], "src": r["src"], "eps_t": r["eps_t"], "f_mc": r["f_mc"], "why": []}
    if not j["valid"] or "E_err" not in j:
        o["why"].append("fit łączny nieważny"); return o
    if j["at_limit"]:
        o["why"].append("granica dmu/sc (łączny)")
    if not s70["valid"] or s70["err"] is None:
        o["why"].append("fit WP70 nieważny"); return o
    if s70["at_limit"]:
        o["why"].append("granica dmu/sc (WP70)")
    lowmu = r["lowmu_frac"] > 0.2
    if not lowmu and j["p"] < 1e-3:
        o["why"].append(f"p(χ²) łączny {j['p']:.1e}")
    if not lowmu and s70["p"] < 1e-3:
        o["why"].append(f"p(χ²) WP70 {s70['p']:.1e}")
    bj, b7, bt = r["boot_joint"], r["boot_wp70"], r["boot_tight"]
    bf, bE, bEe, bNA = arr(bj, "f"), arr(bj, "E"), arr(bj, "E_err"), arr(bj, "NA")
    b7N, b7e = arr(b7, "N"), arr(b7, "err")
    btN = arr(bt, "N")
    fm = r["f_mc"]
    # SF_rel
    sf = j["f"] / fm
    sf_h, sf_b = j["f_err"] / fm, (rms(bf) or 0.0) / fm
    sf_stat = math.hypot(max(sf_h, sf_b), sf * r["f_mc_err"] / fm)
    ov = ordvars(r["joint_var"])
    if not ov:
        o["why"].append("brak wariantów rzędu (łączny)")
    sf_bkg = max(abs(v["f"] / fm - sf) for v in ov) if ov else None
    E = j["E"]
    E_stat = max(j["E_err"], rms(bE) or 0.0)
    E_bkg = max(abs(v["E"] - E) for v in ov) if ov else None
    ov7 = ordvars(r["wp70_var"])
    if not ov7:
        o["why"].append("brak wariantów rzędu (WP70)")
    N7 = s70["N"]
    N7_stat = max(s70["err"], rms(b7N) or 0.0)
    N7_bkg = max(abs(v["N"] - N7) for v in ov7) if ov7 else None
    st_t = math.hypot(E_stat, E_bkg or 0.0) / E
    st_7 = math.hypot(N7_stat, N7_bkg or 0.0) / N7
    g = st_t / st_7
    rr = (bEe / bE) / (b7e / b7N)
    sg1 = rms(rr)
    sg2 = jk_ratio_rms(bE / E, b7N / N7)
    g_h = (j["E_err"] / E) / (s70["err"] / N7)
    g_b = ((rms(bE) or float("nan")) / E) / ((rms(b7N) or float("nan")) / N7)
    rel = [x for x in ((sg1 / g_h) if sg1 else None, (sg2 / g_b) if sg2 and math.isfinite(g_b) else None) if x is not None]
    sg = g * max(rel) if rel else None                   # błąd względny (większy z dwóch oszacowań) przeniesiony na g
    dl = E / N7 - 1
    dl_e = rms(bE / b7N)
    # rozkład
    NA = j["NA"]
    NA_rel = max(j["NA_err"], rms(bNA) or 0.0) / NA
    sf_rel_tot = math.hypot(sf_stat, sf_bkg or 0.0) / sf
    okb = np.isfinite(bNA) & np.isfinite(bf)
    corr_b = float(np.corrcoef(bNA[okb], bf[okb])[0, 1]) if okb.sum() >= 10 else None
    # dryf WPSCAN (dwa niezależne fity jednej kategorii)
    t = r["tight"]
    d = (t["N"] / N7) / fm - 1 if t["valid"] else None
    d_e = (rms(btN / b7N) or float("nan")) / fm if t["valid"] else None
    o.update({"ok": not o["why"], "sf": sf, "sf_stat": sf_stat, "sf_hesse": sf_h, "sf_boot": sf_b, "sf_bkg": sf_bkg, "sf_tot": math.hypot(sf_stat, sf_bkg or 0.0),
              "E": E, "E_stat": E_stat, "E_hesse": j["E_err"], "E_boot": rms(bE), "E_bkg": E_bkg, "tot_t": st_t, "stat_t": E_stat / E, "bkg_t": (E_bkg or 0.0) / E,
              "N70": N7, "N70_stat": N7_stat, "N70_hesse": s70["err"], "N70_boot": rms(b7N), "N70_bkg": N7_bkg, "tot_70": st_7, "stat_70": N7_stat / N7,
              "bkg_70": (N7_bkg or 0.0) / N7, "g": g, "g_err": sg, "g_stat": (E_stat / E) / (N7_stat / N7), "g_hesse": (j["E_err"] / E) / (s70["err"] / N7),
              "delta": dl, "delta_err": dl_e, "NA": NA, "NA_rel": NA_rel, "sf_rel_tot": sf_rel_tot, "naive": math.hypot(NA_rel, sf_rel_tot),
              "corr_h": j.get("corr_NA_f"), "corr_b": corr_b, "N_tight_single": t["N"], "tight_rel": (max(t["err"] or 0.0, rms(btN) or 0.0) / t["N"]) if t["valid"] and t["N"] > 0 else None,
              "drift": d, "drift_err": d_e, "p_joint": j["p"], "p_70": s70["p"], "sc": j["sc"], "sc70": s70["sc"], "dmu": j["dmu"], "lowmu": lowmu,
              "sb_A": NA * 1.0 / max(j["bkg_win"][0], 1e-9), "sb_B": j["NB"] / max(j["bkg_win"][1], 1e-9), "sb_70": N7 / max(s70["bkg_win"], 1e-9),
              "var": {k: (v["f"] / fm if v.get("valid") and v.get("f") else None) for k, v in r["joint_var"].items()},
              "varE": {k: (v["E"] / E - 1 if v.get("valid") and v.get("E") else None) for k, v in r["joint_var"].items()}})
    if "truth" in r:
        tr = r["truth"]
        T7 = tr["NA"] + tr["NB"]
        o.update({"sf_truth": tr["NA"] / T7 / fm, "sf_truth_err": math.sqrt(tr["NA"] * tr["NB"] / T7 ** 3) / fm, "T70": T7, "E_clos": E / T7, "E_clos_err": E_stat / T7,
                  "N70_clos": N7 / T7, "N70_clos_err": N7_stat / T7})
    return o


def verdict(a, gate):
    if not gate["pass"]:
        return "NIEROZSTRZYGNIĘTE (domknięcie: " + gate["why"] + ")"
    if not a.get("ok"):
        return "NIEROZSTRZYGNIĘTE (fit: " + "; ".join(a["why"]) + ")"
    if a["delta_err"] and abs(a["delta"]) > 2.5 * a["delta_err"]:
        return "NIEROZSTRZYGNIĘTE (E ≠ N70)"
    g, sg = a["g"], a["g_err"] or 0.0
    if g <= 0.90 and g + 2 * sg < 1:
        return "OPŁACA SIĘ"
    if g <= 0.90:
        return "zysk niepotwierdzony (2σ)"
    if g <= 1.0:
        return "nie opłaca się (< 10 %)"
    return "szkodzi"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", required=True); ap.add_argument("--data", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    MC = {(r["sp"], r["bin"], r["eps_t"]): analyse(r) for r in json.load(open(a.mc))["res"]}
    DA = {(r["sp"], r["bin"], r["src"], r["eps_t"]): analyse(r) for r in json.load(open(a.data))["res"]}
    gate = {}
    P = print
    P("# WPTIGHT — tabele (generowane: wpt_report.py)\n")
    P("## M1. Domknięcie na MC (pseudo-dane mb_nd, szablony i f_MC z fd_nd) — bramka G1 (|SF−1| ≤ 2,5σ) i G2 (σ ≤ 0,05)\n")
    P("| gatunek | pT | ε | f_MC | SF_rel^MC (fit) ± σ [HESSE / repliki] | SF_rel^prawda (zliczanie) | pull | E/N_prawda | N_70/N_prawda | S/B okno A / B / WP70 | sc | p(χ²) | bramka |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k, m in MC.items():
        if "sf" not in m:
            gate[k] = {"pass": False, "why": "fit MC nieważny"}; P(f"| {k[0]} | {k[1]} | {F(k[2], 2)} | — | fit nieważny | | | | | | | | NIE |"); continue
        s = max(m["sf_hesse"], m["sf_boot"])
        why = []
        if abs(m["sf"] - 1) > G1_NS * s:
            why.append(f"G1 {abs(m['sf'] - 1) / s:.1f}σ")
        if s > G2_MAX:
            why.append(f"G2 σ = {F(s)}")
        if m["why"]:
            why.append("; ".join(m["why"]))
        gate[k] = {"pass": not why, "why": ", ".join(why), "sf": m["sf"], "err": s, "sf_truth": m["sf_truth"]}
        P(f"| {k[0]} | {k[1]} | {F(k[2], 2)} | {F(m['f_mc'])} | {F(m['sf'])} ± {F(s)} [{F(m['sf_hesse'])} / {F(m['sf_boot'])}] | {pm(m['sf_truth'], m['sf_truth_err'])} | {(m['sf'] - 1) / s:+.1f} | "
          f"{pm(m['E_clos'], m['E_clos_err'])} | {pm(m['N70_clos'], m['N70_clos_err'])} | {F(m['sb_A'], 2)} / {F(m['sb_B'], 2)} / {F(m['sb_70'], 2)} | {F(m['sc'], 2)} | {F(m['p_joint'])} | "
          f"{'TAK' if gate[k]['pass'] else 'NIE: ' + gate[k]['why']} |")
    P("\n## D1. Dane: SF_rel i bilans błędu uzysku (E = N_A + N_B z fitu łącznego ≡ N_A/(f_MC·SF_rel); N_70 = fit jednej kategorii)\n")
    P("| gatunek | próbka | pT | ε | SF_rel ± stat ± tło | σ_tot(ciasny) [stat / tło] | σ_tot(WP70) [stat / tło] | g ± σ_g | g sam stat (HESSE) | Δ = E/N_70 − 1 | werdykt |")
    P("|---|---|---|---|---|---|---|---|---|---|---|")
    summ = []
    for k, d in DA.items():
        if k[2] not in MAIN:
            continue
        gt = gate[(k[0], k[1], k[3])]
        if "sf" not in d:
            P(f"| {k[0]} | {k[2].split('|')[0]} | {k[1]} | {F(k[3], 2)} | — | — | — | — | — | — | NIEROZSTRZYGNIĘTE (fit: {'; '.join(d['why'])}) |")
            summ.append({"key": k, "verdict": "NIEROZSTRZYGNIĘTE (fit)"}); continue
        v = verdict(d, gt)
        d["verdict"] = v
        P(f"| {k[0]} | {k[2].split('|')[0]} | {k[1]} | {F(k[3], 2)} | {F(d['sf'])} ± {F(d['sf_stat'])} ± {F(d['sf_bkg'])} | {F(d['tot_t'])} [{F(d['stat_t'])} / {F(d['bkg_t'])}] | "
          f"{F(d['tot_70'])} [{F(d['stat_70'])} / {F(d['bkg_70'])}] | {pm(d['g'], d['g_err'])} | {F(d['g_stat'])} ({F(d['g_hesse'])}) | {d['delta']:+.3f} ± {F(d['delta_err'])} | {v} |")
        summ.append({"key": k, **{x: d[x] for x in d if x not in ("why",)}, "gate": gt})
    P("\n## D2. Rozkład błędu ciasnego punktu: σ(N_A)/N_A, σ(SF_rel)/SF_rel, korelacja; suma naiwna ⊕ wobec poprawnej (z korelacją)\n")
    P("| gatunek | próbka | pT | ε | N_A (łączny) | σ(N_A)/N_A | σ(SF_rel)/SF_rel (stat ⊕ tło) | corr(N_A, f) HESSE / repliki | naiwne ⊕ | poprawne σ_tot(E) | N_ciasny (1 kat.) ± rel | S/B okno A / B / WP70 | sc łączny / WP70 | p łączny / WP70 |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k, d in DA.items():
        if k[2] not in MAIN or "sf" not in d:
            continue
        P(f"| {k[0]} | {k[2].split('|')[0]} | {k[1]} | {F(k[3], 2)} | {F(d['NA'], 0)} | {F(d['NA_rel'])} | {F(d['sf_rel_tot'])} | {F(d['corr_h'], 2)} / {F(d['corr_b'], 2)} | {F(d['naive'])} | {F(d['tot_t'])} | "
          f"{F(d['N_tight_single'], 0)} ± {F(d['tight_rel'])} | {F(d['sb_A'], 2)} / {F(d['sb_B'], 2)} / {F(d['sb_70'], 2)} | {F(d['sc'], 2)} / {F(d['sc70'], 2)} | {F(d['p_joint'])} / {F(d['p_70'])} |")
    P("\n## D3. SF_rel: fit łączny wobec dryfu WPSCAN (dwa fity jednej kategorii), wobec domknięcia MC; stosunek podwójny dane/MC (opisowo)\n")
    P("| gatunek | próbka | pT | ε | SF_rel − 1 (łączny) ± tot | d WPSCAN (N_c/N_70/f_MC − 1) ± repliki | SF_rel^MC (fit) | SF_rel dane / SF_rel^MC | warianty SF_rel: rząd−1 / rząd+1 / pol / zakres / wspólny szablon |")
    P("|---|---|---|---|---|---|---|---|---|")
    for k, d in DA.items():
        if k[2] not in MAIN or "sf" not in d:
            continue
        gt = gate[(k[0], k[1], k[3])]
        dr = d["sf"] / gt["sf"]
        dre = dr * math.hypot(d["sf_tot"] / d["sf"], gt["err"] / gt["sf"])
        vv = " / ".join(F(x) for x in d["var"].values())
        d["dratio"], d["dratio_err"] = dr, dre
        P(f"| {k[0]} | {k[2].split('|')[0]} | {k[1]} | {F(k[3], 2)} | {d['sf'] - 1:+.3f} ± {F(d['sf_tot'])} | {(d['drift'] if d['drift'] is not None else float('nan')):+.3f} ± {F(d['drift_err'])} | "
          f"{pm(gt['sf'], gt['err'])} | {pm(dr, dre)} | {vv} |")
    P("\n## D4. 428770: SF_rel per klasa triggera (opisowo; klasy rozłączne ⇒ niezależne)\n")
    P("| gatunek | pT | ε | all | J12only | notJ12only | J12only − notJ12only [σ] | bramka MC |")
    P("|---|---|---|---|---|---|---|---|")
    cls = []
    for (sp, bn, e), gt in gate.items():
        row = [DA.get((sp, bn, f"428770|{c}", e)) for c in ("all", "J12only", "notJ12only")]
        tx = [pm(x["sf"], x["sf_tot"]) + ("" if x.get("ok") else " (!)") if x and "sf" in x else "—" for x in row]
        z = None
        if all(x and "sf" in x for x in row[1:]):
            z = (row[1]["sf"] - row[2]["sf"]) / math.hypot(row[1]["sf_tot"], row[2]["sf_tot"])
        cls.append({"sp": sp, "bin": bn, "eps": e, "z": z, "sf": [x.get("sf") if x else None for x in row], "err": [x.get("sf_tot") if x else None for x in row]})
        P(f"| {sp} | {bn} | {F(e, 2)} | {tx[0]} | {tx[1]} | {tx[2]} | {'—' if z is None else f'{z:+.1f}'} | {'TAK' if gt['pass'] else 'NIE'} |")
    P("\n## D5. Bilans bezwzględny: σ_abs = σ_tot ⊕ σ[SF_70]; g_abs = σ_abs(ciasny)/σ_abs(WP70) dla σ[SF_70] = 0 / 5 / 10 / 15 %\n")
    P("| gatunek | próbka | pT | ε | σ_tot ciasny / WP70 | g_abs (0 %) | (5 %) | (10 %) | (15 %) | σ_abs ciasny → WP70 przy 10 % |")
    P("|---|---|---|---|---|---|---|---|---|---|")
    for k, d in DA.items():
        if k[2] not in MAIN or "sf" not in d:
            continue
        ga = [math.hypot(d["tot_t"], s) / math.hypot(d["tot_70"], s) for s in SF70_SCEN]
        d["g_abs"] = ga
        P(f"| {k[0]} | {k[2].split('|')[0]} | {k[1]} | {F(k[3], 2)} | {F(d['tot_t'])} / {F(d['tot_70'])} | " + " | ".join(F(x, 2) for x in ga) +
          f" | {F(math.hypot(d['tot_t'], 0.10))} → {F(math.hypot(d['tot_70'], 0.10))} |")
    P("\n## D6. Opisowo: zgodność SF_rel między runami i SF_rel(428770 all) przeniesiony na 427929 (dane niezależne ⇒ suma ⊕)\n")
    P("| gatunek | pT | ε | SF_rel 427929 | SF_rel 428770 | różnica [σ] | σ(N_c)/N_c 427929 (1 kat.) | ⊕ σ(SF_rel 428770) = σ_x | σ_tot(WP70) 427929 | σ_x/σ_tot(WP70) |")
    P("|---|---|---|---|---|---|---|---|---|---|")
    xr = []
    for (sp, bn, e) in gate:
        a7, a8 = DA.get((sp, bn, MAIN[0], e)), DA.get((sp, bn, MAIN[1], e))
        if not a7 or not a8 or "sf" not in a7 or "sf" not in a8 or a7["tight_rel"] is None:
            continue
        z = (a7["sf"] - a8["sf"]) / math.hypot(a7["sf_tot"], a8["sf_tot"])
        sx = math.hypot(a7["tight_rel"], a8["sf_tot"] / a8["sf"])
        xr.append({"sp": sp, "bin": bn, "eps": e, "z": z, "sx": sx, "ratio": sx / a7["tot_70"]})
        P(f"| {sp} | {bn} | {F(e, 2)} | {pm(a7['sf'], a7['sf_tot'])} | {pm(a8['sf'], a8['sf_tot'])} | {z:+.1f} | {F(a7['tight_rel'])} | {F(sx)} | {F(a7['tot_70'])} | {F(sx / a7['tot_70'], 2)} |")

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, float) and not math.isfinite(o):
            return None
        if isinstance(o, np.bool_):
            return bool(o)
        return o
    json.dump(clean({"gate": {"|".join(map(str, k)): v for k, v in gate.items()}, "mc": {"|".join(map(str, k)): v for k, v in MC.items()},
                     "data": {"|".join(map(str, k)): v for k, v in DA.items()}, "classes": cls, "xrun": xr}), open(a.out, "w"), indent=1)


if __name__ == "__main__":
    sys.exit(main())
