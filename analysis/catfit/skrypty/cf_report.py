#!/usr/bin/env python3
"""CATFIT — odczyt: out/{mc,toys,data}.json + XSEC_final/out/final.json (tylko odczyt) -> out/catfit.json, tabele/catfit.md. Wlasna kopia rachunku przekroju (KRYTERIA §4). POST HOC."""
import json
import math
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
BASE = os.path.join(HERE, "..")
CW = os.path.join(BASE, "..")
J = lambda *p: json.load(open(os.path.join(*p)))  # noqa: E731
FINALL = J(CW, "XSEC_final_2026-11-20", "out", "final.json")
FIN = FINALL["cells"]
LUM, LREL = FINALL["L_nb"], FINALL["L_rel"]
HSC = {"D0": 1.3, "Dplus": 1.0}
NM = ["L", "C", "B", "A"]
BINS = ["3-5", "5-8", "8-30"]
COMMON = ["BR", "refl", "nonprompt", "evsel", "tracking", "ip_material", "vertex3", "generator"]
DSTAR = {"5-8": {"0.7": (1.005, 0.126), "0.5": (0.953, 0.130), "0.8": (1.080, 0.128)}, "8-30": {"0.7": (0.965, 0.097), "0.5": (0.967, 0.119), "0.8": (0.935, 0.084)},
         "3-5": {"0.7": (1.309, 0.426), "0.5": (1.282, 0.434), "0.8": (1.399, 0.393)}}        # dnet4_results.json, drzewo H, 427929


def q(*a):
    return math.sqrt(sum(x * x for x in a))


def pc(x, n=1):
    return f"{100 * x:.{n}f}"


def cov(pull):
    pull = pull[np.isfinite(pull)]
    return float((np.abs(pull) < 1).mean()) if len(pull) else float("nan")


def toys_summary(T):
    out = {}
    for key in sorted({(r["sp"], r["bin"], r["model"], r["dord"]) for r in T}):
        sp, bn, mn, dord = key
        rs = [r for r in T if (r["sp"], r["bin"], r["model"], r["dord"]) == key]
        R = np.array([x for r in rs for x in r["rows"]])
        n70, nt, sfr = rs[0]["N70_inj"], rs[0]["Ntot_inj"], rs[0]["SF_rel"]
        o = {"n": len(R), "SF_rel": sfr, "info": rs[0]["info"]}
        for nm, i, tr in (("WP70", 0, n70), ("F0", 3, nt), ("F70", 10, nt), ("Ffix", 13, nt)):
            ok = R[:, i + 2] > 0.5
            x, e = R[ok, i] / tr, R[ok, i + 1] / tr * HSC[sp]
            o[nm] = {"valid": float(ok.mean()), "bias": float(x.mean() - 1), "bias_err": float(x.std(ddof=1) / math.sqrt(max(len(x), 1))), "rms": float(x.std(ddof=1)), "med_err": float(np.nanmedian(e)),
                     "cover": cov((x - 1) / e), "rms68": float(0.5 * (np.percentile(x, 84) - np.percentile(x, 16)))}
        o["tot_cut"] = q(o["WP70"]["rms"], sfr)                       # sciezka nominalna: N70/(eps SF), SF niezalezny
        o["tot_F0"] = o["F0"]["rms"]
        o["tot_F70"] = o["F70"]["rms"]                                # srodek ograniczenia losowany z sigma_SF
        o["tot_Ffix_plus_SF"] = q(o["Ffix"]["rms"], sfr)
        ok = R[:, 5] > 0.5
        o["f_mean"] = dict(zip(NM, R[ok, 6:10].mean(0).tolist())); o["f_rms"] = dict(zip(NM, R[ok, 6:10].std(0, ddof=1).tolist()))
        out["|".join(map(str, key))] = o
    return out


def main():
    MC = {(r["sp"], r["bin"], r["model"]): r for r in J(BASE, "out", "mc.json")["res"]}
    DA = {(r["sp"], r["bin"], r["model"]): r for r in J(BASE, "out", "data.json")["res"]}
    TS = toys_summary(J(BASE, "out", "toys.json")["res"])
    R = {"_meta": {"zadanie": "CATFIT_2026-11-20", "status": "POST HOC wzgledem XSEC_final; PRELIMINARY", "kategorie": "L=[100,70) C=[70,50) B=[50,30) A=[30,0) % wydajnosci sygnalu MC (drzewo H)",
                   "run": "427929|MB427929", "L_nb": LUM}, "cells": {}, "toys": TS}
    for (sp, bn, mn), d in DA.items():
        fin, mc = FIN[f"{sp}|{bn}"], MC[(sp, bn, mn)]
        nom = (mn == {"D0": "M8", "Dplus": "M0"}[sp])
        k = mc["K4"]
        Tt = k["truth_tot"]
        assert abs(mc["WP"]["0.7"]["truth"] - fin["N_reco_truth"]) < 0.5, (sp, bn, mc["WP"]["0.7"]["truth"], fin["N_reco_truth"])
        Ct = Tt / fin["N_gen_excl"]
        common = {r: fin["rel"][r] for r in COMMON}
        c = {"model": mn, "nominal_as_final": nom, "C_tot": Ct, "C_WP70": fin["C"], "N_reco_truth_tot": Tt}
        # ---- obecny (final.json; tylko dla modelu nominalnego jest to scisle ten sam fit)
        c["final"] = {"N": fin["N"], "rel": fin["rel"], "syst_rel": fin["syst_rel"], "tot_rel": fin["tot_rel"], "sigma_ub": fin["sigma_ub"], "stat_ub": fin["stat_ub"], "syst_ub": fin["syst_ub"]}
        # ---- domkniecie MC
        b = k["Ntot"] / Tt - 1
        sb = max(1.3 * (k["Ntot_err"] or 0), k["rep_rms_tot"]) / Tt
        w7 = mc["WP"]["0.7"]
        c["mc"] = {"b_tot": b, "sb_tot": sb, "b_70": w7["N"] / w7["truth"] - 1, "sb_70": max(1.3 * (w7["err"] or 0), w7["rep_rms"]) / w7["truth"], "p": k["p"], "p_cat": k["p_cat"], "valid": k["valid"], "sc": k["sc"],
                   "f_fit": k["f"], "f_true": {n: k["truth"][n] / Tt for n in NM}, "f_err": {n: max((k.get("f_err") or {}).get(n, 0), k["rep_rms_f"][n]) for n in NM},
                   "E70_fit": 1 - k["f"]["L"], "E70_true": 1 - k["truth"]["L"] / Tt}
        # ---- F0
        o = d["F0"]
        N = o["Ntot"]
        rel = {"stat": max(HSC[sp] * o["Ntot_err"], o["boot_tot"]) / N,
               "fit_model": max(abs(v["Ntot"] / N - 1) for v in o["var"].values() if v["valid"]),
               "fit_bias": q(b, sb)}
        if sp == "D0":
            rel["peak_tmpl"] = fin["rel"]["peak_tmpl"]
            rel["peak_norm"] = 0.5 * abs(o["k+5.8"] - o["k-5.8"]) / N
            rel["fit_altmin"] = abs(o["altmin"]["Ntot"] / N - 1) if o.get("altmin") else 0.0
        rel["SF"] = 0.0
        rel["mc_stat"] = math.sqrt((1 - Ct) / Tt)
        rel.update(common)
        syst = q(*[v for kk, v in rel.items() if kk != "stat"])
        sig = N / (Ct * fin["BR"] * LUM) / 1e3
        fd = d["fd_e"]
        fk = {n: {"data": o["f"][n], "err": max(o["f_err"][n], o["boot_f"][n]), "mc_fd": fd[n], "SF_k": o["f"][n] / fd[n], "SF_k_err": max(o["f_err"][n], o["boot_f"][n]) / fd[n],
                  "mcfit_over_truth": c["mc"]["f_fit"][n] / c["mc"]["f_true"][n], "mcfit_err": c["mc"]["f_err"][n] / c["mc"]["f_true"][n]} for n in NM}
        E70, E70e = o["E"]["C"], max(o["E_err"]["C"], o["boot_E70"])
        e70 = 1 - fd["L"]
        c["F0"] = {"Ntot": N, "hesse": o["Ntot_err"], "boot": o["boot_tot"], "N_k": o["N"], "err_k": o["err"], "valid": o["valid"], "at_limit": o["at_limit"], "p": o["p"], "p_cat": o["p_cat"], "sc": o["sc"], "dmu": o["dmu"],
                   "rel": rel, "syst_rel": syst, "tot_rel": q(rel["stat"], syst, LREL), "sigma_ub": sig, "stat_ub": sig * rel["stat"], "syst_ub": sig * syst, "lumi_ub": sig * LREL,
                   "f_k": fk, "E70": E70, "E70_err": E70e, "SF70_fit": E70 / e70, "SF70_fit_err": E70e / e70, "SF70_mcfit": c["mc"]["E70_fit"] / c["mc"]["E70_true"],
                   "Nsub": o["Nsub"], "Nsub_err": o["Nsub_err"], "var": {kk: v["Ntot"] for kk, v in o["var"].items()}, "var_valid": {kk: v["valid"] for kk, v in o["var"].items()},
                   "Lplus2": {"Ntot": o["L+2"]["Ntot"], "err": o["L+2"]["Ntot_err"], "p": o["L+2"]["p"], "p_L": o["L+2"]["p_cat"]["L"], "valid": o["L+2"]["valid"]}, "altmin": o.get("altmin")}
        c["F0"]["g"] = c["F0"]["tot_rel"] / fin["tot_rel"]
        # ---- F70 / F3
        for wh in ("F70", "F3"):
            o7 = d[wh]
            N7 = o7["Ntot"]
            hf = o7["Ntot_err_fixfrac"] or o7["Ntot_err"]
            r7 = {"stat": max(HSC[sp] * hf, o7["boot_tot"]) / N7, "SF": math.sqrt(max(o7["Ntot_err"] ** 2 - hf ** 2, 0)) / N7, "fit_model": max(abs(v / N7 - 1) for v in o7["var"].values()),
                  "fit_bias_lo": min(fin["rel"]["fit_bias"], rel["fit_bias"]), "fit_bias_hi": max(fin["rel"]["fit_bias"], rel["fit_bias"])}
            extra = [rel.get("peak_tmpl", 0), rel.get("peak_norm", 0), rel.get("fit_altmin", 0), rel["mc_stat"]] + list(common.values())
            s_lo, s_hi = q(r7["SF"], r7["fit_model"], r7["fit_bias_lo"], *extra), q(r7["SF"], r7["fit_model"], r7["fit_bias_hi"], *extra)
            s7 = N7 / (Ct * fin["BR"] * LUM) / 1e3
            c[wh] = {"Ntot": N7, "hesse": o7["Ntot_err"], "hesse_fixfrac": hf, "boot": o7["boot_tot"], "valid": o7["valid"], "p": o7["p"], "f": o7["f"], "rel": r7, "tot_rel_lo": q(r7["stat"], s_lo, LREL), "tot_rel_hi": q(r7["stat"], s_hi, LREL),
                     "sigma_ub": s7, "g_lo": q(r7["stat"], s_lo, LREL) / fin["tot_rel"], "g_hi": q(r7["stat"], s_hi, LREL) / fin["tot_rel"]}
        c["Ffix"] = dict(d["Ffix"], rel_hesse=d["Ffix"]["Ntot_err"] / d["Ffix"]["Ntot"] if d["Ffix"]["Ntot_err"] else None)
        c["K3"] = {"Ntot": d["K3"]["Ntot"], "err": d["K3"]["Ntot_err"], "p": d["K3"]["p"]}
        c["K5"] = {"Ntot": d["K5"]["Ntot"], "err": d["K5"]["Ntot_err"], "p": d["K5"]["p"], "p_cat": d["K5"]["p_cat"], "valid": d["K5"]["valid"]}
        w70 = d["WP"]["0.7"]
        c["WP70_refit"] = {"N": w70["N"], "err": w70["err"], "boot": w70["boot"], "p": w70["p"]}
        # rozklad zysku (sama czesc fit: sigma(N)/N, HESSE x skala; sciezka nominalna = N70 (+) SF)
        sfr = fin["SF_err"] / fin["SF"]
        c["decomp"] = {"cut_stat": HSC[sp] * w70["err"] / w70["N"], "cut_stat_plus_SF": q(HSC[sp] * w70["err"] / w70["N"], sfr), "F0": HSC[sp] * d["F0"]["Ntot_err"] / N,
                       "F70": q(HSC[sp] * c["F70"]["hesse_fixfrac"], c["F70"]["rel"]["SF"] * c["F70"]["Ntot"]) / c["F70"]["Ntot"], "F3": q(HSC[sp] * c["F3"]["hesse_fixfrac"], c["F3"]["rel"]["SF"] * c["F3"]["Ntot"]) / c["F3"]["Ntot"], "Ffix": HSC[sp] * d["Ffix"]["Ntot_err"] / d["Ffix"]["Ntot"], "SF_rel": sfr}
        # ---- werdykt F0 / F70
        t = TS.get(f"{sp}|{bn}|{mn}|0", {})
        k1 = c["F0"]["g"] <= 0.85
        k2 = abs(b) <= 0.05 or abs(b) <= sb
        k3 = bool(t) and 0.60 <= t["F0"]["cover"] <= 0.76 and t["F0"]["valid"] >= 0.90
        k4 = o["valid"] and not o["at_limit"] and o["p"] >= 1e-3 and min(o["p_cat"].values()) >= 1e-3
        c["verdict_F0"] = {"g": c["F0"]["g"], "K1_g": k1, "K2_bias": k2, "K3_toys": k3, "K4_gof": bool(k4), "verdict": "ZYSK" if (k1 and k2 and k3 and k4) else ("maly zysk" if (c["F0"]["g"] <= 0.95 and k2 and k3 and k4) else "BRAK ZYSKU")}
        k3b = bool(t) and 0.60 <= t["F70"]["cover"] <= 0.76 and t["F70"]["valid"] >= 0.90
        g7 = c["F70"]["g_lo"]
        c["verdict_F70"] = {"g_lo": g7, "g_hi": c["F70"]["g_hi"], "K3_toys": k3b, "verdict": "ZYSK" if (g7 <= 0.85 and k3b and d["F70"]["valid"] and d["F70"]["p"] >= 1e-3) else ("maly zysk" if (g7 <= 0.95 and k3b) else "BRAK ZYSKU")}
        # ---- B2
        b2 = {}
        base = [fin["rel"].get("peak_tmpl", 0), fin["rel"].get("peak_norm", 0), fin["rel"].get("fit_altmin", 0), fin["rel"]["mc_stat"]] + list(common.values())
        for wp, w in d["WP"].items():
            m = mc["WP"][wp]
            bb, sbb = m["N"] / m["truth"] - 1, max(1.3 * (m["err"] or 0), m["rep_rms"]) / m["truth"]
            st = max(HSC[sp] * w["err"], w["boot"]) / w["N"]
            fm = max(abs(v / w["N"] - 1) for v in w["var"].values())
            rest = q(st, fm, q(bb, sbb), LREL, *base)
            x = (0.85 * fin["tot_rel"]) ** 2 - rest ** 2
            e = {"N": w["N"], "N_over_e": w["N"] / w["e"], "stat": st, "fit_model": fm, "b": bb, "sb": sbb, "rest": rest, "SF_max_for_gain": math.sqrt(x) if x > 0 else None, "p": w["p"], "valid": w["valid"], "at_limit": w["at_limit"]}
            if sp == "D0" and wp in DSTAR[bn]:
                sf, se = DSTAR[bn][wp]
                e["SF_Dstar"], e["SF_Dstar_err"] = sf, se
                e["tot_with_Dstar"] = q(rest, se / sf); e["g_with_Dstar"] = e["tot_with_Dstar"] / fin["tot_rel"]
            b2[wp] = e
        r70 = b2["0.7"]["rest"]
        for wp, e in b2.items():
            e["rest_over_rest70"] = e["rest"] / r70
        c["B2"] = b2
        R["cells"][f"{sp}|{bn}|{mn}"] = c
    json.dump(R, open(os.path.join(BASE, "out", "catfit.json"), "w"), indent=1, default=float)

    # ------------------------------------------------------------------ tabele
    L = ["# CATFIT — tabele (generowane przez `skrypty/cf_report.py`; POST HOC względem XSEC_final; run 427929, piki inkluzywne, bez AFP)", ""]
    keys = list(R["cells"])
    L += ["## T1. Błąd całkowity: obecny (`final.json`, WP70) → fit K = 4 ze swobodnymi ułamkami (F0) → z ograniczeniem E₇₀ (F70)", "",
          "| gat. | bin | model | obecny: stat / syst / razem [%] | F0: N_tot | F0: stat / syst / razem [%] | g(F0) | F70: stat / ułamki (więz SF + swoboda C,B,A) / razem [%] (obc. min–max) | g(F70) | werdykt F0 / F70 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for kk in keys:
        c = R["cells"][kk]; sp, bn, mn = kk.split("|"); f, o, o7 = c["final"], c["F0"], c["F70"]
        L.append(f"| {sp} | {bn} | {mn}{'' if c['nominal_as_final'] else ' (wariant)'} | {pc(f['rel']['stat'])} / {pc(f['syst_rel'])} / {pc(f['tot_rel'])} | {o['Ntot']:.0f} ± {max(HSC[sp] * o['hesse'], o['boot']):.0f} | "
                 f"{pc(o['rel']['stat'])} / {pc(o['syst_rel'])} / {pc(o['tot_rel'])} | {o['g']:.2f} | {pc(o7['rel']['stat'])} / {pc(o7['rel']['SF'])} / {pc(o7['tot_rel_lo'])}–{pc(o7['tot_rel_hi'])} | {o7['g_lo']:.2f}–{o7['g_hi']:.2f} | "
                 f"{c['verdict_F0']['verdict']} / {c['verdict_F70']['verdict']} |")
    L += ["", "## T2. Budżet F0 [%] (wiersze wspólne przepisane z `final.json`)", "", "| komórka | stat | model fitu | obciążenie (|b| ⊕ σ_b) | κ ± 5,8 % | drugie min. | stat MC | wspólne | razem syst | p(χ²) łączne | p per kat. (L / C / B / A) | sc |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for kk in keys:
        o = R["cells"][kk]["F0"]; r = o["rel"]
        L.append(f"| {kk} | {pc(r['stat'])} | {pc(r['fit_model'])} | {pc(r['fit_bias'])} | {pc(r.get('peak_norm', 0), 2)} | {pc(r.get('fit_altmin', 0))} | {pc(r['mc_stat'])} | {pc(q(*[r[x] for x in COMMON]))} | {pc(o['syst_rel'])} | {o['p']:.3g} | "
                 + " / ".join(f"{o['p_cat'][n]:.2g}" for n in NM) + f" | {o['sc']:.3f} |")
    L += ["", "## T3. Domknięcie na mb_nd (pełne tło): N_tot (K = 4, F0) i N₇₀ (jedna kategoria)", "", "| komórka | b(N_tot) ± σ_b [%] | b(N₇₀) ± σ_b [%] | f_k fit / prawda: L | C | B | A | E₇₀ fit / prawda | p(χ²) |", "|---|---|---|---|---|---|---|---|---|"]
    for kk in keys:
        m = R["cells"][kk]["mc"]
        L.append(f"| {kk} | {100 * m['b_tot']:+.1f} ± {pc(m['sb_tot'])} | {100 * m['b_70']:+.1f} ± {pc(m['sb_70'])} | " + " | ".join(f"{m['f_fit'][n]:.3f} ± {m['f_err'][n]:.3f} / {m['f_true'][n]:.3f}" for n in NM)
                 + f" | {m['E70_fit']:.3f} / {m['E70_true']:.3f} | {m['p']:.2g} |")
    L += ["", "## T4. Zabawki przy czystości danych (tło z pasm danych, rząd nominalny; N wstrzyknięte = N₇₀(final)/0,70 × f_k^MC)", "",
          "| komórka | n | WP70: obc. [%] / rms [%] / pokrycie | F0: obc. / rms / pokrycie / ważne | F70: obc. / rms / pokrycie | Ffix: obc. / rms / pokrycie | błąd ścieżki: cięcie ⊕ SF → F0 → F70 → Ffix ⊕ SF [%] |", "|---|---|---|---|---|---|---|"]
    for kk in keys:
        t = TS.get(kk + "|0")
        if not t:
            continue
        s = lambda x: f"{100 * x['bias']:+.1f} / {pc(x['rms'])} / {x['cover']:.2f}"  # noqa: E731
        L.append(f"| {kk} | {t['n']} | {s(t['WP70'])} | {s(t['F0'])} / {t['F0']['valid']:.2f} | {s(t['F70'])} | {s(t['Ffix'])} | {pc(t['tot_cut'])} → {pc(t['tot_F0'])} → {pc(t['tot_F70'])} → {pc(t['tot_Ffix_plus_SF'])} |")
    L += ["", "Opisowo — prawda tła rzędu + 1, fit rzędem nominalnym (60 zabawek): obciążenie [%] WP70 / F0 / F70 / Ffix:", ""]
    for kk in keys:
        t = TS.get(kk + "|1")
        if t:
            L.append(f"- {kk}: {100 * t['WP70']['bias']:+.1f} ± {pc(t['WP70']['bias_err'])} / {100 * t['F0']['bias']:+.1f} ± {pc(t['F0']['bias_err'])} / {100 * t['F70']['bias']:+.1f} ± {pc(t['F70']['bias_err'])} / {100 * t['Ffix']['bias']:+.1f} ± {pc(t['Ffix']['bias_err'])}; pokrycie F0 {t['F0']['cover']:.2f}")
    L += ["", "## T5. Dane: ułamki f_k i SF_k = f_k^dane/f_k^MC(fd_nd) z fitu F0; obok f_k^fit/f_k^prawda z mb_nd", "",
          "| komórka | kat. | N_k ± σ | f_k dane | f_k MC | SF_k (dane/MC) | MC: fit/prawda |", "|---|---|---|---|---|---|---|"]
    for kk in keys:
        o = R["cells"][kk]["F0"]
        for n in NM:
            x = o["f_k"][n]
            L.append(f"| {kk} | {n} | {o['N_k'][n]:.0f} ± {o['err_k'][n]:.0f} | {x['data']:.3f} ± {x['err']:.3f} | {x['mc_fd']:.3f} | {x['SF_k']:.2f} ± {x['SF_k_err']:.2f} | {x['mcfit_over_truth']:.2f} ± {x['mcfit_err']:.2f} |")
    L += ["", "| komórka | SF₇₀ z fitu K = 4 (E₇₀^dane/E₇₀^MC) | to samo na mb_nd (fit/prawda) | SF₇₀ w `final.json` | tag D* 427929 WP70 / WP50 (dnet4) |", "|---|---|---|---|---|"]
    for kk in keys:
        c = R["cells"][kk]; sp, bn, mn = kk.split("|"); o = c["F0"]; fin = FIN[f"{sp}|{bn}"]
        ds = DSTAR[bn]
        L.append(f"| {kk} | {o['SF70_fit']:.2f} ± {o['SF70_fit_err']:.2f} | {o['SF70_mcfit']:.2f} | {fin['SF']:.3f} ± {fin['SF_err']:.3f} | " + (f"{ds['0.7'][0]:.2f} ± {ds['0.7'][1]:.2f} / {ds['0.5'][0]:.2f} ± {ds['0.5'][1]:.2f}" if sp == "D0" else "—") + " |")
    L += ["", "## T6. Rozkład zysku: σ(N)/N z fitu [%] (część stat. = HESSE przy ułamkach zamrożonych × skala reguły stat.; część ułamków/SF bez skali)", "", "| komórka | σ_SF [%] | cięcie WP70 | cięcie ⊕ SF | F0 (swobodne) | F70 (1 więz) | F3 (3 więzy — założenie toya) | Ffix (ułamki z MC, bez błędu SF) |", "|---|---|---|---|---|---|---|---|"]
    for kk in keys:
        x = R["cells"][kk]["decomp"]
        L.append(f"| {kk} | {pc(x['SF_rel'])} | {pc(x['cut_stat'])} | {pc(x['cut_stat_plus_SF'])} | {pc(x['F0'])} | {pc(x['F70'])} | {pc(x['F3'])} | {pc(x['Ffix'])} |")
    L += ["", "## T7. Przekrój [µb]: `final.json` vs CATFIT", "", "| komórka | final: σ ± stat ± syst | F0: σ ± stat ± syst | F0/final | F70: σ (razem %) | F70/final | C_tot [%] (C_WP70) | kontrole: N₇₀ refit / K = 3 / ΣCBA z K = 4 | K = 5 N_tot | F0 z rzędem + 2 w L: N_tot (p_L) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for kk in keys:
        c = R["cells"][kk]; f, o, o7 = c["final"], c["F0"], c["F70"]
        L.append(f"| {kk} | {f['sigma_ub']:.1f} ± {f['stat_ub']:.1f} ± {f['syst_ub']:.1f} | {o['sigma_ub']:.1f} ± {o['stat_ub']:.1f} ± {o['syst_ub']:.1f} | {o['sigma_ub'] / f['sigma_ub']:.2f} | {o7['sigma_ub']:.1f} ({pc(o7['tot_rel_lo'])}–{pc(o7['tot_rel_hi'])}) | {o7['sigma_ub'] / f['sigma_ub']:.2f} | "
                 f"{pc(c['C_tot'], 2)} ({pc(c['C_WP70'], 2)}) | {c['WP70_refit']['N']:.0f} / {c['K3']['Ntot']:.0f} ± {c['K3']['err'] or 0:.0f} / {o['Nsub']:.0f} ± {o['Nsub_err']:.0f} | {c['K5']['Ntot']:.0f} ± {c['K5']['err'] or 0:.0f} | {o['Lplus2']['Ntot']:.0f} ± {o['Lplus2']['err'] or 0:.0f} ({o['Lplus2']['p_L']:.2g}) |")
    L += ["", "## T8. B2 — pojedynczy punkt pracy: σ_reszta (bez SF) i dopuszczalny błąd SF dla zysku ≥ 15 %", "", "| komórka | WP | N | N/ε względem WP70 | stat [%] | model [%] | b ± σ_b (mb_nd) [%] | σ_reszta [%] | reszta/reszta(70) | σ_SF^max [%] | z SF tagu D*: razem [%] (g) |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for kk in keys:
        c = R["cells"][kk]; n70 = c["B2"]["0.7"]["N_over_e"]
        for wp, e in c["B2"].items():
            L.append(f"| {kk} | {wp} | {e['N']:.0f} | {e['N_over_e'] / n70:.2f} | {pc(e['stat'])} | {pc(e['fit_model'])} | {100 * e['b']:+.1f} ± {pc(e['sb'])} | {pc(e['rest'])} | {e['rest_over_rest70']:.2f} | "
                     + (pc(e["SF_max_for_gain"]) if e["SF_max_for_gain"] else "niemożliwe") + " | " + (f"{pc(e['tot_with_Dstar'])} ({e['g_with_Dstar']:.2f})" if "g_with_Dstar" in e else "—") + " |")
    open(os.path.join(BASE, "tabele", "catfit.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
