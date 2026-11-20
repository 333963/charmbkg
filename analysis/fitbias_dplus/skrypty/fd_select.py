#!/usr/bin/env python3
"""FITBIAS_Dplus — tabele wyboru (KRYTERIA §4) z out/mc.json i out/toys.json → tabele/mc.md, tabele/select.md, out/select.json. Bez danych okna sygnału."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fd_lib as D  # noqa: E402

O = os.path.join(D.BASE, "out")
R = json.load(open(os.path.join(O, "mc.json")))
TY = json.load(open(os.path.join(O, "toys.json")))["res"]
MC = {(r["bin"], r["model"]): r for r in R}
OOS = {(r["bin"], r["model"], r["comp"]): r for r in json.load(open(os.path.join(O, "mc_oos.json")))} if os.path.exists(os.path.join(O, "mc_oos.json")) else {}
pc = lambda x, e=None: (f"{100 * x:+.1f}" + (f" ± {100 * e:.1f}" if e is not None else ""))  # noqa: E731
S = {}
T = ["# FITBIAS_Dplus — domknięcia na MC (mb_nd, WP70), b = N_fit/N_prawda − 1 [%]", "",
     "σ_b = max(1,3 × HESSE, rms 40 replik Poissona)/N_prawda. Połówki: szablony związane i κ z jednej połowy plików mb_nd, pseudo-dane z drugiej (poza próbką). Asimov: N_prawda × szablon fd + rzeczywiste DMM + CM + PartComb + wygładzona kombinatoryka.", ""]
T += ["## Domknięcie pełne (szablony w próbce)", "", "| ustawienie | " + " | ".join(f"{b}: b | sc | p(χ²) | b (kształt stały)" for b in D.BINS) + " | K1 pełne |", "|---|" + "---|---|---|---|" * 3 + "---|"]


def ok1(b, e):
    return abs(b) <= 0.05 or abs(b) <= e


for m in D.SET:
    row, k1 = [], True
    for bn in D.BINS:
        f = MC[(bn, m)]["full"]
        row.append(f"{pc(f['b'], f['b_err'])} | {f['sc']:.2f} | {f['p']:.2f} | {pc(f['b_fix'])}")
        k1 = bool(k1 and ok1(f["b"], f["b_err"]))
    S[m] = {"K1_full": k1, "b_full": [MC[(bn, m)]["full"]["b"] for bn in D.BINS], "b_full_err": [MC[(bn, m)]["full"]["b_err"] for bn in D.BINS]}
    T.append(f"| {m} | " + " | ".join(row) + f" | {'tak' if k1 else 'NIE'} |")
T += ["", "## Domknięcie poza próbką (połówki)", "", "| ustawienie | " + " | ".join(f"{b}: A→B | B→A | średnia" for b in D.BINS) + " | K1 połówki |", "|---|" + "---|---|---|" * 3 + "---|"]
for m in D.SET:
    row, k1, bh, bhe = [], True, [], []
    for bn in D.BINS:
        a, b = MC[(bn, m)]["half_AtoB"], MC[(bn, m)]["half_BtoA"]
        av, ae = 0.5 * (a["b"] + b["b"]), 0.5 * np.hypot(a["b_err"], b["b_err"])
        row.append(f"{pc(a['b'], a['b_err'])} | {pc(b['b'], b['b_err'])} | {pc(av, ae)}")
        k1 = bool(k1 and ok1(av, ae)); bh.append(float(av)); bhe.append(float(ae))
    S[m].update({"K1_half": k1, "b_half": bh, "b_half_err": bhe})
    T.append(f"| {m} | " + " | ".join(row) + f" | {'tak' if k1 else 'NIE'} |")
T += ["", "## Asimov (składowe powabowe rzeczywiste, kombinatoryka wygładzona): b_A [%]; błąd = rms replik Poissona składowych powabowych", "",
      "| ustawienie | " + " | ".join(f"{b}: w próbce | poza próbką (średnia A→B, B→A)" for b in D.BINS) + " | b̄ (poza próbką) |", "|---|" + "---|---|" * 3 + "---|"]
for m in D.SET:
    row, ba = [], []
    for bn in D.BINS:
        r = MC[(bn, m)]
        i, a, b = r["asimov_in"], r["asimov_AtoB"], r["asimov_BtoA"]
        av = 0.5 * (a["b"] + b["b"])
        ae = 0.5 * np.hypot(a["b_err"] or 0, b["b_err"] or 0)
        row.append(f"{pc(i['b'], i['b_err'])} | {pc(av, ae)}"); ba.append(float(av))
    S[m].update({"bA": ba, "bbar": float(np.mean(np.abs(ba)))})
    T.append(f"| {m} | " + " | ".join(row) + f" | {100 * S[m]['bbar']:.1f} |")
T += ["", "## Rozkład na kategorie prawdy: sygnał pozorny N_spur/N_prawda [%] (kształt sygnału stały — reżim addytywny; w próbce)", "",
      "„model.” = składowa objęta szablonem związanym (w próbce 0 z konstrukcji); w nawiasie (tylko A2) wersja poza próbką z `fd_mc_oos.py` (średnia A→B, B→A: rzeczywiste zdarzenia jednej połowy wobec szablonu z drugiej, na gładkim cokole tła).",
      "Uwaga: liczby `dec.AtoB/BtoA` dla składowych modelowanych w `mc.json` są nieważne (fit bez tła, N_B na granicy 0) — zastąpione przez `mc_oos.json`.", ""]
for bn in D.BINS:
    names = ["DMM", "CM_Ds", "CM_oth", "PC", "NoCharm", "CombB"]
    T += [f"### {bn} GeV", "", "| ustawienie | " + " | ".join(names) + " | suma | b (kształt stały) |", "|---|" + "---|" * (len(names) + 2)]
    for m in D.SET:
        r = MC[(bn, m)]
        row, tot = [], 0.0
        for cn in names:
            d = r["dec"]["in"][cn]
            if d.get("modelled") and "spur" not in d:
                q = OOS.get((bn, m, cn))
                row.append(f"model. ({pc(q['mean'], q['mean_err'])})" if q else "model.")
            else:
                row.append(pc(d["spur"], d["err"])); tot += d["spur"]
        T.append(f"| {m} | " + " | ".join(row) + f" | {pc(tot)} | {pc(r['full']['b_fix'])} |")
    q = OOS.get((bn, "A2", "pk"))
    if q:
        T += ["", f"A2, cały szablon DMM + CM poza próbką (szablon z połowy statystyki): {pc(q['mean'], q['mean_err'])} % (A→B {pc(q['AtoB']['spur'], q['AtoB']['err'])}, B→A {pc(q['BtoA']['spur'], q['BtoA']['err'])})."]
    T.append("")
    pcn = sorted([k for k in MC[(bn, "A0")]["dec"]["in"] if k.startswith("PC|")], key=lambda k: -MC[(bn, "A0")]["dec"]["in"][k]["n"])
    T += ["PartComb per rodzic (A0, kształt stały): " + "; ".join(f"{k[3:]}: {pc(MC[(bn, 'A0')]['dec']['in'][k]['spur'], MC[(bn, 'A0')]['dec']['in'][k]['err'])} (n = {MC[(bn, 'A0')]['dec']['in'][k]['n']:.0f})" for k in pcn), ""]
open(os.path.join(D.BASE, "tabele", "mc.md"), "w").write("\n".join(T) + "\n")

# ---- zabawki
U = ["# FITBIAS_Dplus — zabawki przy czystości danych i wybór (KRYTERIA §3.4, §4)", "",
     "Tło z pasm danych WP70 (okno zaślepione), N_inj = znany uzysk M0, sc_t ~ N(1,05; 0,10), dmu_t ~ N(0; 0,5 MeV); 200 zabawek na komórkę. pull = (N − N_inj)/(1,3 × HESSE).", ""]
cells = [(r["bin"], r["run"]) for r in TY]
U += ["| komórka | N_inj | S/B w oknie | χ²/ndf pasm |", "|---|---|---|---|"] + [f"| {r['run']} {r['bin']} | {r['info']['N_inj']:.0f} | {r['info']['SB_win']:.3f} | {r['info']['chi2_sb']:.0f}/{r['info']['ndf_sb']} |" for r in TY]
U += ["", "| ustawienie | " + " | ".join(f"{run} {bn}: obc. % / σ_toy % / ⟨pull⟩ / rms pull / ważne" for bn, run in cells) + " | K2 komórek | śr. obc. % | mediana σ_toy % | mediana rms(N)/HESSE |", "|---|" + "---|" * (len(cells) + 4)]
for m in D.SET:
    row, nok, bias, st, rat, cellinfo = [], 0, [], [], [], {}
    for r in TY:
        a = np.array(r["rows"][m], float)
        n = r["info"]["N_inj"]
        good = (a[:, 2] > 0) & (a[:, 3] == 0) & np.isfinite(a[:, 1])
        g = a[good]
        b = float(np.mean(g[:, 0]) / n - 1)
        s = float(np.std(g[:, 0]) / n)
        pl = (g[:, 0] - n) / (1.3 * g[:, 1])
        mp, rp, fr = float(np.mean(pl)), float(np.std(pl)), float(good.mean())
        ratio = float(np.std(g[:, 0]) / np.median(g[:, 1]))
        okc = abs(mp) <= 0.25 and rp <= 1.25 and fr >= 0.90
        nok += int(okc); bias.append(b); st.append(s); rat.append(ratio)
        cellinfo[f"{r['run']}|{r['bin']}"] = {"bias": b, "sig_toy": s, "pull_mean": mp, "pull_rms": rp, "valid": fr, "rms_over_hesse": ratio, "ok": bool(okc), "pull_mean_ok": bool(abs(mp) <= 0.25)}
        row.append(f"{100 * b:+.1f} / {100 * s:.1f} / {mp:+.2f} / {rp:.2f} / {fr:.2f}")
    mb = float(np.mean(bias))
    S[m].update({"toys": cellinfo, "K2_cells": int(nok), "toy_bias": mb, "sig_toy": float(np.median(st)), "K2": bool(nok >= 5 and abs(mb) <= 0.02),
                 "K2_relaxed": bool(sum(c["pull_mean_ok"] and c["valid"] >= 0.9 for c in cellinfo.values()) >= 5 and abs(mb) <= 0.02)})
    S[m]["Q"] = float(np.hypot(S[m]["bbar"], S[m]["sig_toy"]))
    S[m]["K1"] = bool(S[m]["K1_full"] and S[m]["K1_half"])
    U.append(f"| {m} | " + " | ".join(row) + f" | {nok}/6 | {100 * mb:+.2f} | {100 * S[m]['sig_toy']:.1f} | {np.median(rat):.2f} |")
U += ["", "## Zestawienie", "", "| ustawienie | K1 pełne | K1 połówki | K2 | K2 złagodzone (bez rms) | b̄ Asimov poza próbką % | σ_toy % | Q |", "|---|---|---|---|---|---|---|---|"]
t = lambda x: "tak" if x else "NIE"  # noqa: E731
for m in D.SET:
    s = S[m]
    U.append(f"| {m} | {t(s['K1_full'])} | {t(s['K1_half'])} | {t(s['K2'])} ({s['K2_cells']}/6) | {t(s['K2_relaxed'])} | {100 * s['bbar']:.1f} | {100 * s['sig_toy']:.1f} | {s['Q']:.4f} |")
open(os.path.join(D.BASE, "tabele", "select.md"), "w").write("\n".join(U) + "\n")
json.dump(S, open(os.path.join(O, "select.json"), "w"), indent=1)
print("\n".join(U[-12:]))
