#!/usr/bin/env python3
"""FITBIAS — wybór wariantu wg KRYTERIA §2.3 (tylko MC + zabawki). Wejście: out/var_{asimov,mc,toys}.json → tabele/select.md, out/select.json"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
O = lambda n: os.path.join(HERE, "..", "out", n)  # noqa: E731
MS = [f"M{i}" for i in range(10)]
SIMPLE = ["M0", "M2", "M1", "M6", "M7", "M3", "M5", "M4", "M9", "M8"]
DESC = {"M0": "nominał (wolny)", "M1": "(a) Gauss", "M2": "(a′) stały", "M3": "(b) wolny + pik", "M4": "(a+b) Gauss + pik", "M5": "(a′+b) stały + pik", "M6": "(c1) exp(nom+1)",
        "M7": "(c2) zakres alt.", "M8": "(a+b+c2)", "M9": "(a+b+c1)"}
BINS = ["3-5", "5-8", "8-30"]
f = lambda x, n=3: "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:+.{n}f}".replace(".", ",")  # noqa: E731
u = lambda x, n=3: "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{n}f}".replace(".", ",")  # noqa: E731


def wmean(v, e):
    v, e = np.array(v, float), np.array(e, float)
    w = 1 / e ** 2
    return float((v * w).sum() / w.sum()), float(1 / math.sqrt(w.sum()))


A = json.load(open(O("var_asimov.json")))["res"]
MC = json.load(open(O("var_mc.json")))["res"]
T = json.load(open(O("var_toys.json")))["res"]
L = ["# FITBIAS — wybór wariantu modelu (KRYTERIA §2.3): tylko MC i zabawki, bez okna sygnału danych", ""]
sel = {}
for sp in ("D0", "Dplus"):
    L += [f"## {sp}", ""]
    # ---- S1 Asimov
    L += ["### S1. Asimov T (sygnał + odbicie + DMM + CM prawdy mb_nd + wygładzone PartComb + Comb): b_A = N_fit/N_prawda − 1", "",
          "| ustawienie | " + " | ".join(f"{b} {c}" for c in ("WP80", "WP50", "L80-50") for b in BINS) + " | ⟨|b_A|⟩ WP80 | f/f_prawda − 1 (J; 3 biny) |", "|---|" + "---|" * 11]
    bA = {}
    for m in MS:
        rows = [r for r in A if r["sp"] == sp]
        vals = {c: [r["S"][m][c]["N"] / r["truth"][c] - 1 for r in rows] for c in ("WP80", "WP50", "L80-50")}
        bA[m] = float(np.mean(np.abs(vals["WP80"])))
        fj = [r["J"][m]["f"] / (r["truth"]["WP50"] / (r["truth"]["WP50"] + r["truth"]["L80-50"])) - 1 for r in rows]
        L.append(f"| {m} {DESC[m]} | " + " | ".join(f(v) for c in ("WP80", "WP50", "L80-50") for v in vals[c]) + f" | {u(bA[m])} | " + " / ".join(f(v) for v in fj) + " |")
    # ---- S2 MC closure W
    L += ["", "### S2. Domknięcie mb_nd (W, pełne tło, jedna realizacja): b ± σ_HESSE/N", "",
          "| ustawienie | WP70: 3–5 | 5–8 | 8–30 | średnia ważona WP70 | D2 | B50: 3–5 | 5–8 | 8–30 | średnia B50 | WP90 średnia | L90-70 średnia | A50 średnia |", "|---|" + "---|" * 12]
    D2 = {}
    mcs = [r for r in MC if r["sp"] == sp]
    for m in MS:
        cells = {}
        for cn in ("WP70", "B50", "WP90", "L90-70", "WP50"):
            b = [r["S"][cn][m]["N"] / r["S"][cn]["N_true"] - 1 for r in mcs]
            e = [(r["S"][cn][m]["err"] or float("nan")) / r["S"][cn]["N_true"] for r in mcs]
            cells[cn] = (b, e, wmean(b, e) if all(np.isfinite(e)) else (float("nan"), float("nan")))
        b, e, (mb_, me_) = cells["WP70"]
        D2[m] = bool(abs(mb_) <= 2 * me_ and all(abs(x) <= 2.5 * y for x, y in zip(b, e)))
        L.append(f"| {m} | " + " | ".join(f"{f(x)} ± {u(y)}" for x, y in zip(b, e)) + f" | {f(mb_)} ± {u(me_)} | {'TAK' if D2[m] else 'NIE'} | "
                 + " | ".join(f"{f(x)} ± {u(y)}" for x, y in zip(*cells["B50"][:2])) + f" | {f(cells['B50'][2][0])} ± {u(cells['B50'][2][1])} | {f(cells['WP90'][2][0])} ± {u(cells['WP90'][2][1])} | "
                 f"{f(cells['L90-70'][2][0])} ± {u(cells['L90-70'][2][1])} | {f(cells['WP50'][2][0])} ± {u(cells['WP50'][2][1])} |")
    # ---- S3 pass/fail MC
    L += ["", "### S3. Pass/fail na mb_nd (tryb J): SF_rel^MC = f_fit/f_MC ± σ (σ = max(HESSE, rms 30 replik)); prawda = f_prawda/f_MC; E/N_prawda,70", "",
          "| ustawienie | ε_t | 3–5: SF ± σ (prawda) | 5–8 | 8–30 | średnia SF/SF_prawda − 1 | bramka σ ≤ 0,05 ∧ |Δ| ≤ 2,5σ (biny) | E/N_prawda − 1: 3–5 / 5–8 / 8–30 | N_70/N_prawda − 1 |", "|---|" + "---|" * 8]
    gate = {}
    for m in MS:
        for et in ("0.5", "0.4", "0.35"):
            sf, se, st, ok, Eb, Nb = [], [], [], [], [], []
            for r in mcs:
                o = r["J"][et]; x = o[m]
                s = max(x.get("f_err") or 0, x.get("f_rms") or 0) / o["f_mc"]
                sf.append(x["f"] / o["f_mc"]); se.append(s); st.append(o["f_true"] / o["f_mc"])
                ok.append(bool(s <= 0.05 and abs(sf[-1] - st[-1]) <= 2.5 * s and x["valid"]))
                nt = o["NA_true"] + o["NB_true"]
                Eb.append(x["E"] / nt - 1); Nb.append(x["N70"] / nt - 1)
            gate[(m, et)] = ok
            mm, me = wmean([a / b - 1 for a, b in zip(sf, st)], [a / b for a, b in zip(se, st)])
            L.append(f"| {m} | {et.replace('.', ',')} | " + " | ".join(f"{u(a)} ± {u(b)} ({u(c)})" for a, b, c in zip(sf, se, st)) + f" | {f(mm)} ± {u(me)} | " + " / ".join("TAK" if k else "NIE" for k in ok)
                     + " | " + " / ".join(f(v) for v in Eb) + " | " + " / ".join(f(v) for v in Nb) + " |")
    # ---- S4 toys
    L += ["", "### S4. Zabawki: sygnał MC (kształt losowany z priorytetu) wstrzyknięty w tło z pasm danych; 12 komórek (3 biny × 2 runy × 2 poziomy S/B)", "",
          "| ustawienie | D1: komórek OK / 12 | ⟨pull⟩ N_70 (zakres) | rms pull (zakres) | pokrycie 68 % (mediana) | obciążenie wzgl. N_70 (średnia; zakres) | σ_toy = mediana rms(N/N_wstrz.) | mediana σ_HESSE/N | udział OK (min) | "
          "f: ⟨pull⟩ / rms pull (mediany) | mediana rms(SF_rel) | (d): komórek z rms(E) ≤ 0,95·rms(N_70) | D1 dla E |", "|---|" + "---|" * 12]
    toys = [r for r in T if r["sp"] == sp]
    stat, D1, dJ = {}, {}, {}
    per_cell = {}
    for m in MS:
        nok, pm, pr, cov, rb, rms, hes, okf, fpm, fpr, frms, dcount, e_ok = 0, [], [], [], [], [], [], [], [], [], [], 0, 0
        for r in toys:
            a = np.array(r["rows"][m], float)
            Ni = sum(r["N_inj"]); fi = r["N_inj"][0] / Ni
            valid = a[:, 2] > 0
            ok = valid & (a[:, 3] == 0) & np.isfinite(a[:, 1]) & (a[:, 1] > 0)
            pull = (a[ok, 0] - Ni) / a[ok, 1]
            c_ok = abs(pull.mean()) <= 0.25 and pull.std() <= 1.25 and ok.mean() >= 0.90
            nok += c_ok
            pm.append(pull.mean()); pr.append(pull.std()); cov.append((np.abs(pull) < 1).mean()); rb.append(a[valid, 0].mean() / Ni - 1); rms.append(a[valid, 0].std() / Ni)
            hes.append(np.median(a[ok, 1]) / Ni); okf.append(ok.mean())
            jv = (a[:, 11] > 0) & np.isfinite(a[:, 10]) & (a[:, 10] > 0) & (a[:, 12] == 0)
            fp = (a[jv, 9] - fi) / a[jv, 10]
            fpm.append(fp.mean()); fpr.append(fp.std()); frms.append(a[a[:, 11] > 0, 9].std() / fi)
            E = a[:, 6] + a[:, 7]
            jv2 = a[:, 11] > 0
            erms = E[jv2].std() / Ni
            dcount += erms <= 0.95 * rms[-1]
            ep = (E[jv] - Ni) / a[jv, 8]
            ep = ep[np.isfinite(ep)]
            e_ok += len(ep) > 10 and abs(ep.mean()) <= 0.25 and ep.std() <= 1.25 and jv.mean() >= 0.90
            per_cell[(m, r["bin"], r["run"], r["level"])] = {"SB70": r["info"]["SB_win_70"], "SB_B": r["info"]["SB_win_B"], "pull_mean": float(pull.mean()), "pull_rms": float(pull.std()), "bias": float(rb[-1]),
                                                             "rms": float(rms[-1]), "hesse": float(hes[-1]), "ok": float(ok.mean()), "f_pull_mean": float(fp.mean()), "f_pull_rms": float(fp.std()),
                                                             "sf_rms": float(frms[-1]), "E_rms": float(erms), "cov": float(cov[-1])}
        D1[m] = bool(nok >= 10)
        stat[m] = {"D1_cells": int(nok), "sigma_toy": float(np.median(rms)), "bias_mean": float(np.mean(rb)), "pull_mean_rng": [float(min(pm)), float(max(pm))], "pull_rms_rng": [float(min(pr)), float(max(pr))],
                   "sf_rms_med": float(np.median(frms)), "d_cells": int(dcount), "E_D1_cells": int(e_ok)}
        dJ[m] = bool(dcount >= 8 and e_ok >= 10)
        L.append(f"| {m} | {nok} | {f(np.mean(pm), 2)} ({f(min(pm), 2)}…{f(max(pm), 2)}) | {u(np.median(pr), 2)} ({u(min(pr), 2)}…{u(max(pr), 2)}) | {u(np.median(cov), 2)} | {f(np.mean(rb))} ({f(min(rb))}…{f(max(rb))}) | "
                 f"{u(np.median(rms))} | {u(np.median(hes))} | {u(min(okf), 2)} | {f(np.median(fpm), 2)} / {u(np.median(fpr), 2)} | {u(np.median(frms))} | {dcount} | {e_ok} |")
    # ---- S5 ranking
    L += ["", "### S5. Ranking: Q = √(b_A² + σ_toy²)", "", "| ustawienie | D1 | D2 | ⟨|b_A|⟩ WP80 | σ_toy | Q | (d) lepszy od S? |", "|---|---|---|---|---|---|---|"]
    Q = {m: math.sqrt(bA[m] ** 2 + stat[m]["sigma_toy"] ** 2) for m in MS}
    for m in MS:
        L.append(f"| {m} {DESC[m]} | {'TAK' if D1[m] else 'NIE'} ({stat[m]['D1_cells']}/12) | {'TAK' if D2[m] else 'NIE'} | {u(bA[m])} | {u(stat[m]['sigma_toy'])} | {u(Q[m])} | {'TAK' if dJ[m] else 'NIE'} ({stat[m]['d_cells']}/12; D1(E) {stat[m]['E_D1_cells']}/12) |")
    adm = [m for m in MS if D1[m] and D2[m]]
    note = ""
    if not adm:
        adm = [m for m in MS if D1[m]]
        note = " (D2 niespełnione przez żadne ustawienie — wybór spośród spełniających D1)"
        if not adm:
            adm = MS; note = " (ani D1, ani D2 niespełnione przez żadne ustawienie — wybór wg samego Q)"
    best = min(adm, key=lambda m: Q[m])
    choice = best
    for m in SIMPLE:
        if m in adm and Q[m] <= 1.05 * Q[best]:
            choice = m
            break
    L += ["", f"**Dopuszczalne:** {', '.join(adm)}{note}. **Najmniejsze Q:** {best} ({u(Q[best])}). **Wybór (reguła prostoty 5 %):** **{choice}** — {DESC[choice]}; Q = {u(Q[choice])}; tryb (d) łączny: {'TAK' if dJ[choice] else 'NIE (zostaje tryb S)'}.", ""]
    sel[sp] = {"choice": choice, "best": best, "admissible": adm, "note": note, "Q": Q, "bA": bA, "D1": D1, "D2": D2, "toy": stat, "joint_better": dJ,
               "gate_passfail": {f"{m}|{et}": v for (m, et), v in gate.items()}, "per_cell": {"|".join(map(str, k)): v for k, v in per_cell.items()}}
    # ---- S6 per-cell for M0 and choice
    L += [f"### S6. Zabawki per komórka: M0 vs {choice}", "", "| pT | run | poziom | S/B okno (WP70) | S/B (B) | M0: obciąż. / ⟨pull⟩ / rms pull / rms N/N / HESSE | " + choice + ": obciąż. / ⟨pull⟩ / rms pull / rms N/N / HESSE | rms(SF_rel): M0 / " + choice + " | f pull rms: M0 / " + choice + " |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in toys:
        k0, k1 = per_cell[("M0", r["bin"], r["run"], r["level"])], per_cell[(choice, r["bin"], r["run"], r["level"])]
        g = lambda k: f"{f(k['bias'])} / {f(k['pull_mean'], 2)} / {u(k['pull_rms'], 2)} / {u(k['rms'])} / {u(k['hesse'])}"  # noqa: E731
        L.append(f"| {r['bin']} | {r['run']} | {r['level']} | {u(k0['SB70'])} | {u(k0['SB_B'])} | {g(k0)} | {g(k1)} | {u(k0['sf_rms'])} / {u(k1['sf_rms'])} | {u(k0['f_pull_rms'], 2)} / {u(k1['f_pull_rms'], 2)} |")
    L.append("")
open(os.path.join(HERE, "..", "tabele", "select.md"), "w").write("\n".join(L) + "\n")
json.dump(sel, open(O("select.json"), "w"), default=lambda o: o.item() if hasattr(o, "item") else str(o))
print({sp: (v["choice"], v["best"], v["admissible"], v["note"]) for sp, v in sel.items()})
