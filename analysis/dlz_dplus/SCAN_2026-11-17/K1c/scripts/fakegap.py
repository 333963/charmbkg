"""K1c POST HOC: domknięcie interpolacji ε_B(m) w miejscu bez sygnału („fałszywa przerwa” nad pikiem).
Ta sama regresja lokalnie liniowa (h = 25 MeV) co w tle prawdy V1/V2, ale przerwa wycięta wokół m_c = M0 + 130 MeV (V1: ±45 MeV; V2: ±70 MeV),
prawdziwy pik (±45 MeV) też wycięty z dopasowania. Porównanie: przewidziane n_all·ε_B vs zaobserwowane (dane − ogony sygnału) w fałszywej przerwie.
Wyjście: K1c/fakegap.md, out/fakegap.json.  Użycie: python fakegap.py"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fastfit as FF  # noqa: E402

K1C = os.path.dirname(HERE)
P = json.load(open(os.path.join(K1C, "out", "prep.json")))
S = json.load(open(os.path.join(K1C, "out", "k1c.json")))
RUNS, BN, TAGS, WPS = ["427929", "428770"], ["3-5", "5-8", "8-30"], ["H", "DLZ"], [0.5, 0.6, 0.7, 0.8, 0.9]
C = FF.C
MC = FF.M0 + 0.130


def smooth(npass, nall, tails, fit_mask, h=0.025):
    e = (npass - tails)[fit_mask] / nall[fit_mask]
    var = np.maximum(npass[fit_mask], 1.0) / nall[fit_mask] ** 2
    x = C[fit_mask]
    out = np.zeros(len(C))
    for i, x0 in enumerate(C):
        w = np.exp(-0.5 * ((x - x0) / h) ** 2) / var
        A = np.c_[np.ones(len(x)), x - x0]
        out[i] = np.linalg.solve(A.T @ (A * w[:, None]), A.T @ (w * e))[0]
    return np.clip(out, 0, None)


res, L = {}, []
L.append("# K1c POST HOC — domknięcie interpolacji ε_B(m) w fałszywej przerwie (m_c = M0 + 130 MeV)\n")
L.append("δ = (przewidziane − zaobserwowane)/zaobserwowane tło w fałszywej przerwie; „→ N” = δ·B_win/N_true, czyli jakie pozorne obciążenie N_S dałby taki sam błąd interpolacji pod prawdziwym pikiem.\n")
L.append("| run | pT | klas. | WP | V1 (±45): δ [%] ± stat | → N [%] | V2 (±70): δ [%] ± stat | → N [%] |")
L.append("|---|---|---|---|---|---|---|---|")
for run in RUNS:
    for bn in BN:
        S0 = S["S0"][f"{run}|{bn}"]
        nall = np.asarray(P["hall"][f"{run}|{bn}"], float)
        for t in TAGS:
            for wp in WPS:
                npass = np.asarray(P["hists"][f"{run}|{t}|{bn}|{wp}"], float)
                tails = S0 * np.asarray(P["sig"][f"{t}|{bn}|{wp}"])
                ty = S["toys"][f"{run}|{bn}|{t}|{wp}"]
                row = []
                for vn, half in (("V1", 0.045), ("V2", 0.070)):
                    gap = np.abs(C - MC) < half
                    fm = (~gap) & (np.abs(C - FF.M0) > 0.045) & (nall > 0)
                    pred = (nall * smooth(npass, nall, tails, fm))[gap].sum()
                    obs = (npass - tails)[gap].sum()
                    d = pred / obs - 1
                    de = np.sqrt(npass[gap].sum()) / obs
                    toN = d * ty[vn]["B_win"] / ty["N_true"]
                    res[f"{run}|{bn}|{t}|{wp}|{vn}"] = dict(delta=float(d), delta_err=float(de), obs=float(obs), pred=float(pred), implied_bias_N=float(toN), implied_err=float(de * ty[vn]["B_win"] / ty["N_true"]))
                    row.append(f"{100*d:+.2f} ± {100*de:.2f} | {100*toN:+.1f} ± {100*de*ty[vn]['B_win']/ty['N_true']:.1f}")
                L.append(f"| {run} | {bn} | {t} | {int(wp*100)} | " + " | ".join(row) + " |")
json.dump(res, open(os.path.join(K1C, "out", "fakegap.json"), "w"), indent=1)
open(os.path.join(K1C, "fakegap.md"), "w").write("\n".join(L) + "\n")
print("\n".join(l for l in L if "| 70 |" in l or l.startswith("| run")))
