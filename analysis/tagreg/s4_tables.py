#!/usr/bin/env python3
"""s4: tabele części C z s3.json (strony A, C i średnia ważona N) + sprawdzenie klasy `none` na 2 dodatkowych łańcuchach z trg_data_428770.npz."""
import json, os, numpy as np, common as C
S = json.load(open(os.path.join(C.OUT, "s3.json"))); L = []
def comb(acc, sub="all"):
    """suma stron: ε = Σ(N_SR − acc_SR)/Σ(N_tot − acc_tot) per wariant normalizacji."""
    o = {}
    sides = [s for s in ("A", "C") if s in acc and sub in acc[s]]
    if not sides: return None
    N = sum(acc[s][sub]["N_tot"] for s in sides); NS = sum(acc[s][sub]["N_SR"] for s in sides)
    o["N"] = N; o["raw"] = NS / N
    for nm in ("n10", "n045", "n20"):
        a_sr = sum(acc[s][sub][nm]["N_acc_SR"] for s in sides); a_t = sum(acc[s][sub][nm]["N_acc_tot"] for s in sides)
        err = np.sqrt(sum((acc[s][sub][nm]["eps_err"] * (acc[s][sub]["N_tot"] - acc[s][sub][nm]["N_acc_tot"])) ** 2 for s in sides)) / (N - a_t)
        o[nm] = ((NS - a_sr) / (N - a_t), err, a_sr / NS, a_t / N)
    return o
def fmt(o):
    e, de = o["n10"][:2]; lo = min(o[k][0] for k in ("n10", "n045", "n20")); hi = max(o[k][0] for k in ("n10", "n045", "n20"))
    return f"{o['N']:.0f} | {o['raw']:.4f} | {e:.4f} ± {de:.4f} | {o['n045'][0]:.4f} | {o['n20'][0]:.4f} ± {o['n20'][1]:.4f} | {o['n10'][2]:.4f} | {o['n10'][3]:.4f} | [{lo:.3f}; {hi:.3f}]"
H = "| run | strumień | próbka | pula | N tagów (A + C) | surowy ułamek Δ ≤ 0,02 | ε_true (norm. Δ > 0,10) ± stat | ε_true (norm. > 0,045) | ε_true (norm. > 0,20) | f_acc w SR | udział przypadkowych we wszystkich tagach | rozstęp wariantów |\n|---|---|---|---|---|---|---|---|---|---|---|---|"
L.append("## C2. ε_true cięcia ξ_cal ≤ ξ_AFP + 0,02 (strony A + C)\n\n" + H)
for r in ("427929", "428770"):
    for key, acc in S[r].items():
        o = comb(acc)
        if o is None: continue
        st, smp, pool = key.split("|"); L.append(f"| {r} | {st} | {smp} | {pool[5:]} | " + fmt(o) + " |")
L.append("\n## C2a. Strony osobno (próbka „z kand. D”, pierwsza pula)\n\n| run | strumień | strona | N | surowy | ε_true n10 ± stat | n045 | n20 |\n|---|---|---|---|---|---|---|---|")
for r in ("427929", "428770"):
    for key, acc in S[r].items():
        st, smp, pool = key.split("|")
        if smp != "z kand. D" or pool not in ("pool=all", "pool=nonAFP"): continue
        for s in ("A", "C"):
            if s in acc:
                a = acc[s]["all"]; L.append(f"| {r} | {st} | {s} | {a['N_tot']:.0f} | {a['eps_raw']:.4f} | {a['n10']['eps_true']:.4f} ± {a['n10']['eps_err']:.4f} | {a['n045']['eps_true']:.4f} | {a['n20']['eps_true']:.4f} |")
L.append("\n## C3. Zależności (strony A + C; ε_true przy normalizacji Δ > 0,10 [w nawiasie > 0,045 / > 0,20])\n\n| run | strumień | próbka | przekrój | N | surowy | ε_true n10 ± stat | n045 | n20 | f_acc SR |\n|---|---|---|---|---|---|---|---|---|---|")
for r in ("427929", "428770"):
    for key, acc in S[r].items():
        st, smp, pool = key.split("|")
        if pool not in ("pool=all", "pool=nonAFP") or smp not in ("z kand. D", "wszystkie"): continue
        subs = [k for k in acc.get("A", {}) if k not in ("all", "n_prot", "n_cen")]
        for sub in subs:
            o = comb(acc, sub)
            if o is None or o["N"] < 30: continue
            mu = ""
            if sub.startswith("LB"): mu = f" (⟨μ⟩ = {np.mean([acc[s][sub]['mu'] for s in ('A','C') if s in acc]):.4f})"
            L.append(f"| {r} | {st} | {smp} | {sub}{mu} | {o['N']:.0f} | {o['raw']:.4f} | {o['n10'][0]:.4f} ± {o['n10'][1]:.4f} | {o['n045'][0]:.4f} | {o['n20'][0]:.4f} | {o['n10'][2]:.4f} |")
# rozkład Δ: dane i szablon (A + C), próbka z kand. D
E_ = S["edges"]
L.append("\n## C1. Rozkład Δ (zdarzenia tagowane z kandydatem D, A + C): dane | przypadkowe (szablon, norm. Δ > 0,10)\n")
keys = {"427929": "all|z kand. D|pool=all", "428770": "all|z kand. D|pool=nonAFP"}
L.append("| Δ od | Δ do | 427929 dane | 427929 przyp. | 428770 dane | 428770 przyp. |\n|---|---|---|---|---|---|")
hh = {r: (sum(np.array(S[r][k][s]["all"]["data_hist"]) for s in ("A", "C")), sum(np.array(S[r][k][s]["all"]["acc_hist"]) for s in ("A", "C"))) for r, k in keys.items()}
for i in range(len(E_) - 1):
    lo = "−∞" if E_[i] is None or E_[i] == -np.inf or str(E_[i]) == "-Infinity" else f"{float(E_[i]):.3f}"; hi = "+∞" if str(E_[i + 1]) in ("Infinity", "inf") else f"{float(E_[i+1]):.3f}"
    L.append(f"| {lo} | {hi} | {hh['427929'][0][i]:.0f} | {hh['427929'][1][i]:.1f} | {hh['428770'][0][i]:.0f} | {hh['428770'][1][i]:.1f} |")
# klasa none: łańcuchy pt6/pt8 z trg_data
d = C.load(428770); cl, b = C.classes(d); T0, T1, T2 = C.regions(d)
z = np.load(os.path.join(C.BT, "trg_data_428770.npz"), allow_pickle=False)
o = np.argsort(z["event"]); p = np.searchsorted(z["event"][o], d["evno"]); ok = z["event"][o][np.clip(p, 0, len(o) - 1)] == d["evno"]
L.append(f"\n## A dodatek: klasa `none` a łańcuchy sptrk_pt6/pt8_L1MBTS_2 (trg_data_428770.npz; dopasowanie zdarzeń {ok.mean():.4f})\n")
for c in ("HLT_mb_sptrk_pt6_L1MBTS_2", "HLT_mb_sptrk_pt8_L1MBTS_2"):
    v = np.zeros(len(ok), bool); v[ok] = z["tp_" + c][o][p[ok]] == 1
    u = np.unique(z["ps_" + c][z["ps_" + c] > 0])
    L.append(f"- {c}: preskale {u.tolist()}; w klasie none (GRL): {(cl['none'] & v).sum()} z {cl['none'].sum()}; none ∧ T2: {(cl['none'] & T2 & v).sum()} z {(cl['none'] & T2).sum()}; none ∧ T2 ∧ D⁰: {(cl['none'] & T2 & d['hasD0'] & v).sum()} z {(cl['none'] & T2 & d['hasD0']).sum()}")
n = cl["none"] & T2
L.append(f"- none ∧ T2: TBP L1_AFP {(n & d['TBP_L1_AFP_A_OR_C']).mean():.3f}→{(n & d['TBP_L1_AFP_A_OR_C']).sum()/n.sum():.3f}, TBP L1_MBTS_1 {(n & d['TBP_L1_MBTS_1']).sum()/n.sum():.3f}, TBP L1_MBTS_2 {(n & d['TBP_L1_MBTS_2']).sum()/n.sum():.3f}, TBP L1_J12 {(n & d['TBP_L1_J12']).sum()/n.sum():.3f}; ⟨nSelBL⟩ {d['nsel'][n].mean():.1f}; z D⁰: TBP L1_AFP {(n & d['hasD0'] & d['TBP_L1_AFP_A_OR_C']).sum()/(n & d['hasD0']).sum():.3f}, TBP MBTS_2 {(n & d['hasD0'] & d['TBP_L1_MBTS_2']).sum()/(n & d['hasD0']).sum():.3f}, ⟨nSelBL⟩ {d['nsel'][n & d['hasD0']].mean():.1f}")
open(os.path.join(C.OUT, "s4_tables.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
