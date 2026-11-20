"""LOWPT_calib: tabele i werdykt z out/scan.json (+ out/d0_qual.json, R3 verdict.json).  -> tables.md, verdict.json"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
CW = os.path.dirname(TOP)
J = json.load(open(os.path.join(TOP, "out", "scan.json")))
BN = ["3-5", "5-8", "8-30"]
WPS = [0.5, 0.6, 0.7, 0.8, 0.9]
RUNS = ["427929", "428770"]
L = []
P = L.append
V = dict(criteria={}, data={}, mc={}, dstar={}, d0={})


def cell(src, bn, o):
    N = np.array([J["cells"][f"{src}|{bn}|{wp}"]["fit"][str(o)]["N"] for wp in WPS])
    E = np.array([J["cells"][f"{src}|{bn}|{wp}"]["fit"][str(o)]["N_err"] for wp in WPS])
    F = [J["cells"][f"{src}|{bn}|{wp}"]["fit"][str(o)] for wp in WPS]
    if str(o) in J["cells"][f"{src}|{bn}|0.5"]["boot"]:
        B = np.array([J["cells"][f"{src}|{bn}|{wp}"]["boot"][str(o)]["N"] for wp in WPS]).T
        BE = np.array([J["cells"][f"{src}|{bn}|{wp}"]["boot"][str(o)]["err"] for wp in WPS]).T
        OK = np.array([J["cells"][f"{src}|{bn}|{wp}"]["boot"][str(o)]["ok"] for wp in WPS]).T.astype(bool)
    else:
        B = BE = OK = None
    return N, E, F, B, BE, OK


def slope(y, rel, x):
    w = 1 / rel ** 2
    ly = np.log(np.clip(y, 1e-9, None))
    xm = np.sum(w * x) / w.sum()
    return np.sum(w * (x - xm) * ly, axis=-1) / np.sum(w * (x - xm) ** 2)


EPS = {bn: np.array([J["eps"][f"{bn}|{wp}"]["eps"] for wp in WPS]) for bn in BN}
EPSR = {bn: np.array([J["eps"][f"{bn}|{wp}"]["ratio90_relerr"] for wp in WPS]) for bn in BN}
PAIRS = [(0, 4), (1, 4), (2, 4), (3, 4), (2, 3), (0, 3)]
PN = {(a, b): f"{int(100*WPS[a])}/{int(100*WPS[b])}" for a, b in PAIRS}

P("# LOWPT_calib — tabele (generowane przez scripts/report.py)\n")
P(f"Repliki bootstrapu: {J['nrep']}; wariant kształtu per bin: {J['fitvar']}. Dane: {json.dumps(J['info'])}\n")
P("## T1. Dane: y = N_fit/ε_MC per WP (pol3), błąd hesjanu; S/B w oknie ±45 MeV; χ²/ndf\n")
P("| run | pT | WP50 | WP60 | WP70 | WP80 | WP90 |\n|---|---|---|---|---|---|---|")
for run in RUNS:
    for bn in BN:
        N, E, F, B, BE, OK = cell(run, bn, 3)
        P(f"| {run} | {bn} | " + " | ".join(f"{N[i]/EPS[bn][i]:.0f} ± {E[i]/EPS[bn][i]:.0f} (S/B {F[i]['SB']:.3f}; {F[i]['chi2']:.0f}/{F[i]['ndf']}{'' if F[i]['valid'] else ' !'})" for i in range(5)) + " |")
P("\n## T2. Dane: SF względny ρ(a/b) = [N/ε](WPa)/[N/ε](WPb), pol3; błąd = bootstrap sparowany ⊕ stat. MC ε (dla par z WP90); w nawiasie pol4 i pol2 (nominalne)\n")
P("| run | pT | " + " | ".join(f"ρ({PN[p]})" for p in PAIRS) + " | nachylenie d ln y/dε | udane repliki |\n|---|---|" + "---|" * (len(PAIRS) + 2))
for run in RUNS:
    for bn in BN:
        N, E, F, B, BE, OK = cell(run, bn, 3)
        N4, _, _, B4, _, OK4 = cell(run, bn, 4)
        N2 = cell(run, bn, 2)[0]
        y, yb, y4, y2 = N / EPS[bn], B / EPS[bn], N4 / EPS[bn], N2 / EPS[bn]
        row = []
        V["data"][f"{run}|{bn}"] = {}
        for a, b in PAIRS:
            ok = OK[:, a] & OK[:, b]
            rho = y[a] / y[b]
            sb = np.std(yb[ok, a] / yb[ok, b], ddof=1)
            smc = rho * np.hypot(EPSR[bn][a], 0) if b == 4 else 0.0
            s = float(np.hypot(sb, smc))
            ok4 = OK4[:, a] & OK4[:, b]
            s4 = float(np.std(B4[ok4, a] / B4[ok4, b], ddof=1) * EPS[bn][b] / EPS[bn][a])
            bar = rho * np.sqrt(abs((E[a] / N[a]) ** 2 - (E[b] / N[b]) ** 2))
            V["data"][f"{run}|{bn}"][PN[(a, b)]] = dict(rho=float(rho), err=s, err_boot=float(sb), err_mc=float(smc), sens2s=float(2 * s / rho), rho_pol4=float(y4[a] / y4[b]), err_pol4=s4,
                                                    rho_pol2=float(y2[a] / y2[b]), f_ok=float(ok.mean()), err_barlow=float(bar), pull_vs1=float((rho - 1) / s))
            row.append(f"{rho:.3f} ± {s:.3f} ({y4[a]/y4[b]:.3f}; {y2[a]/y2[b]:.3f})")
        ok = OK.all(1)
        rel = E / N
        sl = float(slope(y, rel, EPS[bn]))
        sle = float(np.std(slope(yb[ok], rel, EPS[bn]), ddof=1))
        V["data"][f"{run}|{bn}"]["slope"] = dict(val=sl, err=sle, f_ok=float(ok.mean()))
        V["data"][f"{run}|{bn}"]["nominal_valid"] = bool(all(f["valid"] for f in F))
        P(f"| {run} | {bn} | " + " | ".join(row) + f" | {sl:+.2f} ± {sle:.2f} | {100*ok.mean():.0f} % |")
P("\n## T3. Czułość 2σ na odchylenie względnego SF [%] (pol3): stat.; w nawiasie stat ⊕ |ρ(pol3) − ρ(pol4)|; σ_boot/σ_Barlow\n")
P("| run | pT | ρ(70/90) | ρ(50/90) | ρ(70/80) | ρ(50/80) | σ_boot/σ_Barlow (70/90) |\n|---|---|---|---|---|---|---|")
for run in RUNS:
    for bn in BN:
        d = V["data"][f"{run}|{bn}"]
        f = lambda k: f"{200*d[k]['err']/d[k]['rho']:.1f} ({200*np.hypot(d[k]['err'], d[k]['rho']-d[k]['rho_pol4'])/d[k]['rho']:.1f})"  # noqa: E731
        P(f"| {run} | {bn} | {f('70/90')} | {f('50/90')} | {f('70/80')} | {f('50/80')} | {d['70/90']['err_boot']/max(d['70/90']['err_barlow'],1e-9):.2f} |")

# ---- MC
P("\n## T4. MC mb z pełnym tłem: obciążenie b = N_fit/N_true per WP (pol3; błąd z bootstrapu po zdarzeniach, prawda per replika)\n")
P("| wariant | pT | WP50 | WP60 | WP70 | WP80 | WP90 |\n|---|---|---|---|---|---|---|")
MCV = {}
for var in ("mc_nom", "mc_inj"):
    M = J["mc"][var]
    nt, ntu, ntb = np.array(M["nt"]), np.array(M["nt_unw"]), np.array(M["ntb"])
    for bi, bn in enumerate(BN):
        N, E, F, B, BE, OK = cell(var, bn, 3)
        b = N / nt[bi]
        be = [np.std((B[:, i] / ntb[:, bi, i])[OK[:, i]], ddof=1) for i in range(5)]
        P(f"| {var} | {bn} | " + " | ".join(f"{b[i]:.3f} ± {be[i]:.3f}" for i in range(5)) + " |")
        MCV[(var, bn)] = (N, E, B, BE, OK, nt[bi], ntu[bi], ntb[:, bi])
P("\n## T5. Wstrzyknięcie: ρ_rec vs ρ_true (pol3; w nawiasie pol4). σ_stat = rozrzut ρ_rec przy stałej ε (statystyka MC ≈ 427929); obciążenie = ρ_rec/ρ_true − 1, błąd z prawdą per replika; pull width względem σ_Barlow (hesjany)\n")
P("| wariant | pT | para | ρ_true | ρ_rec ± σ_stat | obciążenie [%] | (ρ_rec − ρ_true)/σ_stat | pull width (odporny) | udane |\n|---|---|---|---|---|---|---|---|---|")
for var in ("mc_nom", "mc_inj"):
    for bn in BN:
        N, E, B, BE, OK, nt, ntu, ntb = MCV[(var, bn)]
        N4, _, _, B4, _, OK4 = cell(var, bn, 4)
        V["mc"][f"{var}|{bn}"] = {}
        for a, b in PAIRS[:4] + [(2, 3)]:
            ok = OK[:, a] & OK[:, b]
            rt = (nt[a] / ntu[a]) / (nt[b] / ntu[b])
            rr = (N[a] / ntu[a]) / (N[b] / ntu[b])
            rr4 = (N4[a] / ntu[a]) / (N4[b] / ntu[b])
            rb = (B[ok, a] / ntu[a]) / (B[ok, b] / ntu[b])
            sst = float(np.std(rb, ddof=1))
            bias_r = (B[ok, a] / ntb[ok, a]) / (B[ok, b] / ntb[ok, b])
            bias = rr / rt - 1
            sbias = float(np.std(bias_r, ddof=1))
            sig_r = rb * np.sqrt(np.abs((BE[ok, a] / B[ok, a]) ** 2 - (BE[ok, b] / B[ok, b]) ** 2))
            pl = (rb - rr) / sig_r
            pw = float(np.std(pl, ddof=1))
            q = np.percentile(pl, [16, 84])
            pwr = float((q[1] - q[0]) / 2)
            V["mc"][f"{var}|{bn}"][PN[(a, b)]] = dict(rho_true=float(rt), rho_rec=float(rr), rho_rec_pol4=float(rr4), sigma_stat=sst, bias=float(bias), bias_err=sbias,
                                                     pull_rec=float((rr - rt) / sst), pull_width=pw, pull_width_robust=pwr, f_ok=float(ok.mean()),
                                                     barlow_nom=float(rr * np.sqrt(abs((E[a] / N[a]) ** 2 - (E[b] / N[b]) ** 2))))
            P(f"| {var} | {bn} | {PN[(a, b)]} | {rt:.3f} | {rr:.3f} ± {sst:.3f} ({rr4:.3f}) | {100*bias:+.1f} ± {100*sbias:.1f} | {(rr-rt)/sst:+.2f} | {pw:.2f} ({pwr:.2f}) | {100*ok.mean():.0f} % |")
P(f"\nWagi wstrzyknięcia (sygnał prawdy ze score ≥ próg WP70): " + ", ".join(f"{bn}: {J['mc']['w_inj|'+bn]:.3f}" for bn in BN) + f"; zdarzeń MC po WP90: {J['mc']['n_events']}\n")

# ---- kontrola D*
P("## T6. Kontrola: SF względny ze skanu D⁺ vs z tagu D* dla D⁰ (R3, fit kształtu; błąd z pełnej kowariancji R3). Różne gatunki i klasyfikatory — tylko rząd wielkości\n")
P("| run | pT | para | D⁺ skan ρ | D⁰ tag D* SF(a)/SF(b) | różnica [σ] |\n|---|---|---|---|---|---|")
R3 = json.load(open(os.path.join(CW, "FIT_scorebins", "R3_2026-11-18", "verdict.json")))["data"]
for run in RUNS:
    for bn in ("5-8", "8-30"):
        d = R3[f"{run}_{bn}"]
        sf, cov, wps = np.array(d["SF"]), np.array(d["SF_cov_tot"]), d["wps"]
        for pa, pb in (("70", "80"), ("50", "80")):
            i, j = wps.index(pa), wps.index(pb)
            r = sf[i] / sf[j]
            g = np.zeros(len(sf))
            g[i], g[j] = 1 / sf[j], -sf[i] / sf[j] ** 2
            e = float(np.sqrt(g @ cov @ g))
            x = V["data"][f"{run}|{bn}"][f"{pa}/{pb}"]
            z = (x["rho"] - r) / np.hypot(x["err"], e)
            V["dstar"][f"{run}|{bn}|{pa}/{pb}"] = dict(dplus=x["rho"], dplus_err=x["err"], dstar=float(r), dstar_err=e, z=float(z))
            P(f"| {run} | {bn} | {pa}/{pb} | {x['rho']:.3f} ± {x['err']:.3f} | {r:.3f} ± {e:.3f} | {z:+.1f} |")

# ---- D0 jakościowo
p0 = os.path.join(TOP, "out", "d0_qual.json")
if os.path.exists(p0):
    D0 = json.load(open(p0))
    P("\n## T7. D⁰ inkluzywne (jakościowo; okno 1,75–2,00 GeV, bez bootstrapu; błąd ρ ze wzoru Barlowa)\n")
    P("MC fd: R = odbicie(w oknie fitu)/sygnał, ε odbicia, odbicie/sygnał w ±40 MeV:\n")
    P("| pT | WP50 | WP70 | WP90 |\n|---|---|---|---|")
    for bn in BN:
        P(f"| {bn} | " + " | ".join(f"R {D0['mc'][f'{bn}|{wp}']['R']:.2f}, ε_odb {D0['mc'][f'{bn}|{wp}']['eps_refl']:.2f}, okno {D0['mc'][f'{bn}|{wp}']['refl_win_over_sig_win']:.2f}" for wp in (0.5, 0.7, 0.9)) + " |")
    P("\n| run | pT | wariant | y(WP70) | y(WP90) | ρ(70/90) ± σ_Barlow | ρ(50/90) |\n|---|---|---|---|---|---|---|")
    for run in RUNS:
        for bn in BN:
            for lab in ("R_MC|pol3", "R_MC|pol2", "R0|pol3", "R1.3|pol3", "Rfree|pol2"):
                f = {wp: D0["data"][f"{run}|{bn}|{wp}"]["fits"][lab] for wp in WPS}
                e = {wp: D0["mc"][f"{bn}|{wp}"]["eps"] for wp in WPS}
                y = {wp: f[wp]["N"] / e[wp] for wp in WPS}
                rel = {wp: f[wp]["N_err"] / max(f[wp]["N"], 1e-9) for wp in WPS}
                r79, r59 = y[0.7] / y[0.9], y[0.5] / y[0.9]
                s79 = r79 * np.sqrt(abs(rel[0.7] ** 2 - rel[0.9] ** 2))
                V["d0"][f"{run}|{bn}|{lab}"] = dict(y70=y[0.7], y70_err=f[0.7]["N_err"] / e[0.7], y90=y[0.9], y90_err=f[0.9]["N_err"] / e[0.9], rho7090=r79, rho7090_barlow=float(s79), rho5090=r59,
                                                    rel70=rel[0.7], rel90=rel[0.9])
                P(f"| {run} | {bn} | {lab} | {y[0.7]:.0f} ± {f[0.7]['N_err']/e[0.7]:.0f} | {y[0.9]:.0f} ± {f[0.9]['N_err']/e[0.9]:.0f} | {r79:.2f} ± {s79:.2f} | {r59:.2f} |")

# ---- kryteria
d = V["data"]["427929|3-5"]["70/90"]
KA = bool(d["sens2s"] <= 0.15 and V["data"]["427929|3-5"]["nominal_valid"] and d["f_ok"] >= 0.8)
m = V["mc"]["mc_inj|3-5"]["70/90"]
KB1 = bool(abs(m["rho_rec"] - 0.90) <= m["sigma_stat"])
KB1s = bool(abs(m["rho_rec"] - 0.90) <= 0.90 * m["bias_err"])
KB2 = bool(abs(m["rho_rec"] / 0.90 - 1) < 0.05)
KB = bool(KB1 and KB2 and m["f_ok"] >= 0.8)
V["criteria"] = dict(K_A=KA, K_A_value=d["sens2s"], K_B=KB, K_B_within_1sigma_stat=KB1, K_B_within_1sigma_bias_err=KB1s, K_B_bias_lt_5pct=KB2, K_B_rho_rec=m["rho_rec"], K_B_sigma_stat=m["sigma_stat"],
                     K_B_bias=m["bias"], K_B_bias_err=m["bias_err"], verdict="UŻYTECZNA" if KA and KB else "NIEUŻYTECZNA")
P("\n## T8. Kryteria (KRYTERIA.md §2)\n")
P(f"- K_A: czułość 2σ ρ(70/90), 427929, 3–5 GeV = {100*d['sens2s']:.1f} % (próg ≤ 15 %) → {'SPEŁNIONE' if KA else 'NIESPEŁNIONE'}")
P(f"- K_B: wstrzyknięcie 3–5 GeV: ρ_rec = {m['rho_rec']:.3f} ± {m['sigma_stat']:.3f} (σ_stat), obciążenie {100*m['bias']:+.1f} ± {100*m['bias_err']:.1f} %; "
  f"w 1σ_stat: {'tak' if KB1 else 'NIE'}; w 1σ błędu obciążenia: {'tak' if KB1s else 'NIE'}; |obciążenie| < 5 %: {'tak' if KB2 else 'NIE'} → {'SPEŁNIONE' if KB else 'NIESPEŁNIONE'}")
P(f"- **Werdykt: {V['criteria']['verdict']}**")
open(os.path.join(TOP, "tables.md"), "w").write("\n".join(L) + "\n")
json.dump(V, open(os.path.join(TOP, "verdict.json"), "w"), indent=1)
print("\n".join(L))
