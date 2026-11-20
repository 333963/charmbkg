"""K1c: tabele i werdykt z out/k1c.json → K1c/tables.md, K1c/verdict.json. Użycie: python report.py"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K1C = os.path.dirname(HERE)
S = json.load(open(os.path.join(K1C, "out", "k1c.json")))
P = json.load(open(os.path.join(K1C, "out", "prep.json")))
RUNS, BN, TAGS, WPS, ORD = ["427929", "428770"], ["3-5", "5-8", "8-30"], ["H", "DLZ"], [0.5, 0.6, 0.7, 0.8, 0.9], [2, 3, 4]
VAR = ["F", "22", "33", "44"]
VN = {"F": "F-test", "22": "pol2/pol2", "33": "pol3/pol3", "44": "pol4/pol4"}
L = []
p = L.append
pc = lambda v, e: f"{100*v:+.1f} ± {100*e:.1f}"  # noqa: E731
pm = lambda v, e: f"{v:.3f} ± {e:.3f}"  # noqa: E731


def gls_const(y, Cv):
    Ci = np.linalg.inv(Cv)
    one = np.ones(len(y))
    v = 1.0 / (one @ Ci @ one)
    c = v * (one @ Ci @ y)
    d = y - c
    return float(c), float(np.sqrt(v)), float(d @ Ci @ d)


def gls_slope(y, Cv, x):
    A = np.c_[np.ones(len(y)), np.asarray(x) - 0.7]
    Ci = np.linalg.inv(Cv)
    cov = np.linalg.inv(A.T @ Ci @ A)
    be = cov @ (A.T @ Ci @ y)
    J = np.array([-be[1] / be[0] ** 2, 1 / be[0]])
    return float(be[1] / be[0]), float(np.sqrt(J @ cov @ J))


def order(run, bn, t, wp, var):
    fv = S["fitvar"][f"{run}|{bn}"]
    return S["nominal"][f"{run}|{bn}|{t}|{wp}"][fv]["sel"] if var == "F" else int(var[0] if t == "H" else var[1])


V = {}
T_gain, T_bias, T_scan, T_R, T_ft, T_bkg = [], [], [], [], [], []
for run in RUNS:
    V[run] = {"bins": {}}
    for bn in BN:
        fv = S["fitvar"][f"{run}|{bn}"]
        cell = {"fitvar": fv, "gain": {}, "scan": {}}
        for var in VAR:
            # --- skan: macierz replik N dla 2 klasyfikatorów × 5 WP
            N, eps, rel, cols, okc, erc = [], [], [], [], [], []
            for t in TAGS:
                for wp in WPS:
                    o = order(run, bn, t, wp, var)
                    f = S["nominal"][f"{run}|{bn}|{t}|{wp}"][fv]["fits"][str(o)]
                    bt = np.array(S["boot"][f"{run}|{bn}|{t}|{wp}"][str(o)], float)
                    N.append(f["N"]); rel.append(f["rel"]); eps.append(P["shapes"][f"{t}|{bn}|{wp}"]["eps_mc"])
                    cols.append(bt[:, 0]); erc.append(bt[:, 1]); okc.append((bt[:, 4] > 0) & (bt[:, 0] > 0) & (bt[:, 1] > 0))
            N, eps, rel = np.array(N), np.array(eps), np.array(rel)
            Nr, Er, ok = np.column_stack(cols), np.column_stack(erc), np.column_stack(okc)
            nW = len(WPS)
            i7 = WPS.index(0.7)
            # zysk WP70 (def. K1), sparowany
            o7 = ok[:, i7] & ok[:, nW + i7]
            g = 1 - rel[nW + i7] / rel[i7]
            grep = 1 - (Er[o7, nW + i7] / Nr[o7, nW + i7]) / (Er[o7, i7] / Nr[o7, i7])
            oH, oD = order(run, bn, "H", 0.7, var), order(run, bn, "DLZ", 0.7, var)
            tH, tD = S["toys"][f"{run}|{bn}|H|0.7"], S["toys"][f"{run}|{bn}|DLZ|0.7"]
            q = {}
            for vn in ("V1", "V2"):
                a, b = tH[vn]["fit"][str(oH)], tD[vn]["fit"][str(oD)]
                q[vn] = dict(bH=a["b"], bH_err=a["b_err"], bD=b["b"], bD_err=b["b_err"], wH=a["pull_w"], wD=b["pull_w"], wH_err=a["pull_w_err"], wD_err=b["pull_w_err"], okH=a["f_ok"], okD=b["f_ok"],
                             pmH=a["pull_mean"], pmD=b["pull_mean"], asH=a["asimov"], asD=b["asimov"],
                             g_corr=1 - (rel[nW + i7] * b["b"]) / (rel[i7] * a["b"]), g_toy=1 - (b["rms"] / b["b"]) / (a["rms"] / a["b"]))
            cell["gain"][var] = dict(oH=oH, oD=oD, g=float(g), g_err=float(np.std(grep)), n_ok=int(o7.sum()), relH=float(rel[i7]), relD=float(rel[nW + i7]), NH=float(N[i7]), ND=float(N[nW + i7]),
                                     Nratio=float(N[nW + i7] / N[i7]), Nratio_err=float(np.std(Nr[o7, nW + i7] / Nr[o7, i7])),
                                     sigboot_over_hesse_H=float(np.std(Nr[ok[:, i7], i7]) / (rel[i7] * N[i7])), sigboot_over_hesse_D=float(np.std(Nr[ok[:, nW + i7], nW + i7]) / (rel[nW + i7] * N[nW + i7])), **q)
            # skan
            oa = ok.all(1)
            y = N / eps
            Y = Nr[oa] / eps
            Cy = np.cov(Y.T)
            sc = {"n_ok": int(oa.sum()), "orders": {t: [order(run, bn, t, wp, var) for wp in WPS] for t in TAGS}}
            for ti, t in enumerate(TAGS):
                sl = slice(ti * nW, (ti + 1) * nW)
                c, ce, chi = gls_const(y[sl], Cy[sl, sl])
                s, se = gls_slope(y[sl], Cy[sl, sl], WPS)
                sc[t] = dict(y=y[sl].tolist(), y_err=np.sqrt(np.diag(Cy[sl, sl])).tolist(), mean=c, mean_err=ce, chi2_flat=chi, slope=s, slope_err=se)
            R = y[nW:] / y[:nW]
            Rrep = Y[:, nW:] / Y[:, :nW]
            CR = np.cov(Rrep.T)
            c, ce, chi = gls_const(R, CR)
            sc["R"] = dict(R=R.tolist(), R_err=np.sqrt(np.diag(CR)).tolist(), mean=c, mean_err=ce, chi2_const=chi)
            # korekta: V1 toye; systematyka |V1 − V2| z Asimova
            r1, r1e, dv, bHl, bDl = [], [], [], [], []
            for wi, wp in enumerate(WPS):
                a = S["toys"][f"{run}|{bn}|H|{wp}"]
                b = S["toys"][f"{run}|{bn}|DLZ|{wp}"]
                oh, od = order(run, bn, "H", wp, var), order(run, bn, "DLZ", wp, var)
                a1, b1 = a["V1"]["fit"][str(oh)], b["V1"]["fit"][str(od)]
                r1.append(b1["b"] / a1["b"]); r1e.append(b1["b"] / a1["b"] * np.hypot(b1["b_err"] / b1["b"], a1["b_err"] / a1["b"]))
                dv.append(b1["asimov"] / a1["asimov"] - b["V2"]["fit"][str(od)]["asimov"] / a["V2"]["fit"][str(oh)]["asimov"])
                bHl.append(a1["b"]); bDl.append(b1["b"])
            r1, r1e, dv = np.array(r1), np.array(r1e), np.array(dv)
            Rc = R / r1
            CRc = CR / np.outer(r1, r1) + np.diag(Rc ** 2 * ((r1e / r1) ** 2 + (dv / r1) ** 2))
            c, ce, chi = gls_const(Rc, CRc)
            sc["Rcorr"] = dict(R=Rc.tolist(), R_err=np.sqrt(np.diag(CRc)).tolist(), mean=c, mean_err=ce, chi2_const=chi, chi2_vs1=float((Rc - 1) @ np.linalg.inv(CRc) @ (Rc - 1)),
                               bH=bHl, bD=bDl, dV=dv.tolist())
            # y skorygowane: płaskość
            for ti, t in enumerate(TAGS):
                sl = slice(ti * nW, (ti + 1) * nW)
                bb = np.array(bHl if t == "H" else bDl)
                yc = y[sl] / bb
                Cc = Cy[sl, sl] / np.outer(bb, bb)
                c2, ce2, chi2 = gls_const(yc, Cc)
                s2, se2 = gls_slope(yc, Cc, WPS)
                sc[t].update(chi2_flat_corr=chi2, slope_corr=s2, slope_corr_err=se2)
            cell["scan"][var] = sc
            gg = cell["gain"][var]
            T_gain.append((run, bn, var, gg))
            T_scan.append(f"| {run} | {bn} | {VN[var]} | " + " | ".join(
                f"{t}: " + " / ".join(f"{a:.0f}±{e:.0f}" for a, e in zip(sc[t]["y"], sc[t]["y_err"])) + f"; χ² {sc[t]['chi2_flat']:.1f} → {sc[t]['chi2_flat_corr']:.1f}; nach. {sc[t]['slope']:+.2f}±{sc[t]['slope_err']:.2f} → {sc[t]['slope_corr']:+.2f}" for t in TAGS)
                + f" | {'/'.join(map(str, sc['orders']['H']))} ; {'/'.join(map(str, sc['orders']['DLZ']))} | {sc['n_ok']} |")
            T_R.append(f"| {run} | {bn} | {VN[var]} | " + " | ".join(f"{a:.3f}±{ae:.3f} → {x:.3f}±{xe:.3f}" for a, ae, x, xe in zip(sc["R"]["R"], sc["R"]["R_err"], sc["Rcorr"]["R"], sc["Rcorr"]["R_err"]))
                       + f" | {pm(sc['R']['mean'], sc['R']['mean_err'])} | {pm(sc['Rcorr']['mean'], sc['Rcorr']['mean_err'])} | {sc['Rcorr']['chi2_vs1']:.1f} |")
        V[run]["bins"][bn] = cell
        for t in TAGS:
            for wp in WPS:
                n = S["nominal"][f"{run}|{bn}|{t}|{wp}"][fv]
                f, ft = n["fits"], n["ftest"]
                T_ft.append(f"| {run} | {bn} | {t} | {int(wp*100)} | " + " | ".join(f"{f[o]['N']:.0f} ± {f[o]['N_err']:.0f} ({100*f[o]['rel']:.1f} %); {f[o]['chi2']/f[o]['ndf']:.2f}{'*' if f[o]['at_limit'] else ''}{'' if f[o]['valid'] else '!'}" for o in "234")
                            + f" | {ft['2->3']['p']:.3g} | {ft['3->4']['p']:.3g} | {ft['2->4']['p']:.3g} | **pol{n['sel']}** |")
                ty = S["toys"][f"{run}|{bn}|{t}|{wp}"]

                def cellb(vn, o):
                    q = ty[vn]["fit"][o]
                    if "b" in q:
                        return f"{q['b']:.3f}±{q['b_err']:.3f} (As {q['asimov']:.3f}; w {q['pull_w']:.2f}; pull {q['pull_mean']:+.2f}; ok {q['f_ok']:.2f})"
                    return f"As {q['asimov']:.3f}"
                T_bias.append(f"| {run} | {bn} | {t} | {int(wp*100)} | pol{n['sel']} | " + " | ".join(cellb("V1", o) for o in "234") + " | " + " | ".join(cellb("V2", o) for o in "234") + " |")
                T_bkg.append(f"| {run} | {bn} | {t} | {int(wp*100)} | {ty['N_true']:.0f} | {ty['V1']['B_win']:.0f} | {ty['V2']['B_win']:.0f} | {ty['V1']['SB_true']:.3f} | {ty['V1']['chi2_sb']:.2f} | {ty['V2']['chi2_sb']:.2f} | "
                             f"{ty['V1']['data_win_minus_SplusB']:+.0f} | {ty['V2']['data_win_minus_SplusB']:+.0f} | {ty['sig_in_window_frac']:.4f} |")
    bb = V[run]["bins"]
    summ = {}
    for var in VAR:
        gs = [bb[b]["gain"][var] for b in BN]
        summ[var] = dict(mean=float(np.mean([x["g"] for x in gs])), err=float(np.sqrt(np.sum([x["g_err"] ** 2 for x in gs])) / 3),
                         mean_corr_V1=float(np.mean([x["V1"]["g_corr"] for x in gs])), mean_corr_V2=float(np.mean([x["V2"]["g_corr"] for x in gs])))
    c1 = summ["F"]["mean"] >= 0.10

    def c23(vn):
        c2 = {b: bool(abs(bb[b]["gain"]["F"][vn]["bH"] - 1) < 0.05 and abs(bb[b]["gain"]["F"][vn]["bD"] - 1) < 0.05 and bb[b]["gain"]["F"][vn]["okH"] >= 0.8 and bb[b]["gain"]["F"][vn]["okD"] >= 0.8) for b in BN}
        c3 = {b: bool(0.8 <= bb[b]["gain"]["F"][vn]["wH"] <= 1.2 and 0.8 <= bb[b]["gain"]["F"][vn]["wD"] <= 1.2 and bb[b]["gain"]["F"][vn]["okH"] >= 0.8 and bb[b]["gain"]["F"][vn]["okD"] >= 0.8) for b in BN}
        return c2, c3
    c2, c3 = c23("V1")
    c2b, c3b = c23("V2")
    c4 = (summ["33"]["mean"] >= 0.10) == c1
    V[run].update(summary=summ, c1=bool(c1), c1_nsig=float((summ["F"]["mean"] - 0.10) / summ["F"]["err"]), c2_bias_V1=c2, c3_pull_V1=c3, c2_bias_V2=c2b, c3_pull_V2=c3b, c4_same_sign=bool(c4),
                  status="NIEROZSTRZYGNIETY" if not c4 else ("POTWIERDZONY" if (c1 and all(c2.values()) and all(c3.values())) else "NIEPOTWIERDZONY"))
st = [V[r]["status"] for r in RUNS]
V["verdict"] = "NIEROZSTRZYGNIETY" if "NIEROZSTRZYGNIETY" in st else ("POTWIERDZONY" if all(s == "POTWIERDZONY" for s in st) else "NIEPOTWIERDZONY")
json.dump(V, open(os.path.join(K1C, "verdict.json"), "w"), indent=1)

p(f"# K1c — wersja na kandydatach ({S['nrep']} sparowanych replik zdarzeniowych, {S['ntoy']} pseudo-eksperymentów)\n")
p("Wariant kształtu sygnału: " + ", ".join(f"{k}: **{v}**" for k, v in S["fitvar"].items()) + "; S0: " + ", ".join(f"{k}: {v:.0f}" for k, v in S["S0"].items()) + "\n")
p("## 1. Werdykt (WP70)\n")
for run in RUNS:
    v = V[run]
    s = v["summary"]
    p(f"- **{run}: {v['status']}** — średni zysk (def. K1): " + ", ".join(f"{VN[x]} {pc(s[x]['mean'], s[x]['err'])} %" for x in VAR) + f"; kryt. 1: {'PASS' if v['c1'] else 'FAIL'} ({v['c1_nsig']:+.1f}σ od progu); "
      f"kryt. 2 (V1): {v['c2_bias_V1']}; kryt. 3 (V1): {v['c3_pull_V1']}; kryt. 4 (pol3/pol3 ten sam znak): {'tak' if v['c4_same_sign'] else 'NIE'}; opisowo V2: kryt. 2 {v['c2_bias_V2']}, kryt. 3 {v['c3_pull_V2']}; "
      f"średni zysk skorygowany F-test: V1 {100*s['F']['mean_corr_V1']:+.1f} %, V2 {100*s['F']['mean_corr_V2']:+.1f} %; pol3/pol3: V1 {100*s['33']['mean_corr_V1']:+.1f} %, V2 {100*s['33']['mean_corr_V2']:+.1f} %; pol4/pol4: V1 {100*s['44']['mean_corr_V1']:+.1f} %")
p(f"\n**Łącznie: zysk {V['verdict']}.**\n")
p("## 2. Zysk przy WP70, definicja K1: 1 − [σ/N]_DLZ/[σ/N]_H [%]; błąd: sparowany bootstrap po zdarzeniach\n")
p("| run | pT | rzędy | rząd H/DLZ | zysk | N_DLZ/N_H | b H / DLZ (V1) | pull w H / DLZ (V1) | średni pull H / DLZ | zysk skor. V1 | b H / DLZ (V2) | w H / DLZ (V2) | zysk skor. V2 | zysk z rozrzutu toyów V1 | σ_boot/σ_hesse H / DLZ |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for run, bn, var, g in T_gain:
    a, b = g["V1"], g["V2"]
    p(f"| {run} | {bn} | {VN[var]} | {g['oH']}/{g['oD']} | {pc(g['g'], g['g_err'])} | {pm(g['Nratio'], g['Nratio_err'])} | {a['bH']:.3f}±{a['bH_err']:.3f} / {a['bD']:.3f}±{a['bD_err']:.3f} | {a['wH']:.2f} / {a['wD']:.2f} | {a['pmH']:+.2f} / {a['pmD']:+.2f} | "
      f"{100*a['g_corr']:+.1f} | {b['bH']:.3f} / {b['bD']:.3f} | {b['wH']:.2f} / {b['wD']:.2f} | {100*b['g_corr']:+.1f} | {100*a['g_toy']:+.1f} | {g['sigboot_over_hesse_H']:.2f} / {g['sigboot_over_hesse_D']:.2f} |")
p("\n## 3. Stosunek R = (N/ε)_DLZ/(N/ε)_H per WP: surowy → skorygowany (V1; błąd: bootstrap ⊕ stat. toyów ⊕ |V1 − V2|)\n")
p("| run | pT | rzędy | WP50 | WP60 | WP70 | WP80 | WP90 | ⟨R⟩ | ⟨R_corr⟩ | χ²(R_corr = 1)/5 |")
p("|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(T_R)
p("\n## 4. Skan N_fit/ε^MC (WP 50/60/70/80/90): wartości, χ²_flat/4 surowe → po korekcie V1, nachylenie [1/ε]\n")
p("| run | pT | rzędy | H | H+DLZ | rzędy per WP (H ; DLZ) | repliki ok |")
p("|---|---|---|---|---|---|---|")
L.extend(T_scan)
p("\n## 5. Fity danych i F-test (wariant kształtu komórki): N ± σ (σ/N); χ²/ndf; `*` granica, `!` nieważny\n")
p("| run | pT | klas. | WP | pol2 | pol3 | pol4 | p(2→3) | p(3→4) | p(2→4) | wybrany |")
p("|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(T_ft)
p("\n## 6. Pseudo-dane: b = ⟨N_fit⟩/N_true ± błąd (As = Asimov; w = szerokość pulla; średni pull; udział udanych); V2 poza WP70 tylko Asimov\n")
p("| run | pT | klas. | WP | wybrany | V1 pol2 | V1 pol3 | V1 pol4 | V2 pol2 | V2 pol3 | V2 pol4 |")
p("|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(T_bias)
p("\n## 7. Tło prawdy: n_all(m)·ε_B(m)\n")
p("| run | pT | klas. | WP | N_true | B_win V1 | B_win V2 | S/B prawdy V1 | χ²/bin SB (dane − ogony − B) V1 | V2 | dane − (S+B) w oknie ±45 MeV, V1 | V2 | udział sygnału MC w oknie fitu |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(T_bkg)
open(os.path.join(K1C, "tables.md"), "w").write("\n".join(L) + "\n")
n3 = 7 + 2 * 1 + 4 + len(T_gain) + 5 + len(T_R)
print("\n".join(L[:n3]))
