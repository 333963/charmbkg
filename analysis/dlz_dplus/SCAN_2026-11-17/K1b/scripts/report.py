"""K1b: tabele i werdykt z out/k1b.json → K1b/tables.md, K1b/verdict.json. Użycie: python report.py"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K1B = os.path.dirname(HERE)
SCAN = os.path.dirname(K1B)
DATA = os.path.join(os.path.dirname(SCAN), "DATA_2026-11-17")
S = json.load(open(os.path.join(K1B, "out", "k1b.json")))
RUNS, BINS, TG, WPS = ["427929", "428770"], ["3-5", "5-8", "8-30"], {"H": "H_s12345", "DLZ": "HDLZ_s12345"}, [0.5, 0.7, 0.8]
L = []
p = L.append
pc = lambda v, e: f"{100*v:+.1f} ± {100*e:.1f}"  # noqa: E731
V = {}

p(f"# K1b — zysk σ(N)/N H+DLZ vs H z tłem exp(pol2/3/4) ({S['nrep']} replik binowych, {S['ntoy']} pseudo-eksperymentów)\n")
p("Wariant kształtu sygnału per komórka: " + ", ".join(f"{k}: **{v}**" for k, v in S["fitvar"].items()) + "\n")

rows_gain, rows_ft, rows_toy, rows_rep = [], [], [], []
for run in RUNS:
    res = json.load(open(os.path.join(DATA, "runs", "lxplus", "ana", f"res_{run}.json")))
    rep = np.load(os.path.join(DATA, "runs", "lxplus", "ana", f"rep_{run}.npz"))
    V[run] = {"bins": {}}
    for b in BINS:
        fv = S["fitvar"][f"{run}|{b}"]
        for wp in WPS:
            nom = {t: S["nominal"][f"{run}|{b}|{t}|{wp}"][fv] for t in TG}
            cel = {t: S["cells"][f"{run}|{b}|{t}|{wp}"] for t in TG}
            sel = {t: nom[t]["sel"] for t in TG}
            for t in TG:
                f, ft = nom[t]["fits"], nom[t]["ftest"]
                oth = S["nominal"][f"{run}|{b}|{t}|{wp}"]["fix" if fv == "free" else "free"]
                rows_ft.append(f"| {run} | {b} | {t} | {int(wp*100)} | " + " | ".join(
                    f"{f[o]['N']:.0f} ± {f[o]['N_err']:.0f} ({100*f[o]['rel']:.1f} %); {f[o]['chi2_ndf']:.2f}{'*' if f[o]['at_limit'] else ''}{'' if f[o]['valid'] else '!'}" for o in "234")
                    + f" | {ft['2->3']['p']:.3g} ({ft['2->3']['p_lr']:.3g}) | {ft['3->4']['p']:.3g} ({ft['3->4']['p_lr']:.3g}) | {ft['2->4']['p']:.3g} | **pol{sel[t]}** | pol{oth['sel']}: "
                    + " / ".join(f"{oth['fits'][o]['N']:.0f} ± {oth['fits'][o]['N_err']:.0f}" for o in "234") + " |")
                ty = cel[t]["toy"]
                tr = cel[t]["truth"]
                rows_toy.append(f"| {run} | {b} | {t} | {int(wp*100)} | {tr['A_chi2sb']:.2f}{'' if tr['A_valid'] else '!'} | " + " | ".join(
                    f"{ty[o]['b']:.3f} ± {ty[o]['b_err']:.3f} (C {ty[o]['asimov_C']:.3f}); pull {ty[o]['pull_mean']:+.2f}, w {ty[o]['pull_w']:.2f} ± {ty[o]['pull_w_err']:.2f}; ok {ty[o]['f_ok']:.2f}" for o in "234") + f" | pol{sel[t]} |")
            # (a) odtworzenie K1 z zapisanych fitów pol2 + repliki sparowane
            F2 = res["fits"] if fv == "free" else res["fits_fixshape"]
            pre = "" if fv == "free" else "fix|"
            kH, kD = f"{TG['H']}|{b}|{wp}", f"{TG['DLZ']}|{b}|{wp}"
            rH, rD = rep[pre + kH], rep[pre + kD]
            ok = (rH[:, 0] > 0) & (rD[:, 0] > 0) & (rH[:, 4] > 0) & (rD[:, 4] > 0)
            relH, relD = F2[kH]["N_err"] / F2[kH]["N"], F2[kD]["N_err"] / F2[kD]["N"]
            gK1, gK1e = 1 - relD / relH, float(np.std(1 - (rD[ok, 1] / rD[ok, 0]) / (rH[ok, 1] / rH[ok, 0])))
            gA, gAe = relH / relD - 1, float(np.std((rH[ok, 1] / rH[ok, 0]) / (rD[ok, 1] / rD[ok, 0]) - 1))
            nrA, nrAe = F2[kD]["N"] / F2[kH]["N"], float(np.std(rD[ok, 0] / rH[ok, 0]))

            def variant(oH, oD):
                fH, fD = nom["H"]["fits"][str(oH)], nom["DLZ"]["fits"][str(oD)]
                bH, bD = np.array(cel["H"]["boot"][str(oH)]), np.array(cel["DLZ"]["boot"][str(oD)])
                o2 = (bH[:, 0] > 0) & (bD[:, 0] > 0) & (bH[:, 4] > 0) & (bD[:, 4] > 0)
                g = fH["rel"] / fD["rel"] - 1
                ge = float(np.std((bH[o2, 1] / bH[o2, 0]) / (bD[o2, 1] / bD[o2, 0]) - 1)) if o2.sum() > 1 else float("nan")
                tH, tD = cel["H"]["toy"][str(oH)], cel["DLZ"]["toy"][str(oD)]
                gc = (1 + g) * tH["b"] / tD["b"] - 1
                gcp = (1 + g) * tH["b"] / tD["b"] * tH["pull_w"] / tD["pull_w"] - 1
                # zysk empiryczny: rozrzut N/N_true w pseudo-danych (nie zależy od hesjanu)
                gt = (tH["rms"] / tH["b"]) / (tD["rms"] / tD["b"]) - 1
                eps = {t: res["shapes"][f"{TG[t]}|{b}|{wp}"]["eps_mc"] for t in TG}
                nr = (fD["N"] / eps["DLZ"]) / (fH["N"] / eps["H"])
                nre = float(np.std(bD[o2, 0] / bH[o2, 0]) / eps["DLZ"] * eps["H"]) if o2.sum() > 1 else float("nan")
                return dict(oH=oH, oD=oD, g=g, g_err=ge, n_ok=int(o2.sum()), g_corr=gc, g_corr_pull=gcp, g_toy=gt, bH=tH["b"], bD=tD["b"], wH=tH["pull_w"], wD=tD["pull_w"],
                            okH=tH["f_ok"], okD=tD["f_ok"], Nratio=nr, Nratio_err=nre, Nratio_corr=nr * tH["b"] / tD["b"],
                            relH=fH["rel"], relD=fD["rel"], NH=fH["N"], ND=fD["N"], validH=fH["valid"], validD=fD["valid"],
                            sigboot_over_hesse_H=float(np.std(bH[bH[:, 4] > 0, 0]) / fH["N_err"]), sigboot_over_hesse_D=float(np.std(bD[bD[:, 4] > 0, 0]) / fD["N_err"]))
            va, vb, vc = variant(2, 2), variant(sel["H"], sel["DLZ"]), variant(4, 4)
            rows_rep.append(f"| {run} | {b} | {int(wp*100)} | {fv} | {pc(gK1, gK1e)} | {pc(gA, gAe)} | {pc(va['g'], va['g_err'])} | {nrA:.3f} ± {nrAe:.3f} | {va['Nratio_err']:.3f} | {int(ok.sum())} |")
            rows_gain.append(f"| {run} | {b} | {int(wp*100)} | {pc(gA, gAe)} | pol{sel['H']}/pol{sel['DLZ']} | {pc(vb['g'], vb['g_err'])} | {pc(vc['g'], vc['g_err'])} | "
                             f"{vb['bH']:.3f} / {vb['bD']:.3f} | {vb['wH']:.2f} / {vb['wD']:.2f} | {100*vb['g_corr']:+.1f} | {100*vb['g_corr_pull']:+.1f} | {100*vb['g_toy']:+.1f} | "
                             f"{100*va['g_corr']:+.1f} | {100*vc['g_corr']:+.1f} ({vc['bH']:.3f}/{vc['bD']:.3f}; w {vc['wH']:.2f}/{vc['wD']:.2f}) | {vb['Nratio']:.3f} ± {vb['Nratio_err']:.3f} → {vb['Nratio_corr']:.3f} | "
                             f"{vb['sigboot_over_hesse_H']:.2f} / {vb['sigboot_over_hesse_D']:.2f} |")
            if wp == 0.7:
                V[run]["bins"][b] = dict(fitvar=fv, a=dict(g=gA, g_err=gAe, g_K1def=gK1, g_K1def_err=gK1e, Nratio=nrA, Nratio_err=nrAe, g_corr=va["g_corr"]), b=vb, c=vc)
    bb = V[run]["bins"]
    mean = {x: float(np.mean([bb[b][x]["g"] for b in BINS])) for x in "abc"}
    mean_err = {x: float(np.sqrt(np.sum([bb[b][x]["g_err"] ** 2 for b in BINS])) / 3) for x in "abc"}
    mean_corr = {x: float(np.mean([bb[b][x]["g_corr"] for b in BINS])) for x in "abc"}
    c1 = mean["b"] >= 0.10
    c2 = {b: bool(abs(bb[b]["b"]["bH"] - 1) < 0.05 and abs(bb[b]["b"]["bD"] - 1) < 0.05 and bb[b]["b"]["okH"] >= 0.8 and bb[b]["b"]["okD"] >= 0.8 and bb[b]["b"]["validH"] and bb[b]["b"]["validD"]) for b in BINS}
    c3 = {b: bool(0.8 <= bb[b]["b"]["wH"] <= 1.2 and 0.8 <= bb[b]["b"]["wD"] <= 1.2 and bb[b]["b"]["okH"] >= 0.8 and bb[b]["b"]["okD"] >= 0.8) for b in BINS}
    V[run].update(mean_gain=mean, mean_gain_err=mean_err, mean_gain_corr=mean_corr, c1_mean_b_ge_10=bool(c1), c2_bias=c2, c3_pull=c3,
                  PASS=bool(c1 and all(c2.values()) and all(c3.values())), mean_corr_b_ge_10=bool(mean_corr["b"] >= 0.10))
V["verdict"] = "POTWIERDZONY" if all(V[r]["PASS"] for r in RUNS) else "NIEPOTWIERDZONY"
json.dump(V, open(os.path.join(K1B, "verdict.json"), "w"), indent=1)

p("## 1. Werdykt (WP70)\n")
for run in RUNS:
    v = V[run]
    p(f"- **{run}: {'PASS' if v['PASS'] else 'FAIL'}** — średni zysk (a) {100*v['mean_gain']['a']:+.1f} ± {100*v['mean_gain_err']['a']:.1f} %, (b) {100*v['mean_gain']['b']:+.1f} ± {100*v['mean_gain_err']['b']:.1f} %, "
      f"(c) {100*v['mean_gain']['c']:+.1f} ± {100*v['mean_gain_err']['c']:.1f} %; skorygowany o obciążenie: (a) {100*v['mean_gain_corr']['a']:+.1f}, (b) {100*v['mean_gain_corr']['b']:+.1f}, (c) {100*v['mean_gain_corr']['c']:+.1f} %; "
      f"kryt. 1 (b ≥ 10 %): {'PASS' if v['c1_mean_b_ge_10'] else 'FAIL'}; kryt. 2 (obciążenie): {v['c2_bias']}; kryt. 3 (pull): {v['c3_pull']}")
p(f"\n**Łącznie: zysk {V['verdict']}.**\n")
p("## 2. Zysk g = [σ/N]_H/[σ/N]_DLZ − 1 [%]\n")
p("| run | pT | WP | (a) pol2/pol2 (repl. sparowane) | rząd F-test H/DLZ | (b) wybrany | (c) pol4/pol4 | b H / DLZ (b) | pull w H / DLZ (b) | (b) skor. o b | (b) skor. o b i pull | (b) z rozrzutu toyów | (a) skor. o b | (c) skor. o b (b; w) | N_DLZ/N_H (b) → skor. | σ_boot/σ_hesse H / DLZ |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(rows_gain)
p("\nBłędy (b), (c): bootstrap binowy niezależny dla H i H+DLZ (zawyżone, KRYTERIA_K1b §1). „z rozrzutu toyów” = [rms(N)/⟨N⟩]_H/[rms(N)/⟨N⟩]_DLZ − 1 z pseudo-danych (bez hesjanu).\n")
p("## 3. Odtworzenie K1 (pol2/pol2, zapisane fity + sparowane repliki zdarzeniowe)\n")
p("| run | pT | WP | kształt | zysk def. K1: 1 − rel_DLZ/rel_H [%] | zysk def. zlecenia [%] | to samo, błąd z bootstrapu binowego niezależnego | N_DLZ/N_H ± (sparowany) | ± (binowy niezależny) | repliki ok |")
p("|---|---|---|---|---|---|---|---|---|---|")
L.extend(rows_rep)
p("\n## 4. Fity danych i F-test (wariant kształtu komórki)\n")
p("Komórka: N ± σ (σ/N); χ²/ndf; `*` na granicy, `!` MIGRAD nieważny. p: F-test (w nawiasie p z ilorazu wiarygodności).\n")
p("| run | pT | klas. | WP | pol2 | pol3 | pol4 | p(2→3) | p(3→4) | p(2→4) | wybrany | drugi wariant kształtu: wybór; N pol2/3/4 |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|")
L.extend(rows_ft)
p("\n## 5. Pseudo-dane (prawda: sygnał MC + exp(pol6) z SB danych, ogony sygnału odjęte), fit rzędu 2/3/4\n")
p("Komórka: b = ⟨N_fit/N_true⟩ ± błąd (Asimov dla tła C); średni pull, szerokość pulla; udział udanych.\n")
p("| run | pT | klas. | WP | χ²/ndf SB prawdy A | fit pol2 | fit pol3 | fit pol4 | wybrany |")
p("|---|---|---|---|---|---|---|---|---|")
L.extend(rows_toy)
open(os.path.join(K1B, "tables.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L[:14 + len(rows_gain) + 4 + len(rows_rep) + 4]))
