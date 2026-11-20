#!/usr/bin/env python3
"""SD_var/predict.py — PRZEWIDYWANIA liczby tagowanych D (domknięcie 19/20.11.2026, projekt B). Nic tu nie jest fitem ani liczbą z okna sygnału regionu tagowanego.
Wejścia (tylko odczyt): out/results.json (R nominalne per bin pT i obwiednia, analyze_var.py), FITBIAS_2026-11-19/out/{var_data,var_data_k1,data_summary}.json
(uzyski INKLUZYWNE: D0 model M8, D+ model M0), TAGREG_2026-11-18/out/s1.json (ułamki kandydatów pasm bocznych w T2 — liczby kandydatów, bez mas w oknie), LUMI (stałe niżej).
Wyników TAGGED_2026-11-19 (out/, tabele/, RAPORT) skrypt NIE czyta.
Użycie: .venv_h5/bin/python SD_var_2026-11-19/predict.py [--res out/results.json] [--out out]  →  <out>/pred.json, <out>/pred_tables.md
"""
import argparse, json, math, os
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.dirname(HERE)
BINS = ["3-5", "5-8", "8-30"]; SPN = {"D0": "D⁰", "Dp": "D⁺", "Ds": "D*⁺"}; FB = {"D0": "D0", "Dp": "Dplus", "Ds": "Dstar"}
# --- stałe (źródła w RAPORT_domkniecie.md §3)
AE, AE_LO, AE_HI = 0.55, 0.35, 0.80          # akceptancja × wydajność AFP (założenie TAGREG §3; niezmierzone)
AE_NEW, AE_NEW_ERR, EXI_NEW = 0.92, 0.05, 1.000   # 20.11: A·ε_AFP zmierzone (AFP_akceptancja_2026-11-20/PLAN.md §0, §5: A 0,97 ± 0,03 × ε_AFP 0,95 ± 0,02, w tym „dokładnie 1 proton”); ε_ξcal 427929 z pomiaru
FNEW = AE_NEW * EXI_NEW / (0.55 * 0.99)            # nowy/stary czynnik detektorowy = 1,690 (A·ε ×1,673; ε_ξcal ×1,010)
EXI = 0.99                                    # wydajność cięcia ξ_cal przyjęta w nocie dla scenariusza M (pomiar: 427929 1,000 +0/−0,005; 428770 0,96–1,00)
FACC = {"427929": 0.0087, "R1": 0.038}       # ułamek tagów przypadkowych (AFP_mix v2; R1: klasa noalg_L1AFP z TAGREG §4.2, zakres 0,038–0,086)
S_LO, S_HI = 0.5, 2.0                         # zakres założenia scenariusza S (TAGREG §3)
L = {"427929": 0.4648, "R1": 0.4687, "R2": 27.65}
# przeniesienie 427929 → R1 (428770): L × ε(L1 AFP | tag offline, zdarzenia z D0: 7139/10438) × weto npu (0,9408/0,9848, nota §6.7) × ξ_cal (0,974/1,000, TAGREG §4.2 klasa noalg)
K_R1 = dict(lumi=L["R1"] / L["427929"], l1afp=7139 / 10438, veto=0.9408 / 0.9848, xical=0.974 / 1.000, afp=0.94 / 0.95)   # afp: ε_AFP 428770/427929 (AFP_akceptancja §3)
K_R1_ALT = 8267 / 10438 * (1 - 0.038) / (1 - 0.0087)                       # alternatywa: stosunek liczby zdarzeń T2 z kandydatem D0 (R1 428770 / 427929) × (1 − f_acc) — górna granica zakresu
K_R1_VLO = 0.8824 / 0.9848                    # weto npu jak w klasie „inne trygery” (88,24 %) — dolna granica zakresu
# uzyski inkluzywne D* (tag D*, WP80; TAGREG §3 za MC_templates WYNIKI §2a) i dawne (nota, fit C0, WP80)
DSTAR = [(513, 38), (265, 18), (107, 11)]
OLD = {"D0": [6773, 2010, 813], "Dp": [4999, 1833, 767], "Ds": [513, 265, 107]}
OLD_NOTE = {"D0": dict(S=(35, 18, 71), M=(18.6, 7.4, 74), M18=(39, 7, 152)), "Dp": dict(S=(18, 9, 36), M=(31, 5, 120)), "Ds": dict(S=(2, 1, 4), M=(4, 0.6, 14))}


def incl(run="427929"):
    # 20.11 (CVER_2026-11-20/RAPORT.md pkt 5): „BR w MC ≠ PDG” WYCOFANE → κ_dane = κ_MC (r = 1). Nominał = M8 BEZ skali κ (var_data_k1.json, kscale 1,0);
    # wariant ze skalą κ × 1,364 (var_data.json, dawny nominał FITBIAS) tylko do pokazania różnicy. Żadnej „poprawki BR” w przewidywaniach.
    d = json.load(open(os.path.join(W, "FITBIAS_2026-11-19/out/var_data_k1.json")))["res"]; k1 = json.load(open(os.path.join(W, "FITBIAS_2026-11-19/out/var_data.json")))["res"]
    ds = json.load(open(os.path.join(W, "FITBIAS_2026-11-19/out/data_summary.json")))["wp70"]
    out = {}
    for wp in ("WP70", "WP80"):
        for sp, m in (("D0", "M8"), ("Dp", "M0")):   # M8 tu = κ_MC
            rows = []; rows1 = []
            for b in BINS:
                r = next(x for x in d if x["sp"] == FB[sp] and x["run"] == run and x["bin"] == b)["S"][wp][m]
                r1 = next(x for x in k1 if x["sp"] == FB[sp] and x["run"] == run and x["bin"] == b)["S"][wp][m]
                e = r["err"] * 1.3                                       # HESSE × 1,3 (FITBIAS: rms pull zabawek)
                if wp == "WP70": e = max(e, ds[f"{FB[sp]}|{run}|{b}"]["s1"])   # albo bootstrap po zdarzeniach, jeśli większy
                rows.append((r["N"], e)); rows1.append(r1["N"])
            out[(wp, sp)] = dict(N=rows, N_k1=rows1, model=m + (", κ_MC" if m == "M8" else ""))
        out[(wp, "Ds")] = dict(N=[(float(a), float(b)) for a, b in DSTAR], N_k1=[float(a) for a, _ in DSTAR], model="tag D* (Δm), WP80 — bez zmian")
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--res", default=os.path.join(HERE, "out/results.json")); ap.add_argument("--out", default=os.path.join(HERE, "out")); a = ap.parse_args()
    R = json.load(open(a.res)); U = R["uncertainty"]; nomname = U["D0"]["nominal"]; NOM = R["scenarios"][nomname]
    S1 = json.load(open(os.path.join(W, "TAGREG_2026-11-18/out/s1.json")))["427929"]["tot"]
    INC = incl(); kR1 = math.prod(K_R1.values()); res = {"K_R1": dict(K_R1, total=kR1, alt_hi=K_R1_ALT, lo=kR1 * K_R1_VLO / K_R1["veto"]), "rows": []}
    T = ["| WP | gatunek | N_incl 3–5 / 5–8 / 8–30 (model) | suma | R nominalne w binach [%] | ⟨R⟩ ważone N [%] | P_M = ⟨R⟩·A·ε·ε_ξ [%] (obwiednia modeli) | **M: N_tag** (obwiednia modeli) [z A·ε 0,35–0,80] | **M po pomiarze AFP: P [%]; N_tag** (A·ε = 0,92, ε_ξ = 1,00) | P_SB [%] | **S: N_tag** (×0,5–2) | przypadkowe |", "|" + "---|" * 12]
    for wp in ("WP70", "WP80"):
        for sp in ("D0", "Dp", "Ds"):
            I = INC[(wp, sp)]; N = [x[0] for x in I["N"]]; Ns = sum(N); eN = math.sqrt(sum(x[1] ** 2 for x in I["N"]))
            rpt = NOM[sp]["R_pt"][1:]; u = U[sp]; fd, fu = u["factor_dn"], u["factor_up"]
            M = sum(n * r for n, r in zip(N, rpt)) * AE * EXI; PM = M / Ns
            # A·ε w kwadraturze (ln) z obwiednią modeli — tylko dodatkowa kolumna
            fd2 = math.exp(-math.hypot(math.log(fd), math.log(AE_LO / AE))); fu2 = math.exp(math.hypot(math.log(fu), math.log(AE_HI / AE)))
            key = FB[sp]; psb = [S1["T2"]["sb_cands"][f"{key}|{b}"] / S1["GRL"]["sb_cands"][f"{key}|{b}"] for b in BINS]
            fa = FACC["427929"]; S = sum(n * p for n, p in zip(N, psb)) * (1 - fa); A = sum(n * p for n, p in zip(N, psb)) * fa
            Mk1 = sum(n * r for n, r in zip(I["N_k1"], rpt)) * AE * EXI
            row = dict(wp=wp, sp=sp, model=I["model"], N=N, N_sum=Ns, N_err=eN, R_pt=rpt, R_mean=M / Ns / AE / EXI, P_M=PM, P_M_lo=PM * fd, P_M_hi=PM * fu, M=M, M_lo=M * fd, M_hi=M * fu, M_lo_AE=M * fd2, M_hi_AE=M * fu2,
                       M_new=M * FNEW, M_new_lo=M * FNEW * fd, M_new_hi=M * FNEW * fu, P_M_new=PM * FNEW, P_M_new_lo=PM * FNEW * fd, P_M_new_hi=PM * FNEW * fu, f_dn=fd, f_up=fu, P_SB=psb, P_S=S / Ns, S=S, S_lo=S * S_LO, S_hi=S * S_HI, acc=A, BR_effect_M=M / Mk1 - 1 if Mk1 else None,
                       R1=dict(M_new=M * kR1 * FNEW, M_new_lo=M * kR1 * FNEW * fd, M_new_hi=M * kR1 * FNEW * fu, acc_M_new=M * kR1 * FNEW * FACC["R1"] / (1 - FACC["R1"]), M=M * kR1, M_lo=M * kR1 * fd, M_hi=M * kR1 * fu, S=S * kR1, S_lo=S * kR1 * S_LO, S_hi=S * kR1 * S_HI, acc_M=M * kR1 * FACC["R1"] / (1 - FACC["R1"]), acc_S=S * kR1 * FACC["R1"] / (1 - FACC["R1"])))
            if wp == "WP80":   # ten sam rachunek na dawnych uzyskach (C0) — rozdzielenie zmiany R od zmiany N_incl
                No = OLD[sp]; row["old_basis"] = dict(N=No, M_newR=sum(n * r for n, r in zip(No, rpt)) * AE * EXI, M_oldR=sum(No) * 0.0075 * AE * EXI, S=sum(n * p for n, p in zip(No, psb)) * (1 - fa))
            res["rows"].append(row)
            T.append(f"| {wp} | {SPN[sp]} | {' / '.join(f'{n:.0f}' for n in N)} ({I['model']}) | {Ns:.0f} ± {eN:.0f} | {' / '.join(f'{100*r:.3f}' for r in rpt)} | {100*row['R_mean']:.3f} | {100*PM:.3f} ({100*PM*fd:.3f}–{100*PM*fu:.2f}) | "
                     f"**{M:.1f}** ({M*fd:.1f}–{M*fu:.0f}) [{M*fd2:.1f}–{M*fu2:.0f}] | {100*PM*FNEW:.3f} ({100*PM*FNEW*fd:.3f}–{100*PM*FNEW*fu:.2f}); **{M*FNEW:.1f}** ({M*FNEW*fd:.1f}–{M*FNEW*fu:.0f}) | {' / '.join(f'{100*p:.3f}' for p in psb)} | **{S:.1f}** ({S*S_LO:.1f}–{S*S_HI:.0f}) | {A:.2f} |")
    T2 = ["| WP | gatunek | 427929 M (stary czynnik 0,55 × 0,99) | 427929 M (nowy 0,92 × 1,00) | 427929 S | R1 M (stary) | R1 M (nowy) | R1 S | R1 przypadkowe (M stary / M nowy / S) |", "|" + "---|" * 9]
    for r in res["rows"]:
        q = r["R1"]; T2.append(f"| {r['wp']} | {SPN[r['sp']]} | {r['M']:.1f} ({r['M_lo']:.1f}–{r['M_hi']:.0f}) | {r['M_new']:.1f} ({r['M_new_lo']:.1f}–{r['M_new_hi']:.0f}) | {r['S']:.1f} ({r['S_lo']:.1f}–{r['S_hi']:.0f}) | {q['M']:.1f} ({q['M_lo']:.1f}–{q['M_hi']:.0f}) | {q['M_new']:.1f} ({q['M_new_lo']:.1f}–{q['M_new_hi']:.0f}) | {q['S']:.1f} ({q['S_lo']:.1f}–{q['S_hi']:.0f}) | {q['acc_M']:.2f} / {q['acc_M_new']:.2f} / {q['acc_S']:.2f} |")
    T3 = ["| gatunek | nota: S | nota: M | dawne N_incl (C0, WP80) z nowym R: M | nowe WP80: S | nowe WP80: M (0,55) | nowe WP80: M (0,92) | nowe WP70: S | nowe WP70: M (0,55) | nowe WP70: M (0,92) | M z κ_MC (nominał) względem M z wycofaną skalą κ × 1,364: WP70 / WP80 |", "|" + "---|" * 11]
    for sp in ("D0", "Dp", "Ds"):
        r8 = next(x for x in res["rows"] if x["wp"] == "WP80" and x["sp"] == sp); r7 = next(x for x in res["rows"] if x["wp"] == "WP70" and x["sp"] == sp); o = OLD_NOTE[sp]
        f = lambda t: f"{t[0]:g} ({t[1]:g}–{t[2]:g})"; g = lambda r, k: f"{r[k]:.1f} ({r[k+'_lo']:.1f}–{r[k+'_hi']:.0f})"
        T3.append(f"| {SPN[sp]} | {f(o['S'])} | {f(o['M'])} | {r8['old_basis']['M_newR']:.1f} | {g(r8, 'S')} | {g(r8, 'M')} | {g(r8, 'M_new')} | {g(r7, 'S')} | {g(r7, 'M')} | {g(r7, 'M_new')} | {100*r7['BR_effect_M']:+.2f} % / {100*r8['BR_effect_M']:+.2f} % |")
    hd = (f"# Przewidywania liczby tagowanych D — tabele (predict.py; R z `{os.path.relpath(a.res, HERE)}`, nominał {nomname})\n\n"
          f"Stałe: czynnik detektorowy stary A·ε_AFP = {AE} ({AE_LO}–{AE_HI}; założenie) × ε_ξcal = {EXI}; nowy A·ε_AFP = {AE_NEW} ± {AE_NEW_ERR} (pomiar AFP_akceptancja 20.11) × ε_ξcal = {EXI_NEW:.2f}; nowy/stary = {FNEW:.3f}; f_acc = {FACC['427929']} (427929), {FACC['R1']} (R1); k(427929 → R1) = "
          + " × ".join(f"{v:.4f}" for v in K_R1.values()) + f" = {kR1:.3f} (zakres {res['K_R1']['lo']:.2f}–{K_R1_ALT:.2f}).\n")
    open(os.path.join(a.out, "pred_tables.md"), "w", encoding="utf-8").write(hd + "\n## 1. Run 427929 (region MBTS, L = 0,4648 nb⁻¹, pT 3–30 GeV): składniki\n\n" + "\n".join(T) + "\n\n## 2. 427929 i R1 runu 428770 (L = 0,4687 nb⁻¹)\n\n" + "\n".join(T2) + "\n\n## 3. Stare → nowe\n\n" + "\n".join(T3) + "\n")
    json.dump(res, open(os.path.join(a.out, "pred.json"), "w"), indent=1, ensure_ascii=False, default=float)
    print(open(os.path.join(a.out, "pred_tables.md"), encoding="utf-8").read())


if __name__ == "__main__":
    main()
