#!/usr/bin/env python3
"""AUDYT_noty 20.11.2026 (projekt B) — przeliczenia arytmetyczne do pkt 1–7 zlecenia. ŻADNYCH fitów; tylko odczyt zapisanych liczb.
Wejścia (tylko odczyt): CVER/out/cver.json, XSEC_final/out/final.json, TAGXSEC/out/tagxsec_excl.json, TAGGED/out/prep.json (tylko uzyski INKLUZYWNE jako wagi),
SD_var/out/results.json. Liczby przepisane z raportów mają w komentarzu źródło.
Użycie: ~/Documents/charm/.venv_h5/bin/python przeliczenia.py  → ../out/przeliczenia.json, ../out/przeliczenia.md"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE); W = os.path.dirname(TOP)
J = lambda p: json.load(open(os.path.join(W, p)))  # noqa: E731
OUT, MD = {}, []
BN = ["3-5", "5-8", "8-30"]
q = lambda *a: math.sqrt(sum(x * x for x in a))  # noqa: E731


def p_chi2(x, ndf):
    if ndf == 1: return math.erfc(math.sqrt(x / 2))
    if ndf == 2: return math.exp(-x / 2)
    if ndf == 3: return math.erfc(math.sqrt(x / 2)) + math.sqrt(2 * x / math.pi) * math.exp(-x / 2)
    raise ValueError


def p_two_sided(z): return math.erfc(abs(z) / math.sqrt(2))


def h(t): MD.append("\n## " + t + "\n")


# ============================================================ 1. BR ekskluzywne MC / PDG
h("1. BR ekskluzywne MC/PDG")
cv = J("CVER_2026-11-20/out/cver.json")["sp"]
PDG = {"D0": (0.03936, 0.00030, 0.000148), "Dplus": (0.0931, 0.0014, 0.0)}      # PDG 2026 (BR_frag_2026-11-18 §1): tryb, błąd, DCS
BKG = {"D0": (3.974, 0.014, 3.985, 3.955), "Dplus": (9.404, 0.030, 9.401, 9.336)}  # BKGCAT tabele/br.md: |η|<2,5 wszystkie pT; ± stat; bez cięcia η; pT ≥ 3 ∧ |η| < 2,2
O = OUT["1_BR"] = {}
MD.append("| gat. | N_ekskl (z FSR) / N_wszystkie, pT ≥ 3, \\|η\\| < 2,2 (CVER) | BR_MC fid. [%] | f_ekskl w decayType [%] | PDG bez DCS [%] | **MC/PDG (fid., bez DCS)** | CVER: /PDG z DCS | BKGCAT: \\|η\\| < 2,5, wszystkie pT [%] → MC/PDG | BKGCAT kolumna kontrolna fid. [%] | (BKGCAT − fid.)/σ |")
MD.append("|---|---|---|---|---|---|---|---|---|---|")
for sp in ("D0", "Dplus"):
    b = cv[sp]["bins"]; ne = sum(b[k]["N_gen_excl"] for k in BN); na = sum(b[k]["N_gen_all"] for k in BN); nd = sum(b[k]["N_gen_decayType"] for k in BN)
    br = ne / na; e = math.sqrt(br * (1 - br) / na); pdg, dp, dcs = PDG[sp]
    r = br / pdg; er = e / pdg; r_dcs = br / (pdg + dcs)
    fx = ne / nd; efx = math.sqrt(fx * (1 - fx) / nd)
    bk = BKG[sp]; nsig = (bk[0] - 100 * br) / math.sqrt(max(bk[1] ** 2 * 0 + (100 * e) ** 2 - bk[1] ** 2, 1e-12))   # podpróbka: var(różnicy) = σ_fid² − σ_all²
    O[sp] = dict(N_excl=ne, N_all=na, N_decayType=nd, BR_mc_fid=br, BR_mc_fid_err=e, f_excl=fx, f_excl_err=efx, ratio_noDCS=r, ratio_noDCS_err=er, ratio_withDCS=r_dcs,
                 ratio_err_incl_PDG=r * q(e / br, dp / pdg), bkgcat_eta25=bk[0] / 100, bkgcat_ratio=bk[0] / (100 * pdg), bkgcat_fid_col=bk[3] / 100, diff_sigma=nsig)
    MD.append(f"| {sp} | {ne} / {na} | {100 * br:.3f} ± {100 * e:.3f} | {100 * fx:.2f} ± {100 * efx:.2f} | {100 * pdg:.3f} ± {100 * dp:.3f} | **{r:.3f} ± {er:.3f}** (z błędem PDG ± {O[sp]['ratio_err_incl_PDG']:.3f}) | {r_dcs:.3f} | "
              f"{bk[0]:.3f} ± {bk[1]:.3f} → {bk[0] / (100 * pdg):.3f} | {bk[3]:.3f} | {nsig:+.1f} |")
MD.append("\nOdczyt: obie analizy dają w obszarze fiducjalnym TĘ SAMĄ liczbę (3,955 % i 9,336 %). 1,010 (BKGCAT) = inny obszar (|η| < 2,5, wszystkie pT; różnica 0,5–0,7 % na poziomie ~1σ stat. podpróbki); "
          "1,001 (CVER, D⁰) = dzielenie przez PDG z DCS (3,951 %), choć licznik MC nie zawiera DCS (DCS ma decayType −1). FSR jest w obu (PDG zawiera FSR). "
          "Spójna definicja dla §6.13 (mianownik C bez DCS, B = 3,936 %): **D⁰ 1,005 ± 0,006, D⁺ 1,003 ± 0,006**.")

# ============================================================ 2. przewidywane uzyski tagowane (§6.11)
h("2. Przewidywane uzyski tagowane w 427929 (PRZEWIDYWANIA; żadnych uzysków zmierzonych w zdarzeniach z protonem)")
sdv = J("SD_var_2026-11-19/out/results.json"); a3 = sdv["scenarios"]["A3off"]; unc = sdv["uncertainty"]
NINC = {"D0": [6773, 2010, 813], "Dp": [4999, 1833, 767], "Ds": [513, 265, 107]}       # TAGREG_2026-11-18/RAPORT.md §3: N_incl WP80 (zmierzone uzyski INKLUZYWNE)
AE_OLD, AE_NEW, dAE, EXI = 0.55, 0.92, 0.05, 0.99                                     # AFP_akceptancja_2026-11-20/PLAN.md §0 pkt 4, §5; ε_ξcal jak w TAGREG §3
O = OUT["2_pred"] = {"formula": "N = Σ_b N_incl,b(WP80, zmierzone inkluzywnie) × R_b(Pythia A3) × A·ε_AFP × ε_ξcal; C nie wchodzi → poprawka 0,734 / 0,939 NIE dotyczy"}
MD.append("N = Σ_b N_incl,b × R_b × A·ε_AFP × ε_ξcal (TAGREG §3 scenariusz M; SD_var `analyze_var.py` w. 21). **N_incl to zmierzone uzyski inkluzywne, więc C ani BR nie wchodzą: czynnik 0,734 / 0,939 nie ma zastosowania.**\n")
MD.append("| gat. | R (2–30 GeV) [%] | R_b 3–5 / 5–8 / 8–30 [%] | obwiednia | N przy A·ε = 0,55 | **N przy A·ε = 0,92** (obwiednia modeli) | ± od A·ε (5,4 %) | P = R × A·ε × ε_ξcal [%] |")
MD.append("|---|---|---|---|---|---|---|---|")
for k, nm in (("D0", "D⁰"), ("Dp", "D⁺"), ("Ds", "D*⁺")):
    rb = a3[k]["R_pt"][1:]; n0 = sum(n * r for n, r in zip(NINC[k], rb)); u = unc[k]; fu, fd = u["factor_up"], u["factor_dn"]
    old = n0 * AE_OLD * EXI; new = n0 * AE_NEW * EXI; P = a3[k]["R"] * AE_NEW * EXI
    O[k] = dict(R=a3[k]["R"], R_b=rb, f_up=fu, f_dn=fd, N_old=old, N_old_lo=old * fd, N_old_hi=old * fu, N_new=new, N_new_lo=new * fd, N_new_hi=new * fu, dN_AE=new * dAE / AE_NEW,
                P=P, P_lo=P * fd, P_hi=P * fu)
    MD.append(f"| {nm} | {100 * a3[k]['R']:.3f} | {' / '.join(f'{100 * x:.3f}' for x in rb)} | × {fd:.2f} … × {fu:.2f} | {old:.1f} ({old * fd:.1f}–{old * fu:.0f}) | **{new:.1f} ({new * fd:.1f}–{new * fu:.0f})** | ± {new * dAE / AE_NEW:.1f} | {100 * P:.2f} ({100 * P * fd:.2f}–{100 * P * fu:.2f}) |")
assert abs(O["D0"]["N_old"] - sdv["scenarios"]["A3off"]["N427929_D0_3_30"]) < 0.05, "nie odtwarzam 18,6 z SD_var"
fin = J("XSEC_final_2026-11-20/out/final.json"); tx = J("TAGXSEC_2026-11-20/out/tagxsec_excl.json")
m = tx["mc"]["D0"]; nom = tx["wynik"]["D0"]["nom"]; var = tx["wynik"]["D0"]["var"]
N70 = [fin["cells"][f"D0|{b}"]["N"] for b in BN]; rb = a3["D0"]["R_pt"][1:]
AEFF, ENPU, EXI1 = nom["AEeff"], 0.985, 1.000
n70 = sum(n * r * rho for n, r, rho in zip(N70, rb, m["rho_b"])) * AEFF * ENPU * EXI1
fu, fd = unc["D0"]["factor_up"], unc["D0"]["factor_dn"]
py_abs = tx["pred"]["D0"]["R_soft"] * tx["pred"]["D0"]["py_incl_ub"]                      # 0,353 % × 1871 µb = 6,6 µb
n_abs_tx = py_abs * nom["den_nb"] * 1e3                                                   # TAGXSEC §4: 6,6 µb × 4,02 µb⁻¹ ≈ 27 (C ekskl., BR 3,951 %, SF(3–5) = 1)
O["D0_WP70"] = dict(N_from_measured_incl=n70, lo=n70 * fd, hi=n70 * fu, N_abs_pythia_TAGXSEC=n_abs_tx, lo_abs=n_abs_tx * tx["pred"]["D0"]["R_soft_lo"] / tx["pred"]["D0"]["R_soft"],
                     hi_abs=n_abs_tx * tx["pred"]["D0"]["R_soft_hi"] / tx["pred"]["D0"]["R_soft"], ratio=n70 / n_abs_tx,
                     data_over_pythia_incl=sum(fin["cells"][f"D0|{b}"]["sigma_ub"] for b in BN) / sum(fin["cells"][f"D0|{b}"]["pythia_ub"] for b in BN))
MD.append(f"\nPorównanie z TAGXSEC (WP70, D⁰ 3–30 GeV): (a) normalizacja do ZMIERZONEGO uzysku inkluzywnego: Σ N_b(WP70, końcowe 4059 / 1536 / 700) × R_b × ρ_b × ⟨A·ε_AFP⟩ {AEFF:.3f} × ε_npu 0,985 = **{n70:.1f} ({n70 * fd:.1f}–{n70 * fu:.0f})**; "
          f"(b) TAGXSEC: bezwzględny przekrój Pythii {py_abs:.2f} µb × mianownik {nom['den_nb'] * 1e3:.3f} µb⁻¹ (C ekskluzywne) = **{n_abs_tx:.1f} ({O['D0_WP70']['lo_abs']:.0f}–{O['D0_WP70']['hi_abs']:.0f})** (raport: ~27 (11–105)). "
          f"Stosunek (a)/(b) = {n70 / n_abs_tx:.2f} ≈ dane/Pythia inkluzywnie ({O['D0_WP70']['data_over_pythia_incl']:.2f}). "
          f"Liczba noty (WP80, {O['D0']['N_new']:.0f}) jest wyższa od (a), bo N_incl(WP80) z TAGREG = 9596 wobec 6295 przy WP70 (fity starego modelu M0, 3–5 GeV z χ² 223/33 — TAGREG §Zastrzeżenia 5).")

# ============================================================ 3. σ_tag^fid przy końcowym σ_incl i SF
h("3. Dodatek tagowany: σ_tag^fid i R^fid przy końcowych wejściach XSEC_final")
g, cN, eN, rhob = m["g_tag"], m["cN"], m["eN"], m["rho_b"]
SF0, dSF0 = m["SF"], m["dSF"]
C = [fin["cells"][f"D0|{b}"] for b in BN]
SF1 = [c["SF"] for c in C]; dSF1 = [c["SF_err"] for c in C]
assert all(abs(c["C"] - x) < 1e-12 for c, x in zip(C, cN)), "C_ekskl TAGXSEC ≠ XSEC_final"
BR0, BR1, dBR = 0.03951, 0.03936, 0.00030
L, dL = 0.4648, 0.0103
prep = J("TAGGED_2026-11-19/out/prep.json")["items"]["D0|427929|MBTS|3-30|0.7|T2"]["incl"]   # tylko uzyski INKLUZYWNE (wagi wiersza „obciążenie fitu”), jak w tx_ana.py
wincl = [x["N"] for x in prep]; wincl = [x / sum(wincl) for x in wincl]
sig0 = m["sig_incl"]; sig1 = [c["sigma_ub"] for c in C]


def cbar(SF): return sum(a * b * c * d for a, b, c, d in zip(g, cN, SF, rhob))


def ratio_with(SFv, cv_, Nb): return (1 / sum(a * b * c * d for a, b, c, d in zip(g, cv_, SFv, rhob))) / sum(n / (c * s) for n, c, s in zip(Nb, cv_, SFv))


Cb0, Cb1 = cbar(SF0), cbar(SF1)
assert abs(Cb0 - nom["Cbar"]) < 1e-9
rest = nom["den_nb"] / (Cb0 * BR0)                                           # L × ⟨A·ε_AFP⟩ × ε_trig × ε_ξcal × ε_npu [nb⁻¹]
den1 = Cb1 * BR1 * rest
contrib1 = [a * b * c * d / Cb1 for a, b, c, d in zip(g, cN, SF1, rhob)]
N1 = [c["N"] for c in C]
# N_incl użyte w TAGXSEC (XSEC 19.11): σ_incl,b = N_b /(C_b SF_b BR0 L) → N_b
N0 = [s * 1e3 * c * sf * BR0 * 0.46478 for s, c, sf in zip(sig0, cN, SF0)]
r0, r1 = ratio_with(SF0, cN, N0), ratio_with(SF1, cN, N1)
r_sf_only, r_N_only = ratio_with(SF1, cN, N0), ratio_with(SF0, cN, N1)
O = OUT["3_tag"] = {"Cbar_old": Cb0, "Cbar_new": Cb1, "contrib_new": contrib1, "den_old_ub-1": nom["den_nb"] * 1e3, "den_new_ub-1": den1 * 1e3, "sigma_incl_old": sum(sig0), "sigma_incl_new": sum(sig1),
                    "R_change_total": r1 / r0 - 1, "R_change_SF_only": r_sf_only / r0 - 1, "R_change_Nincl_only": r_N_only / r0 - 1,
                    "sigma_change_SF": Cb0 / Cb1 - 1, "sigma_change_BR": BR0 / BR1 - 1}
# budżety
bias1 = sum(w * c["rel"]["fit_bias"] for w, c in zip(wincl, C))
sf_sig = q(*[k * d / s for k, d, s in zip(contrib1, dSF1, SF1)])
mcst = q(*[k * e / c for k, e, c in zip(contrib1, eN, cN)])
npr = sum(k * c["rel"]["nonprompt"] for k, c in zip(contrib1, C))
si = sum(sig1)
inc_stat = q(*[s * c["rel"]["stat"] for s, c in zip(sig1, C)]) / si
inc_model = q(*[s * c["rel"]["fit_model"] for s, c in zip(sig1, C)]) / si
inc_alt = q(*[s * c["rel"]["fit_altmin"] for s, c in zip(sig1, C)]) / si
inc_peak = q(*[s * q(c["rel"]["peak_tmpl"], c["rel"]["peak_norm"]) for s, c in zip(sig1, C)]) / si
dsf, dc = [], []
for i in range(3):
    v = list(SF1); v[i] += dSF1[i]; dsf.append(ratio_with(v, cN, N1) / r1 - 1)
    v = list(cN); v[i] += eN[i]; dc.append(ratio_with(SF1, v, N1) / r1 - 1)
sf_R, mc_R = q(*dsf), q(*dc)
b0, bR0 = nom["budget"], nom["budget_R"]
for tag, v in (("nom", nom), ("var", var)):
    Nn = v["N"] - v["N_acc"]; s_old = v["sigma_ub"]; s_new = Nn / den1 / 1e3
    R_old = v["R"]; R_new = s_new / si
    shape = v["budget"]["kształt fitu tagowanego (TAGGED)"] or 0.0
    bud = dict(v["budget"]); bud["obciążenie fitu (domknięcie inkluzywne, liniowo)"] = bias1; bud["SF klasyfikatora (3–5 GeV: scenariusz)"] = sf_sig
    bud["stat. MC C_b (mb_nd)"] = mcst; bud["nieprompt (z XSEC)"] = npr; bud["BR"] = dBR / BR1
    syst = q(*[x for x in bud.values() if x is not None])
    budR = dict(v["budget_R"]); budR["stat. N_incl"] = inc_stat; budR["model tła N_incl"] = inc_model; budR["drugie minimum N_incl (nowy wiersz)"] = inc_alt
    budR["szablon pikujący N_incl (nowy wiersz)"] = inc_peak; budR["SF (reszta po skróceniu)"] = sf_R; budR["stat. MC C_b (reszta)"] = mc_R
    systR = q(*[x for x in budR.values() if x is not None])
    slo, shi = v["stat_lo"] / Nn, v["stat_hi"] / Nn
    P = tx["pred"]["D0"]
    sl = R_new * q(slo, systR)
    O[tag] = dict(N=v["N"], sigma_old=s_old, sigma_new=s_new, stat_lo=s_new * slo, stat_hi=s_new * shi, syst_rel_old=v["syst_rel"], syst_rel_new=syst, syst_new=s_new * syst,
                  R_old=R_old, R_new=R_new, R_stat_lo=R_new * slo, R_stat_hi=R_new * shi, R_syst_rel_old=v["R_syst_rel"], R_syst_rel_new=systR, R_syst_new=R_new * systR,
                  budget_sigma=bud, budget_R=budR, pull_soft=(R_new - P["R_soft"]) / sl, pull_hard=(R_new - P["R_hard"]) / sl, pull_soft_old=v["pull_soft"], pull_hard_old=v["pull_hard"],
                  syst_noshape=q(*[x for k, x in bud.items() if x is not None and not k.startswith("kształt")]),
                  systR_noshape=q(*[x for k, x in budR.items() if x is not None and not k.startswith("kształt")]))
MD.append(f"Wejścia: SF_b 1 ± 0,10 / 1,005 / 0,965 → {SF1[0]:.3f} ± {dSF1[0]:.3f} / {SF1[1]:.3f} / {SF1[2]:.3f}; B 3,951 % → 3,936 % (mianownik C bez DCS, jak XSEC_final); "
          f"N_incl 4105 / 1530 / 681 → {N1[0]:.0f} / {N1[1]:.0f} / {N1[2]:.0f}; σ_incl(3–30) {sum(sig0):.1f} → **{si:.1f} µb**. C_b, g_b, ρ_b, ⟨A·ε_AFP⟩, ε, L bez zmian.\n")
MD.append(f"C̄_tag {Cb0:.4f} → {Cb1:.4f} (udziały binów {' / '.join(f'{100 * x:.0f}' for x in contrib1)} %); mianownik wzoru {nom['den_nb'] * 1e3:.3f} → {den1 * 1e3:.3f} µb⁻¹ "
          f"(SF: σ × {Cb0 / Cb1:.4f}; B: σ × {BR0 / BR1:.4f}).\n")
MD.append("| | σ_tag^fid stare [µb] | **σ_tag^fid nowe [µb] (−stat +stat ± syst)** | R^fid stare [%] | **R^fid nowe [%]** | syst σ: stare → nowe [%] | syst R: stare → nowe [%] | odchylenie od Pythii miękkiej / twardej [σ]: stare → nowe |")
MD.append("|---|---|---|---|---|---|---|---|")
for tag, nm in (("nom", "nominał"), ("var", "okno − pasma")):
    o = O[tag]
    MD.append(f"| {nm} | {o['sigma_old']:.2f} | **{o['sigma_new']:.2f} −{o['stat_lo']:.2f} +{o['stat_hi']:.2f} ± {o['syst_new']:.2f}** | {100 * o['R_old']:.3f} | **{100 * o['R_new']:.3f} −{100 * o['R_stat_lo']:.3f} +{100 * o['R_stat_hi']:.3f} ± {100 * o['R_syst_new']:.3f}** | "
              f"{100 * o['syst_rel_old']:.1f} → {100 * o['syst_rel_new']:.1f} | {100 * o['R_syst_rel_old']:.1f} → {100 * o['R_syst_rel_new']:.1f} | {o['pull_soft_old']:+.2f} / {o['pull_hard_old']:+.2f} → {o['pull_soft']:+.2f} / {o['pull_hard']:+.2f} |")
MD.append(f"\n**R^fid się nie zmienia:** B i L skracają się dokładnie; zmiana SF(3–5) 1 → 0,967 zmienia licznik (1/C̄) o {100 * (Cb0 / Cb1 - 1):+.2f} % i mianownik (σ_incl) prawie tak samo → reszta {100 * O['R_change_SF_only']:+.2f} %; "
          f"nowe uzyski inkluzywne (r = 1, szablon wprost) zmieniają σ_incl o {100 * (1 / (1 + O['R_change_Nincl_only']) - 1):+.2f} % → R o {100 * O['R_change_Nincl_only']:+.2f} %; łącznie R × {1 + O['R_change_total']:.4f} "
          f"(0,986 % → {100 * O['nom']['R_new']:.3f} %, w zaokrągleniu 0,99 % w obu).")
MD.append("\nBudżet (nominał), [%]: stare → nowe\n\n| pozycja | wchodzi do | stare | nowe |\n|---|---|---|---|")
rows = [("obciążenie fitu inkluzywnego (liniowo, wagi uzysków)", "σ", b0["obciążenie fitu (domknięcie inkluzywne, liniowo)"], bias1), ("SF klasyfikatora", "σ", b0["SF klasyfikatora (3–5 GeV: scenariusz)"], sf_sig),
        ("SF, reszta po skróceniu", "R", bR0["SF (reszta po skróceniu)"], sf_R), ("stat. N_incl (D⁰: max(bootstrap, 1,3 × HESSE))", "R", bR0["stat. N_incl"], inc_stat),
        ("model tła N_incl (M4, M9)", "R", bR0["model tła N_incl"], inc_model), ("drugie minimum N_incl", "R", 0.0, inc_alt), ("szablon pikujący N_incl (kształt ⊕ κ)", "R", 0.0, inc_peak),
        ("nieprompt", "σ", b0["nieprompt (z XSEC)"], npr), ("stat. MC C_b", "σ", b0["stat. MC C_b (mb_nd)"], mcst), ("BR", "σ", b0["BR"], dBR / BR1),
        ("**razem syst σ**", "", nom["syst_rel"], O["nom"]["syst_rel_new"]), ("**razem syst R**", "", nom["R_syst_rel"], O["nom"]["R_syst_rel_new"]),
        ("razem syst σ bez kształtu fitu tagowanego", "", nom["syst_rel_noshape"], O["nom"]["syst_noshape"]), ("razem syst R bez kształtu", "", nom["R_syst_rel_noshape"], O["nom"]["systR_noshape"])]
for nm, w, a, b in rows: MD.append(f"| {nm} | {w} | {100 * a:.1f} | {100 * b:.1f} |")
O["budget_rows"] = [(nm, a, b) for nm, w, a, b in rows]

# ============================================================ 4. wycofana skala κ × 1,364
h("4. Efekt wycofanej skali κ × 1,364 na uzyski WP70 (odczyt zapisanych fitów; nic nie fitowano)")
O = OUT["4_kappa"] = {}
MD.append("| run | bin | N(szablon H, κ × 1,364) | N(szablon H, r = 1) | sama skala [%] | N(interp., κ × 1,364) | N(interp., r = 1) | sama skala [%] | w jednostkach σ_stat |\n|---|---|---|---|---|---|---|---|---|")
for b in BN:
    c = fin["cells"][f"D0|{b}"]; d1 = c["N"] / c["N_dir_rFB"] - 1; d2 = c["N_ref_r1"] / c["N_FITBIAS"] - 1; st = max(c["boot"], 1.3 * c["hesse"])
    O[f"427929|{b}"] = dict(dir=d1, ref=d2, in_sigma=(c["N"] - c["N_dir_rFB"]) / st)
    MD.append(f"| 427929 | {b} | {c['N_dir_rFB']:.0f} | {c['N']:.0f} | {100 * d1:+.1f} | {c['N_FITBIAS']:.0f} | {c['N_ref_r1']:.0f} | {100 * d2:+.1f} | {O[f'427929|{b}']['in_sigma']:+.2f} |")
M8_428 = {"3-5": (6386, 6322, 6403, 6336, 613), "5-8": (4822, 4829, 4846, 4834, 302), "8-30": (6489, 6987, 6829, 7402, 238)}   # BKGCAT tabele/m8.md: dir r_FB, dir r=1, ref r_FB, ref r=1, HESSE(dir r=1)
for b, (a1, a2, b1, b2, he) in M8_428.items():
    O[f"428770|{b}"] = dict(dir=a2 / a1 - 1, ref=b2 / b1 - 1, in_sigma=(a2 - a1) / (1.3 * he))
    MD.append(f"| 428770 | {b}{' (komórka z M0; M8 odrzucony)' if b == '8-30' else ''} | {a1} | {a2} | {100 * (a2 / a1 - 1):+.1f} | {b1} | {b2} | {100 * (b2 / b1 - 1):+.1f} | {(a2 - a1) / (1.3 * he):+.2f} |")
MD.append("\nZakresy z XSEC_final RAPORT §1 (trzy rachunki samej skali): −0,5 / +1,4…+1,8 / +2,6…+3,1 % (427929). p(χ²) M8 w 428770 8–30: 1,1·10⁻⁷ (interp., × 1,364) → 3,2·10⁻⁵ (wprost, r = 1); FITBIAS: 3·10⁻⁴ (interp., κ_MC) — czynnik ~300–3000 w p.")

# ============================================================ 5. D+ 3–5 GeV: zgodność trzech liczb
h("5. D⁺ 3–5 GeV: zgodność SF")
lad, tr, cut, comb = (0.741, 0.165), (0.968, q(0.029, 0.012, 0.039)), (1.29, 0.19), (0.9486, 0.0478)   # SFABS §0–1 (drabina 0,741 ± 0,165), §3; XSEC_atlas §B (1/R); SFABS out/final (łączny)
O = OUT["5_sf_dplus"] = {}


def wmean(xs):
    w = [1 / e ** 2 for _, e in xs]; mu = sum(x * wi for (x, _), wi in zip(xs, w)) / sum(w); return mu, 1 / math.sqrt(sum(w)), sum((x - mu) ** 2 * wi for (x, _), wi in zip(xs, w))


MD.append("| zestaw | średnia ważona | χ²/ndf | p | czynnik skali √(χ²/ndf) | błąd × skala |\n|---|---|---|---|---|---|")
for nm, xs in (("A: trzy liczby z zadania (drabina, 1/R, nominał 0,95 ± 0,05) — nominał zawiera drabinę z wagą 0,08, więc nie są niezależne", [lad, cut, comb]),
               ("B: trzy NIEZALEŻNE wejścia (drabina, przeniesienie 0,968 ± 0,050, 1/R)", [lad, tr, cut]),
               ("C: drabina + przeniesienie (= nominał SFABS)", [lad, tr]), ("D: dwie liczby bezpośrednie (drabina, 1/R)", [lad, cut])):
    mu, e, c2 = wmean(xs); nd = len(xs) - 1; p = p_chi2(c2, nd); S = math.sqrt(c2 / nd)
    O[nm[0]] = dict(mean=mu, err=e, chi2=c2, ndf=nd, p=p, scale=S, err_scaled=e * max(S, 1))
    MD.append(f"| {nm} | {mu:.3f} ± {e:.3f} | {c2:.2f}/{nd} | {p:.3f} | {S:.2f} | {e * max(S, 1):.3f} |")
z = (cut[0] - lad[0]) / q(cut[1], lad[1]); O["lad_vs_cut"] = dict(z=z, p=p_two_sided(z))
zl, zc = (lad[0] - comb[0]) / q(lad[1], 0), (cut[0] - comb[0]) / q(cut[1], comb[1])
O["pulls_vs_nominal"] = dict(ladder=(lad[0] - 0.9486) / lad[1], cut=zc, half_diff_direct=(cut[0] - lad[0]) / 2, unweighted_mean_direct=(cut[0] + lad[0]) / 2)
w_lad = (1 / lad[1] ** 2) / (1 / lad[1] ** 2 + 1 / tr[1] ** 2)
c = fin["cells"]["Dplus|3-5"]; v8 = c["variants"]["SFABS_pm008"]
O["xsec_effect"] = dict(syst_005=c["syst_rel"], syst_008=v8["syst_rel"], syst_ub_005=c["syst_ub"], syst_ub_008=c["sigma_ub"] * v8["syst_rel"], weight_ladder_in_nominal=w_lad)
R = fin["ratio"]["3-5"]; sfr = 0.042 / 0.981
sf_new = q(sfr, math.sqrt(max((0.08 / 0.9486) ** 2 - (0.0478 / 0.9486) ** 2, 0)))
syst_R_new = R["R"] * q(R["syst_noSF"] / R["R"], sf_new)
O["ratio_effect"] = dict(syst_old=R["syst"], syst_new_approx=syst_R_new, sf_term_old=sfr, sf_term_new=sf_new)
MD.append(f"\nDrabina wobec 1/R wprost: różnica {cut[0] - lad[0]:.2f} = {z:.2f}σ (p = {p_two_sided(z):.3f}; dwie liczby korzystają z tego samego uzysku WP70 i są dodatnio skorelowane — fluktuacja N₇₀ przesuwa obie w tę samą stronę — więc {z:.1f}σ jest dolnym oszacowaniem napięcia). "
          f"Odchylenia od nominału 0,949 ± 0,048: drabina {(lad[0] - 0.9486) / lad[1]:+.2f}σ (własny błąd), 1/R {zc:+.2f}σ. Waga drabiny w nominale: {w_lad:.2f}.")
MD.append(f"\nSkutek wariantu ± 0,08 (z `final.json`, bez zmiany wartości centralnej): σ(D⁺, 3–5) syst {100 * c['syst_rel']:.1f} % → {100 * v8['syst_rel']:.1f} % ({c['syst_ub']:.0f} → {c['sigma_ub'] * v8['syst_rel']:.0f} µb); "
          f"D⁺/D⁰ (3–5): syst {R['syst']:.3f} → ≈ {syst_R_new:.3f} (człon SF⁺/SF⁰ {100 * sfr:.1f} → {100 * sf_new:.1f} %, przybliżenie: dodany w kwadraturze nadmiar błędu D⁺).")

# ============================================================ 6. rozrzut estymatorów tagu D* w 427929
h("6. Tag D* w runie 427929: zależność od estymatora (izolacja wobec SF klasyfikatora)")
O = OUT["6_dstar_spread"] = {}
# AUDYT_Dstar_tlo §3: ε w danych, ptcone40 < 2 GeV, 427929: PF, PF+expo, SBS-próg, sb_count ±1,0/±1,5/±2,0
ISO = {"3-30": (0.649, 0.029, 0.655, 0.654, (0.697, 0.710, 0.744), 0.98), "3-5": (0.701, 0.046, 0.712, 0.708, (0.784, 0.827, 0.881), 0.93),
       "5-8": (0.661, 0.036, 0.650, 0.673, (0.666, 0.634, 0.642), 1.06), "8-30": (0.350, 0.039, 0.355, 0.357, (0.353, 0.356, 0.364), 0.88)}
MD.append("Izolacja (ptcone40 < 2 GeV), SF = SF_nota × ε/ε_PF:\n\n| pT | S/B | PF (nota) | PF + exp | SBS-próg | zliczanie pasm ±1,0 / ±1,5 / ±2,0 MeV | maks. odchylenie fitów progowych | zliczanie − PF |\n|---|---|---|---|---|---|---|---|")
SB = {"3-30": 0.93, "3-5": 0.57, "5-8": 4.2, "8-30": 14}
for b, (pf, e, ex, sbs, sbc, sf) in ISO.items():
    f = sf / pf; O[f"iso|{b}"] = dict(SF_pf=sf, SF_expo=ex * f, SF_sbs=sbs * f, SF_sbc=[x * f for x in sbc], err=e * f)
    MD.append(f"| {b} | {SB[b]} | {sf:.2f} ± {e * f:.2f} | {ex * f:.2f} | {sbs * f:.2f} | {' / '.join(f'{x * f:.2f}' for x in sbc)} | {max(abs(ex - pf), abs(sbs - pf)) * f:.3f} | {(min(sbc) - pf) * f:+.2f} … {(max(sbc) - pf) * f:+.2f} |")
# AUDYT_Dstar_tlo §3 druga tabela: ε(WP70) klasyfikatora H, 427929: PF, expo, SBS-próg, sb_count ±2; SF nominalne z XSEC_final (dnet4); fit 2D z SFABS §2
CLF = {"5-8": (0.687, 0.084, 0.706, 0.640, 0.650, 1.005, 0.126, 1.10, (0.05, 0.07)), "8-30": (0.692, 0.068, 0.678, 0.723, 0.653, 0.965, 0.097, 1.074, (0.03, 0.03))}
MD.append("\nSF klasyfikatora H przy WP70 (wejście do σ(D⁰) ≥ 5 GeV), SF = SF_nom × ε/ε_PF:\n\n| pT | S/B | PF (nominał) | PF + exp | SBS-próg | zliczanie ±2 MeV | połowa rozpiętości 4 estymatorów Δm | rms 4 estymatorów | fit 2D (SFABS) | obciążenie fitu 2D z zabawek | 2D po odjęciu obciążenia |\n|---|---|---|---|---|---|---|---|---|---|---|")
for b, (pf, e, ex, sbs, sbc, sf, dsf_, s2d, bias) in CLF.items():
    f = sf / pf; vals = [sf, ex * f, sbs * f, sbc * f]; half = (max(vals) - min(vals)) / 2; mu = sum(vals) / 4; rms = math.sqrt(sum((x - mu) ** 2 for x in vals) / 3)
    c = fin["cells"][f"D0|{b}"]; sf_rel_new = q(c["rel"]["SF"], half / sf); syst_new = q(c["syst_rel"], half / sf)
    O[f"clf|{b}"] = dict(vals=vals, half_range=half, rms=rms, twoD=s2d, twoD_minus_pf=s2d - sf, twoD_bias=bias, SF_rel_old=c["rel"]["SF"], SF_rel_with_spread=sf_rel_new, syst_old=c["syst_rel"], syst_with_spread=syst_new,
                         syst_ub_old=c["syst_ub"], syst_ub_with_spread=c["sigma_ub"] * syst_new)
    MD.append(f"| {b} | {'0,42' if b == '5-8' else '1,7'} | {sf:.3f} ± {dsf_:.3f} | {ex * f:.3f} | {sbs * f:.3f} | {sbc * f:.3f} | {half:.3f} | {rms:.3f} | {s2d:.2f} | +{bias[0]:.2f}…+{bias[1]:.2f} | {s2d - bias[1]:.2f}…{s2d - bias[0]:.2f} |")
MD.append("\nOdczyt: (a) izolacja — fity progowe zgodne w ≤ 0,02; odstaje tylko zliczanie pasm, tylko przy 3–5 GeV (S/B 0,57), rośnie z szerokością okna; 5–8 i 8–30 GeV stabilne. "
          "(b) SF klasyfikatora — różnica fitu 2D (+0,10) leży w binach 5–8 i 8–30 GeV, tam gdzie estymatory izolacji są stabilne; ½–⅔ tej różnicy (5–8) i ~⅓ (8–30) to zmierzone obciążenie fitu 2D; cztery estymatory Δm rozrzucone ± 0,05 po OBU stronach nominału. "
          "→ To nie jest wykazany wspólny efekt (inne biny, inny znak wzoru); wspólne jest tylko to, że przy statystyce 427929 wybór estymatora tła w Δm zmienia wydajność z tagu D* o 0,05–0,10, czego błąd statystyczny fitu nie zawiera.")
for b in ("5-8", "8-30"):
    o = O[f"clf|{b}"]
    MD.append(f"\nGdyby dodać połowę rozpiętości estymatorów Δm do błędu SF D⁰ {b} GeV: SF {100 * o['SF_rel_old']:.1f} → {100 * o['SF_rel_with_spread']:.1f} %, syst σ(D⁰) {100 * o['syst_old']:.1f} → {100 * o['syst_with_spread']:.1f} % ({o['syst_ub_old']:.1f} → {o['syst_ub_with_spread']:.1f} µb). NIE wprowadzone (do decyzji).")

# ============================================================ 7. 428770: błędy × 1,3; D+/D0
h("7. Run 428770: konwencja błędu uzysków (§5.13) i kolumna D⁺/D⁰")
O = OUT["7_428"] = {}
Y = {"3-5": (6322, 613, 639, "M8"), "5-8": (4829, 302, 302, "M8"), "8-30": (7716, 332, 353, "M0 (zostaje)")}   # N, HESSE (BKGCAT tabele/m8.md; M0: 332), max(HESSE, bootstrap) z FITBIAS data.md (wariant 19.11)
MD.append("FITBIAS: zabawki dla 3 binów × 2 runy × 2 poziomy S/B (D1, 12 komórek): rms pull M8 = 1,28, M0 = 1,46; „HESSE zaniża błąd o 20–30 %” dla wszystkich ustawień, bez rozróżnienia runu → × 1,3 uzasadnione także dla 428770.\n")
MD.append("| bin | model | N | HESSE | max(HESSE, bootstrap) z 19.11 (wariant κ × 1,364) | **max(bootstrap, 1,3 × HESSE)** | błąd wzgl. |\n|---|---|---|---|---|---|---|")
for b, (n, he, mx, mod) in Y.items():
    boot = mx if mx > he else None; new = max(1.3 * he, mx if boot else 0)
    O[b] = dict(N=n, hesse=he, err_new=new, rel=new / n, bootstrap_known=boot)
    MD.append(f"| {b} | {mod} | {n} | {he} | {mx} | **{new:.0f}** | {100 * new / n:.1f} % |")
MD.append("\nBootstrap dla wariantu końcowego (szablon wprost, r = 1) w 428770 nie istnieje; w wariancie 19.11 bootstrap ≤ HESSE przy 3–5 i 5–8 GeV (max = HESSE), więc 1,3 × HESSE jest większą z dwóch. 8–30 GeV (M0): bootstrap 353 < 1,3 × 332 = 432.")
EJ = {"3-5": ((2.1, 1.5), (2.7, 1.0)), "5-8": ((2.8, 2.0), (2.7, 1.0)), "8-30": ((12.2, 3.6), (12.0, 1.9))}   # REG428 tabele T3: ε_J12(D⁰), ε_J12(D⁺) [%], odniesienie MBTS 427929
SH = {"3-5": ((73, 23), (70, 24)), "5-8": ((74, 22), (75, 20)), "8-30": ((79, 14), (87, 12))}                # REG428 §2: udział R3x / reszty w uzysku, D⁰ i D⁺ [%]
MD.append("\n| bin | ε_J12(D⁰) [%] | ε_J12(D⁺) [%] | ε(D⁺)/ε(D⁰) | błąd wzgl. stosunku | udział R3x / reszty (bez świetlności i bez ε): D⁰ | D⁺ |\n|---|---|---|---|---|---|---|")
for b, ((a, ea), (c_, ec)) in EJ.items():
    r = c_ / a; er = r * q(ea / a, ec / c_); O[f"epsJ12|{b}"] = dict(ratio=r, err=er, rel=er / r)
    MD.append(f"| {b} | {a} ± {ea} | {c_} ± {ec} | {r:.2f} ± {er:.2f} | {100 * er / r:.0f} % | {SH[b][0][0]} / {SH[b][0][1]} % | {SH[b][1][0]} / {SH[b][1][1]} % |")
MD.append("\nWerdykt: równość ε_J12 dla D⁰ i D⁺ jest sprawdzona tylko z dokładnością 34–81 % (stosunek 0,96–1,29, błąd 0,33–1,0); 12–24 % uzysku pochodzi z łańcuchów bez znanej wydajności, w różnych proporcjach dla gatunków (8–30 GeV: R3x 79 % wobec 87 %); "
          "D⁰ 8–30 GeV używa innego modelu fitu (M0) niż 427929. Wydajność triggera NIE jest wykazana jako wspólna → kolumny D⁺/D⁰ dla 428770 nie przywracam.")

os.makedirs(f"{TOP}/out", exist_ok=True)
json.dump(OUT, open(f"{TOP}/out/przeliczenia.json", "w"), indent=1, ensure_ascii=False)
open(f"{TOP}/out/przeliczenia.md", "w").write("# AUDYT_noty — przeliczenia (generowane przez skrypty/przeliczenia.py; nie edytować)\n" + "\n".join(MD) + "\n")
print("\n".join(MD))
