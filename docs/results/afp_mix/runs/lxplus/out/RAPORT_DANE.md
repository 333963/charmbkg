# AFP-MIX na danych — WSTĘPNE

**STATUS: WSTĘPNE.** Wygenerowano 2026-11-17 02:57 przez `AFP_mix/run_data.py` (czas 425 s). Wejście: `/home/k/Documents/charm/charm_work/AFP_mix/runs/lxplus`.
Metoda: `charmbkg/afp_mix.py`, `docs/afp_mixing.md` (charmbkg). Kandydaci **bez score** (luźna selekcja v2.0, okna masy poniżej) — liczby pod pikiem to wszyscy kandydaci w oknie (sygnał + kombinatoryka), bez fitu. Raport zawiera wyłącznie frakcje tła przypadkowego tagu i liczności tagów; uzysków sygnału dyfrakcyjnego ani przekrojów nie liczono (region sygnału głównego łańcucha).
„Frakcja kontrolna” = (1 − f_acc)/(1 − f̄) — podana bez interpretacji.

Konfiguracja: okno ξ [0,030; 0,100], δ = 0,02, dawcy z ξ_cal > 0.14 (po stronach regionu), warstwy (run, LB) z fallbackiem μ (n_min_pool = 200), bootstrap dawców B = 20, definicja tagu dla f̄: `nside`, tryger: `nonafp`, weto nPileupVertices == 0: biorcy tak, dawcy tak; D⁰: L_xy ≥ 0,25 mm (w ekstrakcji), D*: L_xy ≥ 0 mm.
Próg G5: nominalny sprzed zmiany p ≥ 0,01 (wspólny `p_min`; odtworzony z `docs/afp_mixing.md` §3 — p = 0,006 dawało FAIL) oraz p ≥ 0,002 (po toy). **Raportowane są oba; wybór progu pozostaje decyzją otwartą.**

## 0. Werdykt bramek (G0–G5, region A∨C)

„moduł” = bramki dokładnie jak w `afp_mix` (χ² z wariancją z obserwowanych liczb); „Pearson” = G1 i G5 z wariancją z liczb oczekiwanych (odporne na puste przedziały; reszta bramek bez zmian). W nawiasie bramki, które nie przeszły.

| próbka | moduł, G5 p>=0.002 (po toy) | moduł, G5 p>=0.01 (nominalny sprzed zmiany) | Pearson, G5 p>=0.002 (po toy) | Pearson, G5 p>=0.01 (nominalny sprzed zmiany) |
|---|---|---|---|---|
| data_427929 | **FAIL** (G1_donor_flatness, G5_delta_shape_A, G5_delta_shape_C) | **FAIL** (G1_donor_flatness, G5_delta_shape_A, G5_delta_shape_C) | **FAIL** (G5_delta_shape_C) | **FAIL** (G5_delta_shape_C) |
| data_428770 | **FAIL** (G2_anti_xical_closure, G5_delta_shape_A, G5_delta_shape_C) | **FAIL** (G2_anti_xical_closure, G5_delta_shape_A, G5_delta_shape_C) | **FAIL** (G2_anti_xical_closure) | **FAIL** (G2_anti_xical_closure, G5_delta_shape_C) |

## data_427929 (run 427929)

### Wejścia i selekcja

- zdarzeń 28,962,571; GRL 23,874,370; LB 81–437 (249 bloków); μ (5/50/95 %) = 0.0024/0.0048/0.0250; zdarzeń z μ ≤ 0: 28,318
- tryger: HLT_noalg_L1AFP_A_OR_C JEST w tabeli; zdarzenia wyłącznie z łańcuchów L1AFP: 101 (0.00 %; częstość tagu A∨C w nich 0.941 vs 0.030 w wybranych) — **wykluczone** u biorców i dawców; bez żadnego zapisanego łańcucha: 129,534
- po GRL ∧ tryger: 23,850,567; z nPileupVertices == 0: 23,626,322; dawcy 23,626,322; biorcy (≥ 1 kandydat slim) 3,410,405; kandydaci w biorcach 13,062,594
- udział zdarzeń per łańcuch: HLT_j20_L1MBTS_2 3.11 %, HLT_mb_mbts_L1MBTS_1 99.53 %, HLT_mb_sptrk_L1MBTS_1 99.55 %, HLT_mb_sptrk_pt2_L1MBTS_2 37.03 %, HLT_mb_sptrk_pt4_L1AFP_A_OR_C 0.04 %, HLT_mb_sptrk_pt4_L1MBTS_2 6.98 %, HLT_noalg_L1AFP_A_OR_C 2.17 %, HLT_noalg_L1J12 0.26 %, HLT_noalg_L1RD0_FILLED 0.00 %, HLT_noalg_mb_L1MBTS_1 99.55 %, HLT_noalg_mb_L1MBTS_2 98.63 %
- liczności tagów (przed cięciem ξ_cal): dawcy A 363,795 / C 347,165 / A∧C 7,713 z 23,626,322; biorcy A 10,067 / C 9,541 / A∧C 69 z 3,410,405
- zgodność `afp_tag` z definicją: A: tagów 464,222, różnic vs n_side==1 0, vs „1 w oknie” 2,739; C: tagów 443,977, różnic vs n_side==1 0, vs „1 w oknie” 2,443

### Bramki

| bramka | status | kluczowe |
|---|---|---|
| G0_inputs | PASS | coverage=1, frac_lb_strata=1, pool=1.717e+06, n_afp_only_excluded=98 |
| G1_donor_flatness [AorC] | FAIL (Pearson: PASS) | p=1.851e-95, f_pool=1.226e-07, f_below_veto=0.03208, chi2=441.9, ndf=3, p_pearson=0.3058, chi2_pearson=2.369, ndf_pearson=2, min_expected=9.153 |
| G1_donor_flatness [A] | PASS (Pearson: PASS) | p=0.8395, f_pool=0.0001366, f_below_veto=0.01813, chi2=0.8415, ndf=3, p_pearson=0.7368, chi2_pearson=0.6108, ndf_pearson=2, min_expected=34.18 |
| G1_donor_flatness [C] | PASS (Pearson: PASS) | p=0.4856, f_pool=0.0001098, f_below_veto=0.01742, chi2=2.444, ndf=3, p_pearson=0.2419, chi2_pearson=2.839, ndf_pearson=2, min_expected=30.23 |
| G2_anti_xical_closure [AorC] | PASS | n_obs=296, n_pred=287.8, ratio=1.029, ratio_err=0.08071, z=0.3544, ratio_selfexcl=1.029, z_selfexcl=0.3544 |
| G2_anti_xical_closure [A] | PASS | n_obs=247, n_pred=238.3, ratio=1.036, ratio_err=0.07263, z=0.5021, ratio_selfexcl=1.036, z_selfexcl=0.5021 |
| G2_anti_xical_closure [C] | PASS | n_obs=191, n_pred=196.8, ratio=0.9704, ratio_err=0.08308, z=-0.3558, ratio_selfexcl=0.9704, z_selfexcl=-0.3559 |
| G3_mu_slope | PASS | slope_data=0.001334, slope_mix=6.659e-05, z_slope=0.2265, diffraction_mu0=0.005474 |
| G4_xi_control | PASS | n_obs=387, n_pred=6.957, ratio=55.63, ratio_err=2.828 |
| G5_delta_shape_A | moduł: FAIL przy p>=0.002 (po toy); FAIL przy p>=0.01 (nominalny sprzed zmiany) / Pearson: PASS przy p>=0.002 (po toy); PASS przy p>=0.01 (nominalny sprzed zmiany) | p=0, norm_ratio=0.8285, z_norm=-0.7743, chi2=6.757e+05, ndf=13, eps_xical_mix=0.6355, p_pearson=0.02552, chi2_pearson=4.988, ndf_pearson=1, n_obs_ctl=14, n_tagged_donors_used=491 |
| G5_delta_shape_C | moduł: FAIL przy p>=0.002 (po toy); FAIL przy p>=0.01 (nominalny sprzed zmiany) / Pearson: FAIL przy p>=0.002 (po toy); FAIL przy p>=0.01 (nominalny sprzed zmiany) | p=0, norm_ratio=0.5094, z_norm=-2.548, chi2=4.046e+05, ndf=13, eps_xical_mix=0.6347, p_pearson=1, chi2_pearson=0, ndf_pearson=0, n_obs_ctl=7, n_tagged_donors_used=411 |
| G6_eta_asymmetry (wszyscy kandydaci slim) | PASS (informacyjna) | z=-1.431, A_obs=0.2321, A_err=0.005104, f_acc_eta=-0.1603, f_acc_eta_err=0.1188, f_acc_mix=0.00968, n_cand=3.634e+04 |
| G6_eta_asymmetry (Dplus pik) | PASS (informacyjna) | z=-1.989, A_obs=0.2574, A_err=0.01516, f_acc_eta=-0.2869, f_acc_eta_err=0.1493, f_acc_mix=0.01009, n_cand=4064 |
| G6_eta_asymmetry (D0 pik) | PASS (informacyjna) | z=0.297, A_obs=0.1907, A_err=0.0164, f_acc_eta=0.04636, f_acc_eta_err=0.1258, f_acc_mix=0.009001, n_cand=3581 |
| G6_eta_asymmetry (Dstar pik) | PASS (informacyjna) | z=0.4832, A_obs=0.1673, A_err=0.0615, f_acc_eta=0.1634, f_acc_eta_err=0.3187, f_acc_mix=0.009456, n_cand=257 |

### f_acc inkluzywnie (zdarzenia-biorcy, po cięciu ξ_cal)

| region | N_tag | N_acc (mix) | f_acc | ε_ξcal(przyp.) | f̄ (n_side) | f̄ (1 w oknie) | frakcja kontrolna | pula dawców |
|---|---|---|---|---|---|---|---|---|
| AorC | 18871 | 163.7 ± 6.5 | 0.0087 ± 0.0004 | 0.1875 | 0.0002 | 0.0002 | 0.9916 | 1,717,305 |
| A | 9706 | 89.5 ± 2.7 | 0.0092 ± 0.0003 | 0.1926 | 0.0002 | 0.0002 | 0.9910 | 3,588,801 |
| C | 9229 | 73.2 ± 2.8 | 0.0079 ± 0.0003 | 0.1948 | 0.0001 | 0.0001 | 0.9922 | 3,720,572 |

Systematyka względna N_acc (A∨C): stat. puli 0.0398, warianty puli 0.3240 (margin_0 1.156, margin_0.05 0.871, margin_0.1 0.676, strata_mu 0.997; podpróbka 200,000 biorców / 1,000,000 dawców), informacyjnie (nie w systematyce — inna selekcja): dawcy bez weta nPileupVertices → ×3.75; niedomknięcie G2 0.0000 → **razem 0.3264**. Fallback f_acc × r_G2 = 0.0089.

**Naiwne mieszanie** (dowolny dawca, bez wymogu dużego ξ_cal): N_acc = 22010.2, f_acc = 1.1663, ε = 0.2166; naiwne/nominalne = **134.43** (f_acc > 1 — niefizyczne).

### f_acc(strona, μ (actual_mu))

| μ (actual_mu) | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C | ε A∨C | f̄ A∨C |
|---|---|---|---|---|---|---|---|---|---|
| 0.001395–0.003201 | 0.002644 | 1636 | 0.0091 ± 0.0004 | 1576 | 0.0077 ± 0.0003 | 3200 | 0.0085 ± 0.0004 | 0.1867 | 0.0003 |
| 0.003201–0.003952 | 0.003577 | 1592 | 0.0093 ± 0.0004 | 1507 | 0.0082 ± 0.0004 | 3091 | 0.0089 ± 0.0004 | 0.1875 | 0.0003 |
| 0.003952–0.004842 | 0.004365 | 1628 | 0.0090 ± 0.0004 | 1546 | 0.0079 ± 0.0004 | 3163 | 0.0085 ± 0.0004 | 0.1861 | 0.0002 |
| 0.004842–0.007677 | 0.005909 | 1635 | 0.0090 ± 0.0004 | 1495 | 0.0082 ± 0.0004 | 3118 | 0.0086 ± 0.0004 | 0.1879 | 0.0002 |
| 0.007677–0.01627 | 0.01121 | 1607 | 0.0091 ± 0.0004 | 1534 | 0.0078 ± 0.0004 | 3129 | 0.0085 ± 0.0004 | 0.1870 | 0.0002 |
| 0.01627–0.0298 | 0.02325 | 1608 | 0.0098 ± 0.0004 | 1571 | 0.0079 ± 0.0004 | 3170 | 0.0090 ± 0.0004 | 0.1898 | 0.0003 |

### f_acc(strona, lumiBlock)

| lumiBlock | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C | ε A∨C | f̄ A∨C |
|---|---|---|---|---|---|---|---|---|---|
| 119–151 | 129.9 | 1200 | 0.0052 ± 0.0007 | 1130 | 0.0070 ± 0.0012 | 2318 | 0.0051 ± 0.0009 | 0.1775 | 0.0002 |
| 151–178 | 164.6 | 1152 | 0.0090 ± 0.0012 | 1147 | 0.0092 ± 0.0013 | 2294 | 0.0087 ± 0.0010 | 0.1983 | 0.0002 |
| 178–201 | 188.7 | 1182 | 0.0100 ± 0.0011 | 1071 | 0.0089 ± 0.0012 | 2248 | 0.0109 ± 0.0012 | 0.1784 | 0.0003 |
| 201–229 | 213.8 | 1284 | 0.0111 ± 0.0013 | 1200 | 0.0085 ± 0.0014 | 2472 | 0.0097 ± 0.0014 | 0.1992 | 0.0003 |
| 229–250 | 238.9 | 1153 | 0.0102 ± 0.0015 | 1141 | 0.0082 ± 0.0011 | 2288 | 0.0090 ± 0.0016 | 0.1908 | 0.0003 |
| 250–276 | 261.8 | 1284 | 0.0104 ± 0.0010 | 1234 | 0.0065 ± 0.0008 | 2513 | 0.0095 ± 0.0011 | 0.1799 | 0.0003 |
| 276–300 | 287.8 | 1226 | 0.0084 ± 0.0009 | 1125 | 0.0082 ± 0.0012 | 2342 | 0.0084 ± 0.0012 | 0.1841 | 0.0002 |
| 300–324 | 311.4 | 1225 | 0.0092 ± 0.0011 | 1181 | 0.0071 ± 0.0011 | 2396 | 0.0081 ± 0.0008 | 0.1918 | 0.0002 |

### Gatunki × pT: tag przypadkowy pod pikiem i w sidebandach (A∨C; kandydaci, bez fitu)

N_tag = kandydaci w oknie w zdarzeniach z tagiem po cięciu ξ_cal; N_acc = Σ p_acc po kandydatach (z poprawką własnego okna).

**Dplus** — Dmeson_mass: pik (1.84, 1.9), SB ((1.75, 1.8), (1.94, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | ε pik | f̄ pik | f_acc naiwne pik |
|---|---|---|---|---|---|---|---|---|---|---|
| 3–150 | 4077 | 41.1 ± 1.7 | 0.0101 ± 0.0005 | 0.0107 ± 0.0004 | 0.0092 ± 0.0004 | 0.0101 ± 0.0005 | 1.002 ± 0.064 | 0.0939 | 0.0003 | 1.3870 |
| 3–4 | 3015 | 30.7 ± 1.3 | 0.0102 ± 0.0005 | 0.0108 ± 0.0005 | 0.0094 ± 0.0005 | 0.0100 ± 0.0005 | 1.017 ± 0.066 | 0.0929 | 0.0003 | 1.4026 |
| 4–5 | 690 | 6.9 ± 0.3 | 0.0100 ± 0.0006 | 0.0102 ± 0.0007 | 0.0094 ± 0.0007 | 0.0102 ± 0.0005 | 0.977 ± 0.078 | 0.0947 | 0.0003 | 1.3670 |
| 5–6 | 205 | 2.0 ± 0.1 | 0.0096 ± 0.0008 | 0.0109 ± 0.0012 | 0.0081 ± 0.0009 | 0.0114 ± 0.0008 | 0.845 ± 0.096 | 0.0989 | 0.0003 | 1.3129 |
| 6–8 | 104 | 1.1 ± 0.1 | 0.0104 ± 0.0012 | 0.0131 ± 0.0020 | 0.0082 ± 0.0012 | 0.0096 ± 0.0009 | 1.082 ± 0.158 | 0.1043 | 0.0002 | 1.4263 |
| 8–12 | 46 | 0.4 ± 0.0 | 0.0085 ± 0.0015 | 0.0103 ± 0.0029 | 0.0068 ± 0.0016 | 0.0079 ± 0.0011 | 1.065 ± 0.247 | 0.1124 | 0.0003 | 1.1188 |
| 12–150 | 17 | 0.1 ± 0.0 | 0.0059 ± 0.0016 | 0.0116 ± 0.0052 | 0.0037 ± 0.0013 | 0.0096 ± 0.0026 | 0.621 ± 0.242 | 0.1086 | 0.0003 | 0.8178 |

**D0** — Dmeson_mass: pik (1.835, 1.895), SB ((1.75, 1.8), (1.93, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | ε pik | f̄ pik | f_acc naiwne pik |
|---|---|---|---|---|---|---|---|---|---|---|
| 3–150 | 3592 | 32.3 ± 1.3 | 0.0090 ± 0.0004 | 0.0098 ± 0.0004 | 0.0081 ± 0.0004 | 0.0091 ± 0.0004 | 0.993 ± 0.062 | 0.1323 | 0.0003 | 1.2207 |
| 3–4 | 2497 | 22.5 ± 0.9 | 0.0090 ± 0.0004 | 0.0098 ± 0.0004 | 0.0081 ± 0.0004 | 0.0091 ± 0.0004 | 0.987 ± 0.063 | 0.1338 | 0.0003 | 1.2223 |
| 4–5 | 722 | 6.5 ± 0.3 | 0.0090 ± 0.0005 | 0.0094 ± 0.0006 | 0.0084 ± 0.0006 | 0.0091 ± 0.0005 | 0.989 ± 0.075 | 0.1296 | 0.0002 | 1.2147 |
| 5–6 | 203 | 1.9 ± 0.1 | 0.0095 ± 0.0008 | 0.0102 ± 0.0011 | 0.0087 ± 0.0009 | 0.0090 ± 0.0006 | 1.056 ± 0.113 | 0.1276 | 0.0003 | 1.2878 |
| 6–8 | 121 | 1.0 ± 0.0 | 0.0084 ± 0.0009 | 0.0104 ± 0.0016 | 0.0068 ± 0.0009 | 0.0092 ± 0.0008 | 0.913 ± 0.122 | 0.1286 | 0.0002 | 1.1391 |
| 8–12 | 37 | 0.3 ± 0.0 | 0.0090 ± 0.0016 | 0.0104 ± 0.0026 | 0.0069 ± 0.0015 | 0.0066 ± 0.0008 | 1.352 ± 0.287 | 0.1306 | 0.0002 | 1.2244 |
| 12–150 | 12 | 0.1 ± 0.0 | 0.0069 ± 0.0020 | 0.0063 ± 0.0024 | 0.0082 ± 0.0037 | 0.0067 ± 0.0016 | 1.040 ± 0.397 | 0.1312 | 0.0003 | 0.9193 |

**Dstar** — Dmeson_delta_mass: pik (0.144, 0.147), SB ((0.15, 0.16),)

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | ε pik | f̄ pik | f_acc naiwne pik |
|---|---|---|---|---|---|---|---|---|---|---|
| 3–150 | 257 | 2.4 ± 0.1 | 0.0095 ± 0.0007 | 0.0096 ± 0.0009 | 0.0094 ± 0.0010 | 0.0104 ± 0.0006 | 0.908 ± 0.087 | 0.1000 | 0.0003 | 1.2979 |
| 3–4 | 155 | 1.5 ± 0.1 | 0.0100 ± 0.0009 | 0.0105 ± 0.0013 | 0.0095 ± 0.0012 | 0.0104 ± 0.0006 | 0.953 ± 0.107 | 0.0989 | 0.0003 | 1.3730 |
| 4–5 | 68 | 0.6 ± 0.0 | 0.0088 ± 0.0012 | 0.0081 ± 0.0014 | 0.0097 ± 0.0019 | 0.0108 ± 0.0008 | 0.813 ± 0.127 | 0.0989 | 0.0003 | 1.1932 |
| 5–6 | 20 | 0.2 ± 0.0 | 0.0089 ± 0.0020 | 0.0139 ± 0.0053 | 0.0061 ± 0.0017 | 0.0102 ± 0.0013 | 0.868 ± 0.226 | 0.1077 | 0.0003 | 1.2182 |
| 6–8 | 10 | 0.1 ± 0.0 | 0.0080 ± 0.0026 | 0.0055 ± 0.0019 | 0.0183 ± 0.0130 | 0.0075 ± 0.0013 | 1.070 ± 0.392 | 0.1049 | 0.0003 | 1.0843 |
| 8–12 | 3 | 0.0 ± 0.0 | 0.0085 ± 0.0049 | 0.0051 ± 0.0030 | — ± — | 0.0135 ± 0.0064 | 0.628 ± 0.472 | 0.1162 | 0.0003 | 1.1950 |
| 12–150 | 1 | 0.0 ± 0.0 | 0.0077 ± 0.0078 | 0.0037 ± 0.0037 | — ± — | 0.0139 ± 0.0139 | 0.557 ± 0.792 | 0.1601 | 0.0002 | 0.8217 |

## data_428770 (run 428770)

### Wejścia i selekcja

- zdarzeń 7,989,496; GRL 7,989,496; LB 4–323 (136 bloków); μ (5/50/95 %) = 0.0224/0.0512/0.0849; zdarzeń z μ ≤ 0: 695
- tryger: HLT_noalg_L1AFP_A_OR_C JEST w tabeli; zdarzenia wyłącznie z łańcuchów L1AFP: 888,801 (11.12 %; częstość tagu A∨C w nich 0.830 vs 0.007 w wybranych) — **wykluczone** u biorców i dawców; bez żadnego zapisanego łańcucha: 1,364,623
- po GRL ∧ tryger: 5,736,072; z nPileupVertices == 0: 5,449,549; dawcy 5,449,549; biorcy (≥ 1 kandydat slim) 4,475,334; kandydaci w biorcach 44,263,181
- udział zdarzeń per łańcuch: HLT_j20_L1MBTS_2 15.56 %, HLT_mb_mbts_L1MBTS_1 0.00 %, HLT_mb_sptrk_L1MBTS_1 0.00 %, HLT_mb_sptrk_pt2_L1MBTS_2 0.75 %, HLT_mb_sptrk_pt4_L1AFP_A_OR_C 2.77 %, HLT_mb_sptrk_pt4_L1MBTS_2 10.48 %, HLT_noalg_L1AFP_A_OR_C 8.87 %, HLT_noalg_L1J12 56.14 %, HLT_noalg_L1RD0_FILLED 0.00 %, HLT_noalg_mb_L1MBTS_1 0.02 %, HLT_noalg_mb_L1MBTS_2 0.82 %
- liczności tagów (przed cięciem ξ_cal): dawcy A 17,445 / C 16,558 / A∧C 411 z 5,449,549; biorcy A 10,609 / C 10,188 / A∧C 93 z 4,475,334
- zgodność `afp_tag` z definicją: A: tagów 953,504, różnic vs n_side==1 0, vs „1 w oknie” 3,616; C: tagów 654,531, różnic vs n_side==1 0, vs „1 w oknie” 2,117

### Bramki

| bramka | status | kluczowe |
|---|---|---|
| G0_inputs | PASS | coverage=1, frac_lb_strata=1, pool=2.275e+06, n_afp_only_excluded=8.888e+05 |
| G1_donor_flatness [AorC] | PASS (Pearson: PASS) | p=0.2029, f_pool=0.001276, f_below_veto=0.009665, chi2=4.607, ndf=3, p_pearson=0.1233, chi2_pearson=4.187, ndf_pearson=2, min_expected=191.9 |
| G1_donor_flatness [A] | FAIL (Pearson: FAIL) | p=0.001294, f_pool=0.0006304, f_below_veto=0.007069, chi2=15.72, ndf=3, p_pearson=0.001313, chi2_pearson=15.69, ndf_pearson=3, min_expected=17.79 |
| G1_donor_flatness [C] | FAIL (Pearson: FAIL) | p=0.000596, f_pool=0.0006235, f_below_veto=0.006849, chi2=17.36, ndf=3, p_pearson=0.0006938, chi2_pearson=17.04, ndf_pearson=3, min_expected=19.81 |
| G2_anti_xical_closure [AorC] | FAIL | n_obs=2431, n_pred=2662, ratio=0.9132, ratio_err=0.02591, z=-3.351, ratio_selfexcl=0.9132, z_selfexcl=-3.352 |
| G2_anti_xical_closure [A] | FAIL | n_obs=1651, n_pred=1827, ratio=0.9036, ratio_err=0.0296, z=-3.258, ratio_selfexcl=0.9036, z_selfexcl=-3.258 |
| G2_anti_xical_closure [C] | FAIL | n_obs=1664, n_pred=1839, ratio=0.9049, ratio_err=0.0278, z=-3.42, ratio_selfexcl=0.9049, z_selfexcl=-3.42 |
| G3_mu_slope | PASS | slope_data=0.0005679, slope_mix=0.0001733, z_slope=0.2791, diffraction_mu0=0.003399 |
| G4_xi_control | PASS | n_obs=251, n_pred=20.22, ratio=12.41, ratio_err=0.7834 |
| G5_delta_shape_A | moduł: FAIL przy p>=0.002 (po toy); FAIL przy p>=0.01 (nominalny sprzed zmiany) / Pearson: PASS przy p>=0.002 (po toy); PASS przy p>=0.01 (nominalny sprzed zmiany) | p=0, norm_ratio=1.005, z_norm=0.07805, chi2=7783, ndf=13, eps_xical_mix=0.1629, p_pearson=0.2709, chi2_pearson=9.917, ndf_pearson=8, n_obs_ctl=232, n_tagged_donors_used=2081 |
| G5_delta_shape_C | moduł: FAIL przy p>=0.002 (po toy); FAIL przy p>=0.01 (nominalny sprzed zmiany) / Pearson: PASS przy p>=0.002 (po toy); FAIL przy p>=0.01 (nominalny sprzed zmiany) | p=0, norm_ratio=1.163, z_norm=2.291, chi2=2.224e+04, ndf=13, eps_xical_mix=0.1598, p_pearson=0.008836, chi2_pearson=20.43, ndf_pearson=8, n_obs_ctl=268, n_tagged_donors_used=2099 |
| G6_eta_asymmetry (wszyscy kandydaci slim) | FAIL (informacyjna) | z=-2.866, A_obs=0.2664, A_err=0.003553, f_acc_eta=-0.3319, f_acc_eta_err=0.1344, f_acc_mix=0.05317, n_cand=7.363e+04 |
| G6_eta_asymmetry (Dplus pik) | FAIL (informacyjna) | z=-2.582, A_obs=0.2619, A_err=0.01031, f_acc_eta=-0.3097, f_acc_eta_err=0.1407, f_acc_mix=0.05373, n_cand=8770 |
| G6_eta_asymmetry (D0 pik) | FAIL (informacyjna) | z=-2.585, A_obs=0.2649, A_err=0.012, f_acc_eta=-0.3243, f_acc_eta_err=0.1454, f_acc_mix=0.05147, n_cand=6458 |
| G6_eta_asymmetry (Dstar pik) | PASS (informacyjna) | z=-0.745, A_obs=0.2246, A_err=0.04082, f_acc_eta=-0.1228, f_acc_eta_err=0.2329, f_acc_mix=0.05072, n_cand=570 |

### f_acc inkluzywnie (zdarzenia-biorcy, po cięciu ξ_cal)

| region | N_tag | N_acc (mix) | f_acc | ε_ξcal(przyp.) | f̄ (n_side) | f̄ (1 w oknie) | frakcja kontrolna | pula dawców |
|---|---|---|---|---|---|---|---|---|
| AorC | 16065 | 736.9 ± 12.2 | 0.0459 ± 0.0008 | 0.1293 | 0.0012 | 0.0012 | 0.9553 | 2,274,651 |
| A | 8271 | 377.2 ± 9.7 | 0.0456 ± 0.0013 | 0.1329 | 0.0008 | 0.0007 | 0.9552 | 3,276,118 |
| C | 7863 | 364.2 ± 8.4 | 0.0463 ± 0.0012 | 0.1297 | 0.0008 | 0.0007 | 0.9544 | 3,338,441 |

Systematyka względna N_acc (A∨C): stat. puli 0.0166, warianty puli 0.0543 (margin_0 1.038, margin_0.05 0.978, margin_0.1 0.946, strata_mu 1.010; podpróbka 200,000 biorców / 1,000,000 dawców), informacyjnie (nie w systematyce — inna selekcja): dawcy bez weta nPileupVertices → ×2.41; niedomknięcie G2 0.0609 → **razem 0.0833**. Fallback f_acc × r_G2 = 0.0419.

**Naiwne mieszanie** (dowolny dawca, bez wymogu dużego ξ_cal): N_acc = 4020.3, f_acc = 0.2503, ε = 0.1461; naiwne/nominalne = **5.46**.

### f_acc(strona, μ (actual_mu))

| μ (actual_mu) | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C | ε A∨C | f̄ A∨C |
|---|---|---|---|---|---|---|---|---|---|
| 0.01037–0.03488 | 0.02555 | 1380 | 0.0446 ± 0.0017 | 1343 | 0.0439 ± 0.0016 | 2709 | 0.0441 ± 0.0011 | 0.1262 | 0.0012 |
| 0.03488–0.04455 | 0.04052 | 1371 | 0.0451 ± 0.0017 | 1309 | 0.0459 ± 0.0016 | 2669 | 0.0452 ± 0.0012 | 0.1278 | 0.0012 |
| 0.04455–0.05092 | 0.04786 | 1314 | 0.0474 ± 0.0018 | 1318 | 0.0458 ± 0.0016 | 2623 | 0.0464 ± 0.0012 | 0.1285 | 0.0013 |
| 0.05092–0.05803 | 0.05421 | 1375 | 0.0460 ± 0.0017 | 1274 | 0.0479 ± 0.0017 | 2638 | 0.0469 ± 0.0012 | 0.1299 | 0.0012 |
| 0.05803–0.06925 | 0.0634 | 1427 | 0.0444 ± 0.0016 | 1288 | 0.0474 ± 0.0017 | 2707 | 0.0458 ± 0.0012 | 0.1300 | 0.0013 |
| 0.06925–0.2772 | 0.08923 | 1404 | 0.0463 ± 0.0017 | 1331 | 0.0471 ± 0.0017 | 2719 | 0.0468 ± 0.0012 | 0.1333 | 0.0012 |

### f_acc(strona, lumiBlock)

| lumiBlock | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C | ε A∨C | f̄ A∨C |
|---|---|---|---|---|---|---|---|---|---|
| 4–20 | 11.89 | 1021 | 0.0467 ± 0.0039 | 949 | 0.0534 ± 0.0041 | 1962 | 0.0495 ± 0.0029 | 0.1272 | 0.0013 |
| 20–221 | 147.1 | 1030 | 0.0426 ± 0.0031 | 1003 | 0.0420 ± 0.0030 | 2026 | 0.0428 ± 0.0028 | 0.1252 | 0.0012 |
| 221–238 | 229 | 1017 | 0.0476 ± 0.0035 | 996 | 0.0452 ± 0.0034 | 2007 | 0.0446 ± 0.0024 | 0.1261 | 0.0012 |
| 238–256 | 246.5 | 1072 | 0.0466 ± 0.0034 | 973 | 0.0445 ± 0.0031 | 2032 | 0.0437 ± 0.0026 | 0.1302 | 0.0012 |
| 256–273 | 264 | 1006 | 0.0408 ± 0.0032 | 999 | 0.0468 ± 0.0033 | 1999 | 0.0466 ± 0.0029 | 0.1299 | 0.0013 |
| 273–289 | 280.5 | 1035 | 0.0421 ± 0.0033 | 934 | 0.0488 ± 0.0041 | 1960 | 0.0451 ± 0.0035 | 0.1349 | 0.0012 |
| 289–306 | 297 | 1043 | 0.0491 ± 0.0032 | 948 | 0.0477 ± 0.0030 | 1987 | 0.0508 ± 0.0023 | 0.1341 | 0.0013 |
| 306–323 | 314.4 | 1047 | 0.0491 ± 0.0038 | 1061 | 0.0430 ± 0.0036 | 2092 | 0.0440 ± 0.0023 | 0.1268 | 0.0012 |

### Gatunki × pT: tag przypadkowy pod pikiem i w sidebandach (A∨C; kandydaci, bez fitu)

N_tag = kandydaci w oknie w zdarzeniach z tagiem po cięciu ξ_cal; N_acc = Σ p_acc po kandydatach (z poprawką własnego okna).

**Dplus** — Dmeson_mass: pik (1.84, 1.9), SB ((1.75, 1.8), (1.94, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | ε pik | f̄ pik | f_acc naiwne pik |
|---|---|---|---|---|---|---|---|---|---|---|
| 3–150 | 8797 | 472.7 ± 8.1 | 0.0537 ± 0.0013 | 0.0526 ± 0.0019 | 0.0553 ± 0.0020 | 0.0540 ± 0.0013 | 0.994 ± 0.034 | 0.0638 | 0.0013 | 0.2999 |
| 3–4 | 4045 | 237.0 ± 4.1 | 0.0586 ± 0.0016 | 0.0576 ± 0.0023 | 0.0601 ± 0.0024 | 0.0582 ± 0.0015 | 1.006 ± 0.038 | 0.0544 | 0.0013 | 0.3288 |
| 4–5 | 1818 | 95.9 ± 1.6 | 0.0527 ± 0.0018 | 0.0499 ± 0.0024 | 0.0561 ± 0.0027 | 0.0554 ± 0.0016 | 0.951 ± 0.043 | 0.0664 | 0.0013 | 0.2937 |
| 5–6 | 932 | 45.7 ± 0.8 | 0.0490 ± 0.0020 | 0.0476 ± 0.0027 | 0.0510 ± 0.0031 | 0.0490 ± 0.0017 | 1.000 ± 0.054 | 0.0763 | 0.0013 | 0.2725 |
| 6–8 | 857 | 42.9 ± 0.8 | 0.0501 ± 0.0022 | 0.0490 ± 0.0029 | 0.0513 ± 0.0032 | 0.0485 ± 0.0017 | 1.031 ± 0.058 | 0.0861 | 0.0013 | 0.2771 |
| 8–12 | 673 | 31.4 ± 0.5 | 0.0467 ± 0.0023 | 0.0478 ± 0.0033 | 0.0462 ± 0.0032 | 0.0458 ± 0.0019 | 1.020 ± 0.065 | 0.0977 | 0.0013 | 0.2581 |
| 12–150 | 472 | 19.8 ± 0.3 | 0.0419 ± 0.0026 | 0.0420 ± 0.0036 | 0.0419 ± 0.0036 | 0.0444 ± 0.0023 | 0.944 ± 0.076 | 0.1049 | 0.0013 | 0.2311 |

**D0** — Dmeson_mass: pik (1.835, 1.895), SB ((1.75, 1.8), (1.93, 2.0))

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | ε pik | f̄ pik | f_acc naiwne pik |
|---|---|---|---|---|---|---|---|---|---|---|
| 3–150 | 6471 | 333.1 ± 5.6 | 0.0515 ± 0.0012 | 0.0504 ± 0.0017 | 0.0529 ± 0.0018 | 0.0510 ± 0.0011 | 1.009 ± 0.032 | 0.0833 | 0.0013 | 0.2846 |
| 3–4 | 2789 | 147.8 ± 2.5 | 0.0530 ± 0.0015 | 0.0527 ± 0.0022 | 0.0538 ± 0.0021 | 0.0531 ± 0.0013 | 0.998 ± 0.037 | 0.0734 | 0.0013 | 0.2944 |
| 4–5 | 1492 | 76.1 ± 1.3 | 0.0510 ± 0.0017 | 0.0478 ± 0.0023 | 0.0550 ± 0.0027 | 0.0516 ± 0.0014 | 0.989 ± 0.043 | 0.0835 | 0.0013 | 0.2820 |
| 5–6 | 703 | 37.9 ± 0.6 | 0.0539 ± 0.0023 | 0.0543 ± 0.0033 | 0.0538 ± 0.0034 | 0.0500 ± 0.0017 | 1.077 ± 0.059 | 0.0927 | 0.0013 | 0.2968 |
| 6–8 | 649 | 33.8 ± 0.6 | 0.0520 ± 0.0023 | 0.0531 ± 0.0034 | 0.0512 ± 0.0032 | 0.0496 ± 0.0017 | 1.048 ± 0.059 | 0.1007 | 0.0012 | 0.2862 |
| 8–12 | 500 | 23.0 ± 0.4 | 0.0461 ± 0.0023 | 0.0427 ± 0.0029 | 0.0503 ± 0.0037 | 0.0472 ± 0.0019 | 0.977 ± 0.063 | 0.1107 | 0.0013 | 0.2531 |
| 12–150 | 338 | 14.5 ± 0.2 | 0.0428 ± 0.0026 | 0.0433 ± 0.0037 | 0.0420 ± 0.0035 | 0.0419 ± 0.0020 | 1.022 ± 0.077 | 0.1183 | 0.0013 | 0.2335 |

**Dstar** — Dmeson_delta_mass: pik (0.144, 0.147), SB ((0.15, 0.16),)

| pT [GeV] | N_tag pik | N_acc pik | f_acc pik | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | ε pik | f̄ pik | f_acc naiwne pik |
|---|---|---|---|---|---|---|---|---|---|---|
| 3–150 | 572 | 29.0 ± 0.5 | 0.0507 ± 0.0025 | 0.0493 ± 0.0034 | 0.0523 ± 0.0036 | 0.0551 ± 0.0017 | 0.920 ± 0.052 | 0.0678 | 0.0013 | 0.2834 |
| 3–4 | 226 | 12.1 ± 0.2 | 0.0537 ± 0.0038 | 0.0536 ± 0.0055 | 0.0547 ± 0.0055 | 0.0590 ± 0.0025 | 0.910 ± 0.076 | 0.0585 | 0.0013 | 0.3009 |
| 4–5 | 166 | 7.8 ± 0.1 | 0.0471 ± 0.0040 | 0.0435 ± 0.0050 | 0.0513 ± 0.0064 | 0.0555 ± 0.0026 | 0.849 ± 0.082 | 0.0672 | 0.0012 | 0.2642 |
| 5–6 | 71 | 3.7 ± 0.1 | 0.0517 ± 0.0065 | 0.0482 ± 0.0084 | 0.0570 ± 0.0102 | 0.0495 ± 0.0030 | 1.044 ± 0.145 | 0.0770 | 0.0013 | 0.2874 |
| 6–8 | 60 | 2.9 ± 0.1 | 0.0480 ± 0.0064 | 0.0534 ± 0.0104 | 0.0429 ± 0.0076 | 0.0522 ± 0.0036 | 0.920 ± 0.137 | 0.0873 | 0.0012 | 0.2672 |
| 8–12 | 26 | 1.6 ± 0.0 | 0.0624 ± 0.0128 | 0.0543 ± 0.0150 | 0.0667 ± 0.0193 | 0.0558 ± 0.0054 | 1.119 ± 0.253 | 0.1027 | 0.0013 | 0.3430 |
| 12–150 | 23 | 0.9 ± 0.0 | 0.0385 ± 0.0087 | 0.0410 ± 0.0135 | 0.0366 ± 0.0114 | 0.0414 ± 0.0052 | 0.930 ± 0.241 | 0.1181 | 0.0012 | 0.2124 |

## Uwagi i ograniczenia (WSTĘPNE)

- Brak score CharmNet/XGB w h5 v2.0 → kandydaci po luźnej selekcji ntupli; f_acc pod pikiem dotyczy mieszanki sygnał + kombinatoryka. Błędy f_acc: stat. puli dawców (bootstrap) ⊕ √Σw² liczby tagów; systematyka (warianty puli, G2) podana tylko inkluzywnie.
- Definicja tagu w h5 v2.0: dokładnie 1 proton po stronie i w oknie (`afp_tag_def`). `afp_mix.spoil_probability` zakłada „dokładnie 1 w oknie” (toy) → f̄ nominalne tutaj liczone wariantem n_side (`run_data.spoil_nside`); obie wartości w tabeli. Poprawka własnego okna C i overlay 'replace' w p_acc pozostają jak w `afp_mix` (biorcy z własnymi protonami bez tagu nie są osobno traktowani).
- `afp_mix.COLMAP['m']` = `Dmeson_m`, w h5 v2.0 kolumna to `Dmeson_mass` → nadpisane w `run_data.COLMAP_DATA` (repo charmbkg nie zmieniane).
- G2 z bootstrapem nie wyklucza biorcy z puli (przesunięcie ~1/n otagowanych dawców w warstwie) — obok podano `ratio_selfexcl`.
- G1 i G5 w wersji modułu używają wariancji z obserwowanych liczb: przedział z 0 zdarzeń dostaje wariancję ≈ 0 (G1: dominuje średnią ważoną; G5: χ² rośnie do tysięcy). Przy częstości tagu przypadkowego ~1e-4–1e-3 (dane) to artefakt — obok podano wariant Pearsona. G5 traktuje pary jako niezależne, choć otagowanych dawców jest mało (`n_tagged_donors_used`) — wariancja mieszania zaniżona.
- G4, G5 i warianty puli liczone na losowych podpróbkach (rozmiary w tabeli bramek / w tekście); G3 w przedziałach kwantylowych μ próbki.
- 428770 (physics_Main) jest zdominowany przez L1_J12; 427929 (physics_MinBias) ma preskale ≈ 1 — próbek nie łączono.
