# TAGXSEC — tabele (generowane przez skrypty/tx_ana.py; nie edytować ręcznie)

**Definicja mianownika C: sam decayType 10 / 4 (jak XSEC_2026-11-19)** (DEF = xsec).


## D0

### X1. C_D per bin pT (WP70, drzewo H): mb_nd, mb_sd (wszystkie), mb_sd z prawdziwym protonem 0,03 ≤ ξ ≤ 0,10

| bin pT | mb_nd k/n | C_ND | mb_sd k/n | C_SD | ρ₁ = C_SD/C_ND | mb_sd+p k/n | C_SD+p | g_b (SD+p) | g_b (SD) | g_b (inkl. zmierzone) | SF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3-5 | 3653/22306 | 0,164 ± 0,002 | 58/327 | 0,177 ± 0,021 | 1,08 ± 0,13 | 16/87 | 0,184 ± 0,042 | 0,713 | 0,706 | 0,729 | 1,000 ± 0,100 |
| 5-8 | 1919/8706 | 0,220 ± 0,004 | 28/111 | 0,252 ± 0,041 | 1,14 ± 0,19 | 9/27 | 0,333 ± 0,091 | 0,227 | 0,224 | 0,201 | 1,005 ± 0,126 |
| 8-30 | 848/2909 | 0,292 ± 0,008 | 12/42 | 0,286 ± 0,070 | 0,98 ± 0,24 | 0/7 | 0,000 ± 0,000 | 0,060 | 0,070 | 0,070 | 0,965 ± 0,097 |

Udział prawdziwego trybu w mianowniku po decayType: mb_nd 1,000 / 1,000 / 1,000 (łącznie 1,000); mb_sd 1,000 / 1,000 / 1,000 (łącznie 1,000).

Prawdziwe D0 (wszystkie mody, obszar fid.) w mb_sd: 8450, z protonem w oknie: 2191 (f_win = 25,9 %).

### X2. ρ łącznie (3–30 GeV, wagi g_b z MC SD z protonem)

| oszacowanie | ρ | uwagi |
|---|---|---|
| 1a: mb_sd wprost, wszystkie zdarzenia SD | 1,090 ± 0,100 | 98 zrekonstruowanych |
| 1b: mb_sd wprost, proton w oknie | 1,122 ± 0,196 | 25 zrekonstruowanych |
| 1c: mb_sd wprost, zliczenie łączne (ND przeważone do widma mianownika SD) | 1,086 ± 0,110 | |
| 2: krotność nSelBL, dane — okno sygnału (NOMINAŁ) | 0,956 ± 0,012 (bootstrap) | per bin: 0,962 / 0,939 / 0,966; kandydatów w oknie per bin [221, 38, 15] |
| 2: dane — pasma boczne | 0,962 | per bin: 0,960 / 0,950 / 1,016; kandydatów [455, 43, 18] |
| 2: dane — okno − pasma | 0,962 | per bin: 0,983 / 0,929 / 0,915; ujemne przedziały: tak (tylko opisowo) |
| domknięcie: metoda 2 na krotności sygnału reco mb_sd (nSelBL) | 0,954 | per bin: 0,953 / 0,929 / 1,030 |
| domknięcie: to samo, nTracksPt500Eta25 | 0,965 | per bin: 0,964 / 0,950 / 1,019 |
| niedomknięcie ρ₁(1a)/ρ₂(SD) − 1 | 14,3 % ± 10,4 % | |
| **przyjęte ρ** | **0,956 ± 0,137** | stat 0,012, warianty 0,006, niedomknięcie 0,137; ρ₁(1a) wobec ρ₂(dane): 1,3σ |

Mediana nSelBL: dane tagowane okno 19,0, pasma 18,0; sygnał reco mb_nd 28,0, mb_sd 19,0. Strony kandydatów (tylko A / tylko C / obie): okno [137, 137, 0], pasma [269, 247, 0].
Prawdziwe D w mb_sd z protonem: mediana nTracksPt500Eta25 = 19,0 (mb_nd: 33,0); udział zdarzeń z < 5 torami 1,8 %, ≤ 2 torami 0,4 %.

### X3. C_b(n) w mb_nd (nSelBL; przedziały o ≥ 100 zrekonstruowanych)

- 3-5 GeV: n ≥ 0: 0,112; n ≥ 7: 0,157; n ≥ 9: 0,150; n ≥ 11: 0,159; n ≥ 13: 0,140; n ≥ 15: 0,163; n ≥ 17: 0,172; n ≥ 19: 0,169; n ≥ 21: 0,160; n ≥ 23: 0,164; n ≥ 25: 0,167; n ≥ 27: 0,177; n ≥ 29: 0,181; n ≥ 31: 0,163; n ≥ 33: 0,167; n ≥ 35: 0,163; n ≥ 37: 0,155; n ≥ 39: 0,174; n ≥ 41: 0,179; n ≥ 43: 0,171; n ≥ 45: 0,173; n ≥ 48: 0,198; n ≥ 51: 0,172; n ≥ 55: 0,165; n ≥ 62: 0,163
- 5-8 GeV: n ≥ 0: 0,187; n ≥ 10: 0,195; n ≥ 13: 0,215; n ≥ 16: 0,212; n ≥ 19: 0,219; n ≥ 22: 0,205; n ≥ 25: 0,224; n ≥ 28: 0,237; n ≥ 31: 0,200; n ≥ 34: 0,254; n ≥ 37: 0,236; n ≥ 40: 0,241; n ≥ 43: 0,210; n ≥ 48: 0,245; n ≥ 53: 0,238
- 8-30 GeV: n ≥ 0: 0,263; n ≥ 16: 0,296; n ≥ 23: 0,313; n ≥ 28: 0,308; n ≥ 34: 0,293; n ≥ 40: 0,286; n ≥ 48: 0,288

### X4-nom. D0: nominał (fit TAGGED)

N = 54,7 −18,7 +19,4 ± 12,6 (kształt); N_acc = 0,20; C̄_tag = 0,176 (z g_b inkluzywnym 0,176; udziały binów w C̄: 63,9 / 26,9 / 9,2 %); ⟨A·ε_AFP⟩ = 0,922; mianownik = 2,934 µb⁻¹.

**σ_tag^fid = 18,6 −6,4 +6,6 (stat) ± 6,5 (syst) µb** (syst 34,9 %; bez kształtu fitu tagowanego 26,0 %); L 7 % zamiast 2,2 %: syst 35,5 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,99 −0,34 +0,35 (stat) ± 0,29 (syst) %** (σ_incl = 1872,9 µb; syst 29,6 %, bez kształtu 18,4 %); przedział 68 % [0,54; 1,45] %, 95 % [0,11; 1,89] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | 23,2 | 23,2 |
| N_acc | 0,2 | 0,2 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 19,4 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 7,3 | – |
| ρ = C_D(tag)/C_D(ND) | 14,4 | 14,4 |
| widmo pT (g_b: MC SD z protonem wobec inkluzywnego) | 0,1 | – |
| stat. MC C_b (mb_nd) | 1,1 | – |
| nieprompt (z XSEC) | 2,5 | – |
| A (geometria) | 3,1 | 3,1 |
| ε_AFP | 2,1 | 2,1 |
| ustawienie / optyka | 3,4 | 3,4 |
| ε_trig | 0,1 | – |
| ε_ξcal (jednostronnie +) | 0,5 | 0,5 |
| ε_npu | 0,7 | 0,7 |
| ε_wierzch (udział ntrk < 5) | 1,8 | 1,8 |
| BR | 0,8 | – |
| L (2,2 %) | 2,2 | – |
| stat. N_incl | – | 8,4 |
| model tła N_incl | – | 5,5 |
| SF (reszta po skróceniu) | – | 1,3 |
| stat. MC C_b (reszta) | – | 0,2 |
| widmo pT (g_b) | – | 0,1 |
| ε_trig (stosunek) | – | 0,5 |

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,14–1,40) → zgodny (+1.4σ od centralnej); twarda R = 0,66 % (0,58–1,41) → zgodny (+0.7σ). σ_pred = R × σ_incl zmierzone: 6,6 (2,6–26,2) µb; × σ_incl Pythii A3 (1871,0 µb): 6,6 µb; twarda: 12,3 (10,9–26,3) µb.

### X4-var. D0: wariant okno − pasma

N = 29,9 −12,5 +12,5; N_acc = 0,20; C̄_tag = 0,176 (z g_b inkluzywnym 0,176; udziały binów w C̄: 63,9 / 26,9 / 9,2 %); ⟨A·ε_AFP⟩ = 0,922; mianownik = 2,934 µb⁻¹.

**σ_tag^fid = 10,1 −4,3 +4,3 (stat) ± 2,6 (syst) µb** (syst 26,0 %; bez kształtu fitu tagowanego 26,0 %); L 7 % zamiast 2,2 %: syst 26,9 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,54 −0,23 +0,23 (stat) ± 0,10 (syst) %** (σ_incl = 1872,9 µb; syst 18,4 %, bez kształtu 18,4 %); przedział 68 % [0,29; 0,79] %, 95 % [0,05; 1,03] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | – | – |
| N_acc | 0,3 | 0,3 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 19,4 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 7,3 | – |
| ρ = C_D(tag)/C_D(ND) | 14,4 | 14,4 |
| widmo pT (g_b: MC SD z protonem wobec inkluzywnego) | 0,1 | – |
| stat. MC C_b (mb_nd) | 1,1 | – |
| nieprompt (z XSEC) | 2,5 | – |
| A (geometria) | 3,1 | 3,1 |
| ε_AFP | 2,1 | 2,1 |
| ustawienie / optyka | 3,4 | 3,4 |
| ε_trig | 0,1 | – |
| ε_ξcal (jednostronnie +) | 0,5 | 0,5 |
| ε_npu | 0,7 | 0,7 |
| ε_wierzch (udział ntrk < 5) | 1,8 | 1,8 |
| BR | 0,8 | – |
| L (2,2 %) | 2,2 | – |
| stat. N_incl | – | 8,4 |
| model tła N_incl | – | 5,5 |
| SF (reszta po skróceniu) | – | 1,3 |
| stat. MC C_b (reszta) | – | 0,2 |
| widmo pT (g_b) | – | 0,1 |
| ε_trig (stosunek) | – | 0,5 |

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,14–1,40) → zgodny (+0.8σ od centralnej); twarda R = 0,66 % (0,58–1,41) → zgodny (-0.5σ). σ_pred = R × σ_incl zmierzone: 6,6 (2,6–26,2) µb; × σ_incl Pythii A3 (1871,0 µb): 6,6 µb; twarda: 12,3 (10,9–26,3) µb.

Przeliczenia do literatury z prawdy mb_sd (D0): N(pT ≥ 3)/N(pT ≥ 3,5) przy |η| < 2,5 = 1,364; N(|η| < 2,2)/N(|η| < 2,5) przy pT ≥ 3 = 0,889.

## Dplus

### X1. C_D per bin pT (WP70, drzewo H): mb_nd, mb_sd (wszystkie), mb_sd z prawdziwym protonem 0,03 ≤ ξ ≤ 0,10

| bin pT | mb_nd k/n | C_ND | mb_sd k/n | C_SD | ρ₁ = C_SD/C_ND | mb_sd+p k/n | C_SD+p | g_b (SD+p) | g_b (SD) | g_b (inkl. zmierzone) | SF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3-5 | 3245/21458 | 0,151 ± 0,002 | 38/309 | 0,123 ± 0,019 | 0,81 ± 0,12 | 7/86 | 0,081 ± 0,029 | 0,679 | 0,681 | 0,701 | 1,000 ± 0,100 |
| 5-8 | 1827/8440 | 0,216 ± 0,004 | 17/99 | 0,172 ± 0,038 | 0,79 ± 0,18 | 6/25 | 0,240 ± 0,085 | 0,250 | 0,237 | 0,228 | 1,000 ± 0,100 |
| 8-30 | 843/2913 | 0,289 ± 0,008 | 9/42 | 0,214 ± 0,063 | 0,74 ± 0,22 | 3/9 | 0,333 ± 0,157 | 0,070 | 0,081 | 0,071 | 1,000 ± 0,100 |

Udział prawdziwego trybu w mianowniku po decayType: mb_nd 1,000 / 1,000 / 1,000 (łącznie 1,000); mb_sd 1,000 / 1,000 / 1,000 (łącznie 1,000).

Prawdziwe Dplus (wszystkie mody, obszar fid.) w mb_sd: 4494, z protonem w oknie: 1210 (f_win = 26,9 %).

### X2. ρ łącznie (3–30 GeV, wagi g_b z MC SD z protonem)

| oszacowanie | ρ | uwagi |
|---|---|---|
| 1a: mb_sd wprost, wszystkie zdarzenia SD | 0,799 ± 0,093 | 64 zrekonstruowanych |
| 1b: mb_sd wprost, proton w oknie | 0,783 ± 0,177 | 16 zrekonstruowanych |
| 1c: mb_sd wprost, zliczenie łączne (ND przeważone do widma mianownika SD) | 0,797 ± 0,100 | |
| 2: krotność nSelBL, dane — okno sygnału (NOMINAŁ) | 0,980 ± 0,018 (bootstrap) | per bin: 0,978 / 1,012 / 0,905; kandydatów w oknie per bin [81, 15, 3] |
| 2: dane — pasma boczne | 0,973 | per bin: 0,966 / 1,007 / 0,920; kandydatów [227, 27, 13] |
| 2: dane — okno − pasma | 1,014 | per bin: 1,016 / 1,017 / 1,000; ujemne przedziały: tak (tylko opisowo) |
| domknięcie: metoda 2 na krotności sygnału reco mb_sd (nSelBL) | 0,943 | per bin: 0,906 / 1,018 / 0,932 |
| domknięcie: to samo, nTracksPt500Eta25 | 0,960 | per bin: 0,930 / 1,023 / 0,943 |
| niedomknięcie ρ₁(1a)/ρ₂(SD) − 1 | -15,3 % ± 9,9 % | |
| **przyjęte ρ** | **0,980 ± 0,151** | stat 0,018, warianty 0,007, niedomknięcie 0,150; ρ₁(1a) wobec ρ₂(dane): 1,9σ |

Mediana nSelBL: dane tagowane okno 21,0, pasma 20,0; sygnał reco mb_nd 29,0, mb_sd 19,0. Strony kandydatów (tylko A / tylko C / obie): okno [56, 42, 1], pasma [133, 134, 0].
Prawdziwe D w mb_sd z protonem: mediana nTracksPt500Eta25 = 19,0 (mb_nd: 33,0); udział zdarzeń z < 5 torami 1,8 %, ≤ 2 torami 0,2 %.

### X3. C_b(n) w mb_nd (nSelBL; przedziały o ≥ 100 zrekonstruowanych)

- 3-5 GeV: n ≥ 0: 0,101; n ≥ 8: 0,143; n ≥ 10: 0,147; n ≥ 12: 0,167; n ≥ 14: 0,141; n ≥ 16: 0,134; n ≥ 18: 0,143; n ≥ 20: 0,145; n ≥ 22: 0,137; n ≥ 24: 0,156; n ≥ 26: 0,155; n ≥ 28: 0,153; n ≥ 30: 0,151; n ≥ 32: 0,164; n ≥ 34: 0,146; n ≥ 36: 0,168; n ≥ 38: 0,150; n ≥ 40: 0,146; n ≥ 42: 0,188; n ≥ 44: 0,173; n ≥ 46: 0,160; n ≥ 49: 0,150; n ≥ 53: 0,164; n ≥ 58: 0,184
- 5-8 GeV: n ≥ 0: 0,182; n ≥ 10: 0,219; n ≥ 14: 0,228; n ≥ 17: 0,201; n ≥ 20: 0,256; n ≥ 23: 0,220; n ≥ 26: 0,213; n ≥ 29: 0,231; n ≥ 31: 0,171; n ≥ 34: 0,235; n ≥ 37: 0,222; n ≥ 40: 0,242; n ≥ 43: 0,204; n ≥ 47: 0,211; n ≥ 53: 0,227
- 8-30 GeV: n ≥ 0: 0,268; n ≥ 16: 0,249; n ≥ 23: 0,270; n ≥ 29: 0,305; n ≥ 34: 0,290; n ≥ 40: 0,333; n ≥ 47: 0,316

### X4-nom. Dplus: nominał (fit TAGGED) — OPISOWO (kontrolki D⁺ niespełnione)

N = 29,1 −10,8 +11,5 ± 8,8 (kształt); N_acc = 0,17; C̄_tag = 0,174 (z g_b inkluzywnym 0,172; udziały binów w C̄: 57,8 / 31,6 / 10,6 %); ⟨A·ε_AFP⟩ = 0,921; mianownik = 6,822 µb⁻¹.

**σ_tag^fid = 4,2 −1,6 +1,7 (stat) ± 1,5 (syst) µb** (syst 36,4 %; bez kształtu fitu tagowanego 20,0 %); L 7 % zamiast 2,2 %: syst 37,0 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,55 −0,20 +0,22 (stat) ± 0,19 (syst) %** (σ_incl = 773,1 µb; syst 35,2 %, bez kształtu 17,6 %); przedział 68 % [0,27; 0,84] %, 95 % [-0,00; 1,12] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | 30,4 | 30,4 |
| N_acc | 0,3 | 0,3 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 8,8 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 6,7 | – |
| ρ = C_D(tag)/C_D(ND) | 15,4 | 15,4 |
| widmo pT (g_b: MC SD z protonem wobec inkluzywnego) | 0,9 | – |
| stat. MC C_b (mb_nd) | 1,2 | – |
| nieprompt (z XSEC) | 1,2 | – |
| A (geometria) | 3,1 | 3,1 |
| ε_AFP | 2,1 | 2,1 |
| ustawienie / optyka | 3,4 | 3,4 |
| ε_trig | 0,1 | – |
| ε_ξcal (jednostronnie +) | 0,5 | 0,5 |
| ε_npu | 0,7 | 0,7 |
| ε_wierzch (udział ntrk < 5) | 1,8 | 1,8 |
| BR | 1,5 | – |
| L (2,2 %) | 2,2 | – |
| stat. N_incl | – | 5,6 |
| model tła N_incl | – | 2,9 |
| SF (reszta po skróceniu) | – | 1,5 |
| stat. MC C_b (reszta) | – | 0,3 |
| widmo pT (g_b) | – | 0,9 |
| ε_trig (stosunek) | – | 0,5 |

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,13–1,32) → zgodny (+0.7σ od centralnej); twarda R = 0,68 % (0,60–1,43) → zgodny (-0.5σ). σ_pred = R × σ_incl zmierzone: 2,7 (1,0–10,2) µb; × σ_incl Pythii A3 (980,0 µb): 3,4 µb; twarda: 5,3 (4,7–11,0) µb.

### X4-var. Dplus: wariant okno − pasma — OPISOWO (kontrolki D⁺ niespełnione)

N = 31,7 −12,1 +12,1; N_acc = 0,17; C̄_tag = 0,174 (z g_b inkluzywnym 0,172; udziały binów w C̄: 57,8 / 31,6 / 10,6 %); ⟨A·ε_AFP⟩ = 0,921; mianownik = 6,822 µb⁻¹.

**σ_tag^fid = 4,6 −1,8 +1,8 (stat) ± 0,9 (syst) µb** (syst 20,0 %; bez kształtu fitu tagowanego 20,0 %); L 7 % zamiast 2,2 %: syst 21,1 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,60 −0,23 +0,23 (stat) ± 0,11 (syst) %** (σ_incl = 773,1 µb; syst 17,6 %, bez kształtu 17,6 %); przedział 68 % [0,35; 0,85] %, 95 % [0,10; 1,09] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | – | – |
| N_acc | 0,3 | 0,3 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 8,8 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 6,7 | – |
| ρ = C_D(tag)/C_D(ND) | 15,4 | 15,4 |
| widmo pT (g_b: MC SD z protonem wobec inkluzywnego) | 0,9 | – |
| stat. MC C_b (mb_nd) | 1,2 | – |
| nieprompt (z XSEC) | 1,2 | – |
| A (geometria) | 3,1 | 3,1 |
| ε_AFP | 2,1 | 2,1 |
| ustawienie / optyka | 3,4 | 3,4 |
| ε_trig | 0,1 | – |
| ε_ξcal (jednostronnie +) | 0,5 | 0,5 |
| ε_npu | 0,7 | 0,7 |
| ε_wierzch (udział ntrk < 5) | 1,8 | 1,8 |
| BR | 1,5 | – |
| L (2,2 %) | 2,2 | – |
| stat. N_incl | – | 5,6 |
| model tła N_incl | – | 2,9 |
| SF (reszta po skróceniu) | – | 1,5 |
| stat. MC C_b (reszta) | – | 0,3 |
| widmo pT (g_b) | – | 0,9 |
| ε_trig (stosunek) | – | 0,5 |

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,13–1,32) → zgodny (+1.0σ od centralnej); twarda R = 0,68 % (0,60–1,43) → zgodny (-0.3σ). σ_pred = R × σ_incl zmierzone: 2,7 (1,0–10,2) µb; × σ_incl Pythii A3 (980,0 µb): 3,4 µb; twarda: 5,3 (4,7–11,0) µb.

Przeliczenia do literatury z prawdy mb_sd (Dplus): N(pT ≥ 3)/N(pT ≥ 3,5) przy |η| < 2,5 = 1,367; N(|η| < 2,2)/N(|η| < 2,5) przy pT ≥ 3 = 0,891.
