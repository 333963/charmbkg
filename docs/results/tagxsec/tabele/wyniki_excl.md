# TAGXSEC — tabele (generowane przez skrypty/tx_ana.py; nie edytować ręcznie)

**Definicja mianownika C: prawdziwy tryb K∓π±[π±] (decayType sygnału ∧ klasa brakującej masy 0/1, jak BKGCAT)** (DEF = excl).


## D0

### X1. C_D per bin pT (WP70, drzewo H): mb_nd, mb_sd (wszystkie), mb_sd z prawdziwym protonem 0,03 ≤ ξ ≤ 0,10

| bin pT | mb_nd k/n | C_ND | mb_sd k/n | C_SD | ρ₁ = C_SD/C_ND | mb_sd+p k/n | C_SD+p | g_b (SD+p) | g_b (SD) | g_b (inkl. zmierzone) | SF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3-5 | 3653/16366 | 0,223 ± 0,003 | 58/243 | 0,239 ± 0,027 | 1,07 ± 0,12 | 16/63 | 0,254 ± 0,055 | 0,713 | 0,706 | 0,728 | 1,000 ± 0,100 |
| 5-8 | 1919/6394 | 0,300 ± 0,006 | 28/83 | 0,337 ± 0,052 | 1,12 ± 0,17 | 9/18 | 0,500 ± 0,118 | 0,227 | 0,224 | 0,201 | 1,005 ± 0,126 |
| 8-30 | 848/2143 | 0,396 ± 0,011 | 12/30 | 0,400 ± 0,089 | 1,01 ± 0,23 | 0/3 | 0,000 ± 0,000 | 0,060 | 0,070 | 0,071 | 0,965 ± 0,097 |

Udział prawdziwego trybu w mianowniku po decayType: mb_nd 0,734 / 0,734 / 0,737 (łącznie 0,734); mb_sd 0,743 / 0,748 / 0,714 (łącznie 0,742).

Prawdziwe D0 (wszystkie mody, obszar fid.) w mb_sd: 8450, z protonem w oknie: 2191 (f_win = 25,9 %).

### X2. ρ łącznie (3–30 GeV, wagi g_b z MC SD z protonem)

| oszacowanie | ρ | uwagi |
|---|---|---|
| 1a: mb_sd wprost, wszystkie zdarzenia SD | 1,079 ± 0,094 | 98 zrekonstruowanych |
| 1b: mb_sd wprost, proton w oknie | 1,174 ± 0,189 | 25 zrekonstruowanych |
| 1c: mb_sd wprost, zliczenie łączne (ND przeważone do widma mianownika SD) | 1,077 ± 0,109 | |
| 2: krotność nSelBL, dane — okno sygnału (NOMINAŁ) | 0,962 ± 0,012 (bootstrap) | per bin: 0,963 / 0,950 / 0,984; kandydatów w oknie per bin [221, 38, 15] |
| 2: dane — pasma boczne | 0,966 | per bin: 0,962 / 0,957 / 1,018; kandydatów [455, 43, 18] |
| 2: dane — okno − pasma | 0,967 | per bin: 0,980 / 0,944 / 0,949; ujemne przedziały: tak (tylko opisowo) |
| domknięcie: metoda 2 na krotności sygnału reco mb_sd (nSelBL) | 0,960 | per bin: 0,956 / 0,943 / 1,029 |
| domknięcie: to samo, nTracksPt500Eta25 | 0,967 | per bin: 0,964 / 0,957 / 1,012 |
| niedomknięcie ρ₁(1a)/ρ₂(SD) − 1 | 12,4 % ± 9,8 % | |
| **przyjęte ρ** | **0,962 ± 0,120** | stat 0,012, warianty 0,004, niedomknięcie 0,119; ρ₁(1a) wobec ρ₂(dane): 1,2σ |

Mediana nSelBL: dane tagowane okno 19,0, pasma 18,0; sygnał reco mb_nd 28,0, mb_sd 19,0. Strony kandydatów (tylko A / tylko C / obie): okno [137, 137, 0], pasma [269, 247, 0].
Prawdziwe D w mb_sd z protonem: mediana nTracksPt500Eta25 = 19,0 (mb_nd: 33,0); udział zdarzeń z < 5 torami 1,8 %, ≤ 2 torami 0,4 %.

### X3. C_b(n) w mb_nd (nSelBL; przedziały o ≥ 100 zrekonstruowanych)

- 3-5 GeV: n ≥ 0: 0,158; n ≥ 7: 0,217; n ≥ 9: 0,204; n ≥ 11: 0,219; n ≥ 13: 0,193; n ≥ 15: 0,221; n ≥ 17: 0,235; n ≥ 19: 0,226; n ≥ 21: 0,221; n ≥ 23: 0,227; n ≥ 25: 0,222; n ≥ 27: 0,238; n ≥ 29: 0,240; n ≥ 31: 0,220; n ≥ 33: 0,223; n ≥ 35: 0,225; n ≥ 37: 0,213; n ≥ 39: 0,234; n ≥ 41: 0,240; n ≥ 43: 0,231; n ≥ 45: 0,229; n ≥ 48: 0,271; n ≥ 51: 0,235; n ≥ 55: 0,230; n ≥ 62: 0,233
- 5-8 GeV: n ≥ 0: 0,270; n ≥ 10: 0,266; n ≥ 13: 0,293; n ≥ 16: 0,283; n ≥ 19: 0,294; n ≥ 22: 0,275; n ≥ 25: 0,301; n ≥ 28: 0,329; n ≥ 31: 0,282; n ≥ 34: 0,331; n ≥ 37: 0,324; n ≥ 40: 0,318; n ≥ 43: 0,292; n ≥ 48: 0,326; n ≥ 53: 0,323
- 8-30 GeV: n ≥ 0: 0,369; n ≥ 16: 0,415; n ≥ 23: 0,413; n ≥ 28: 0,396; n ≥ 34: 0,414; n ≥ 40: 0,379; n ≥ 48: 0,390

### X4-nom. D0: nominał (fit TAGGED)

N = 54,7 −18,7 +19,4 ± 12,6 (kształt); N_acc = 0,20; C̄_tag = 0,241 (z g_b inkluzywnym 0,241; udziały binów w C̄: 63,7 / 27,0 / 9,3 %); ⟨A·ε_AFP⟩ = 0,922; mianownik = 4,020 µb⁻¹.

**σ_tag^fid = 13,6 −4,7 +4,8 (stat) ± 4,6 (syst) µb** (syst 34,1 %; bez kształtu fitu tagowanego 25,0 %); L 7 % zamiast 2,2 %: syst 34,8 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,99 −0,34 +0,35 (stat) ± 0,28 (syst) %** (σ_incl = 1374,8 µb; syst 28,8 %, bez kształtu 17,0 %); przedział 68 % [0,54; 1,44] %, 95 % [0,12; 1,87] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | 23,2 | 23,2 |
| N_acc | 0,2 | 0,2 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 19,4 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 7,3 | – |
| ρ = C_D(tag)/C_D(ND) | 12,5 | 12,5 |
| widmo pT (g_b: MC SD z protonem wobec inkluzywnego) | 0,0 | – |
| stat. MC C_b (mb_nd) | 1,1 | – |
| nieprompt (z XSEC) | 2,4 | – |
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
| widmo pT (g_b) | – | 0,0 |
| ε_trig (stosunek) | – | 0,5 |

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,14–1,40) → zgodny (+1.4σ od centralnej); twarda R = 0,66 % (0,58–1,41) → zgodny (+0.7σ). σ_pred = R × σ_incl zmierzone: 4,8 (1,9–19,3) µb; × σ_incl Pythii A3 (1871,0 µb): 6,6 µb; twarda: 9,0 (8,0–19,3) µb.

### X4-var. D0: wariant okno − pasma

N = 29,9 −12,5 +12,5; N_acc = 0,20; C̄_tag = 0,241 (z g_b inkluzywnym 0,241; udziały binów w C̄: 63,7 / 27,0 / 9,3 %); ⟨A·ε_AFP⟩ = 0,922; mianownik = 4,020 µb⁻¹.

**σ_tag^fid = 7,4 −3,1 +3,1 (stat) ± 1,9 (syst) µb** (syst 25,0 %; bez kształtu fitu tagowanego 25,0 %); L 7 % zamiast 2,2 %: syst 25,9 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,54 −0,23 +0,23 (stat) ± 0,09 (syst) %** (σ_incl = 1374,8 µb; syst 17,0 %, bez kształtu 17,0 %); przedział 68 % [0,29; 0,78] %, 95 % [0,06; 1,02] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | – | – |
| N_acc | 0,3 | 0,3 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 19,4 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 7,3 | – |
| ρ = C_D(tag)/C_D(ND) | 12,5 | 12,5 |
| widmo pT (g_b: MC SD z protonem wobec inkluzywnego) | 0,0 | – |
| stat. MC C_b (mb_nd) | 1,1 | – |
| nieprompt (z XSEC) | 2,4 | – |
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
| widmo pT (g_b) | – | 0,0 |
| ε_trig (stosunek) | – | 0,5 |

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,14–1,40) → zgodny (+0.8σ od centralnej); twarda R = 0,66 % (0,58–1,41) → zgodny (-0.5σ). σ_pred = R × σ_incl zmierzone: 4,8 (1,9–19,3) µb; × σ_incl Pythii A3 (1871,0 µb): 6,6 µb; twarda: 9,0 (8,0–19,3) µb.

Przeliczenia do literatury z prawdy mb_sd (D0): N(pT ≥ 3)/N(pT ≥ 3,5) przy |η| < 2,5 = 1,364; N(|η| < 2,2)/N(|η| < 2,5) przy pT ≥ 3 = 0,889.

## Dplus

### X1. C_D per bin pT (WP70, drzewo H): mb_nd, mb_sd (wszystkie), mb_sd z prawdziwym protonem 0,03 ≤ ξ ≤ 0,10

| bin pT | mb_nd k/n | C_ND | mb_sd k/n | C_SD | ρ₁ = C_SD/C_ND | mb_sd+p k/n | C_SD+p | g_b (SD+p) | g_b (SD) | g_b (inkl. zmierzone) | SF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3-5 | 3245/20141 | 0,161 ± 0,003 | 38/294 | 0,129 ± 0,020 | 0,80 ± 0,12 | 7/82 | 0,085 ± 0,031 | 0,679 | 0,681 | 0,701 | 1,000 ± 0,100 |
| 5-8 | 1827/7923 | 0,231 ± 0,005 | 17/94 | 0,181 ± 0,040 | 0,78 ± 0,17 | 6/25 | 0,240 ± 0,085 | 0,250 | 0,237 | 0,228 | 1,000 ± 0,100 |
| 8-30 | 843/2734 | 0,308 ± 0,009 | 9/41 | 0,220 ± 0,065 | 0,71 ± 0,21 | 3/9 | 0,333 ± 0,157 | 0,070 | 0,081 | 0,071 | 1,000 ± 0,100 |

Udział prawdziwego trybu w mianowniku po decayType: mb_nd 0,939 / 0,939 / 0,939 (łącznie 0,939); mb_sd 0,951 / 0,949 / 0,976 (łącznie 0,953).

Prawdziwe Dplus (wszystkie mody, obszar fid.) w mb_sd: 4494, z protonem w oknie: 1210 (f_win = 26,9 %).

### X2. ρ łącznie (3–30 GeV, wagi g_b z MC SD z protonem)

| oszacowanie | ρ | uwagi |
|---|---|---|
| 1a: mb_sd wprost, wszystkie zdarzenia SD | 0,786 ± 0,092 | 64 zrekonstruowanych |
| 1b: mb_sd wprost, proton w oknie | 0,749 ± 0,169 | 16 zrekonstruowanych |
| 1c: mb_sd wprost, zliczenie łączne (ND przeważone do widma mianownika SD) | 0,783 ± 0,098 | |
| 2: krotność nSelBL, dane — okno sygnału (NOMINAŁ) | 0,982 ± 0,019 (bootstrap) | per bin: 0,978 / 1,017 / 0,910; kandydatów w oknie per bin [81, 15, 3] |
| 2: dane — pasma boczne | 0,974 | per bin: 0,966 / 1,007 / 0,927; kandydatów [227, 27, 13] |
| 2: dane — okno − pasma | 1,018 | per bin: 1,017 / 1,027 / 1,000; ujemne przedziały: tak (tylko opisowo) |
| domknięcie: metoda 2 na krotności sygnału reco mb_sd (nSelBL) | 0,943 | per bin: 0,906 / 1,017 / 0,939 |
| domknięcie: to samo, nTracksPt500Eta25 | 0,960 | per bin: 0,929 / 1,023 / 0,948 |
| niedomknięcie ρ₁(1a)/ρ₂(SD) − 1 | -16,6 % ± 9,7 % | |
| **przyjęte ρ** | **0,982 ± 0,165** | stat 0,019, warianty 0,008, niedomknięcie 0,163; ρ₁(1a) wobec ρ₂(dane): 2,1σ |

Mediana nSelBL: dane tagowane okno 21,0, pasma 20,0; sygnał reco mb_nd 29,0, mb_sd 19,0. Strony kandydatów (tylko A / tylko C / obie): okno [56, 42, 1], pasma [133, 134, 0].
Prawdziwe D w mb_sd z protonem: mediana nTracksPt500Eta25 = 19,0 (mb_nd: 33,0); udział zdarzeń z < 5 torami 1,8 %, ≤ 2 torami 0,2 %.

### X3. C_b(n) w mb_nd (nSelBL; przedziały o ≥ 100 zrekonstruowanych)

- 3-5 GeV: n ≥ 0: 0,108; n ≥ 8: 0,153; n ≥ 10: 0,156; n ≥ 12: 0,177; n ≥ 14: 0,149; n ≥ 16: 0,142; n ≥ 18: 0,152; n ≥ 20: 0,154; n ≥ 22: 0,146; n ≥ 24: 0,169; n ≥ 26: 0,166; n ≥ 28: 0,165; n ≥ 30: 0,160; n ≥ 32: 0,171; n ≥ 34: 0,156; n ≥ 36: 0,178; n ≥ 38: 0,161; n ≥ 40: 0,155; n ≥ 42: 0,200; n ≥ 44: 0,185; n ≥ 46: 0,169; n ≥ 49: 0,162; n ≥ 53: 0,175; n ≥ 58: 0,195
- 5-8 GeV: n ≥ 0: 0,193; n ≥ 10: 0,232; n ≥ 14: 0,245; n ≥ 17: 0,213; n ≥ 20: 0,272; n ≥ 23: 0,239; n ≥ 26: 0,222; n ≥ 29: 0,246; n ≥ 31: 0,181; n ≥ 34: 0,253; n ≥ 37: 0,231; n ≥ 40: 0,261; n ≥ 43: 0,219; n ≥ 47: 0,227; n ≥ 53: 0,241
- 8-30 GeV: n ≥ 0: 0,285; n ≥ 16: 0,271; n ≥ 23: 0,287; n ≥ 29: 0,322; n ≥ 34: 0,311; n ≥ 40: 0,350; n ≥ 47: 0,335

### X4-nom. Dplus: nominał (fit TAGGED) — OPISOWO (kontrolki D⁺ niespełnione)

N = 29,1 −10,8 +11,5 ± 8,8 (kształt); N_acc = 0,17; C̄_tag = 0,186 (z g_b inkluzywnym 0,184; udziały binów w C̄: 57,7 / 31,7 / 10,6 %); ⟨A·ε_AFP⟩ = 0,921; mianownik = 7,286 µb⁻¹.

**σ_tag^fid = 4,0 −1,5 +1,6 (stat) ± 1,5 (syst) µb** (syst 37,0 %; bez kształtu fitu tagowanego 21,1 %); L 7 % zamiast 2,2 %: syst 37,6 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,55 −0,20 +0,22 (stat) ± 0,20 (syst) %** (σ_incl = 725,7 µb; syst 35,8 %, bez kształtu 18,8 %); przedział 68 % [0,26; 0,84] %, 95 % [-0,01; 1,12] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | 30,4 | 30,4 |
| N_acc | 0,3 | 0,3 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 8,8 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 6,7 | – |
| ρ = C_D(tag)/C_D(ND) | 16,8 | 16,8 |
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

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,13–1,32) → zgodny (+0.7σ od centralnej); twarda R = 0,68 % (0,60–1,43) → zgodny (-0.5σ). σ_pred = R × σ_incl zmierzone: 2,5 (1,0–9,6) µb; × σ_incl Pythii A3 (980,0 µb): 3,4 µb; twarda: 4,9 (4,4–10,3) µb.

### X4-var. Dplus: wariant okno − pasma — OPISOWO (kontrolki D⁺ niespełnione)

N = 31,7 −12,1 +12,1; N_acc = 0,17; C̄_tag = 0,186 (z g_b inkluzywnym 0,184; udziały binów w C̄: 57,7 / 31,7 / 10,6 %); ⟨A·ε_AFP⟩ = 0,921; mianownik = 7,286 µb⁻¹.

**σ_tag^fid = 4,3 −1,7 +1,7 (stat) ± 0,9 (syst) µb** (syst 21,1 %; bez kształtu fitu tagowanego 21,1 %); L 7 % zamiast 2,2 %: syst 22,1 %.

**R^fid = σ_tag^fid/σ_incl(3–30) = 0,60 −0,23 +0,23 (stat) ± 0,11 (syst) %** (σ_incl = 725,7 µb; syst 18,8 %, bez kształtu 18,8 %); przedział 68 % [0,34; 0,85] %, 95 % [0,10; 1,09] %.

| pozycja budżetu | σ_tag^fid [%] | R^fid [%] |
|---|---|---|
| kształt fitu tagowanego (TAGGED) | – | – |
| N_acc | 0,3 | 0,3 |
| obciążenie fitu (domknięcie inkluzywne, liniowo) | 8,8 | – |
| SF klasyfikatora (3–5 GeV: scenariusz) | 6,7 | – |
| ρ = C_D(tag)/C_D(ND) | 16,8 | 16,8 |
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

Pythia (poziom cząstek, 3–30 GeV): miękka A3 R = 0,35 % (0,13–1,32) → zgodny (+1.0σ od centralnej); twarda R = 0,68 % (0,60–1,43) → zgodny (-0.3σ). σ_pred = R × σ_incl zmierzone: 2,5 (1,0–9,6) µb; × σ_incl Pythii A3 (980,0 µb): 3,4 µb; twarda: 4,9 (4,4–10,3) µb.

Przeliczenia do literatury z prawdy mb_sd (Dplus): N(pT ≥ 3)/N(pT ≥ 3,5) przy |η| < 2,5 = 1,367; N(|η| < 2,2)/N(|η| < 2,5) przy pT ≥ 3 = 0,891.
