# Nowe cechy dla D⁰→Kπ (odpowiednik DLZ) — raport (17.11.2026, projekt B)

**Wniosek: dla D⁰ nie ma zysku ≥ 10 % bez rzeźbienia.** Najlepszy bezpieczny zestaw to H + parametry zderzenia córek
(„NOKIN” = grupy ip + time, w praktyce sama grupa ip): **+8…+10 % Zproxy na MC, +5…+9 % po przeniesieniu na SB danych**, przy WP70.
Zestawy z kinematyką córek dają na danych pozornie +20…+46 %, ale odtwarzają masę (R² = 0,80) — odrzucone.
Tylko MC i pasma boczne danych 428770; pik D⁰ w danych i pary D–D nie były oglądane.

Tabele: `tables.md`; rysunki: `fig/fig1_zproxy_d0.png`, `fig/fig2_sculpt_d0.png`; dane: `../DLZ_dplus/DATA_2026-11-17/runs/lxplus/d0/`.

## 1. Co zrobiono

- Liczone na lxplus (węzeł `d0` etapu ORCH `dlz_data`), h5 v1.2: fd_nd 24/24 plików, mb_nd 40/80, dane 428770 6/39.
  Protokół jak w `FEAT_charmnet` (study / transfer / sculpt / modelling) + własny `d0_sculpt.py`.
- Baza: H (18 cech). Okno 1,75–2,00 GeV, SR ±55 MeV. 5 foldów po zdarzeniach; Zproxy = ε_S/√ε_B przy WP z OOF sygnału MC per bin pT;
  bootstrap sparowany, 200 replik.
- Testowane grupy nowych cech (`scripts/d0_feats.py`):
  - **d0prod**: d0_K·d0_π;
  - **ip**: istotności d0 córek (K, π, min, max; σ(d0) przybliżone 0,012 ⊕ 0,075/(pT√sinθ) mm), χ²(d0), Δz0, d0prod;
  - **time**: pseudo-czas własny cτ = L_xy·m/pT;
  - **kin2**: cosθ*(K), asymetria pT córek;
  - **kinall**: kin2 + min/max pT córek, stosunki pT, ΔR;
  - **pidbb**: pulle Bethe-Blocha dE/dx (K jako K/π, π jako π/K), dE/dx i p pionu, liczba punktów;
  - złożenia: NOKIN = ip + time, NOKIN_pid, ALLnew; odniesienia: H_seed2 (szum ziarna), H+ALL_noiso (z FEAT).
- Nie testowano (brak w h5): DCA między torami, zmienne z drzewa Tracks. Wymagałyby nowej ekstrakcji — otwarte.

## 2. Zysk ΔZproxy/Z_H przy WP70 [%]

| zestaw | MC 3–5 | MC 5–8 | MC 8–30 | transfer na SB danych 3–5 / 5–8 / 8–30 | trening na SB danych 3–5 / 5–8 / 8–30 |
|---|---|---|---|---|---|
| H_seed2 (szum) | +0,3 ± 0,3 | −0,7 ± 0,7 | +1,8 ± 1,9 | 0,0 / −0,1 / −0,3 | −0,3 ± 0,3 / −0,4 ± 0,6 / −0,7 ± 0,8 |
| H+d0prod | +3,5 ± 0,3 | +1,7 ± 0,6 | +4,8 ± 1,9 | +3,4 / +3,2 / +2,4 | — |
| **H+ip** | +10,6 ± 0,5 | +8,1 ± 1,0 | +8,0 ± 2,5 | +9,5 ± 0,6 / +5,3 ± 1,3 / +5,7 ± 1,7 | — |
| H+time | −0,3 ± 0,2 | +0,7 ± 0,5 | +4,2 ± 2,1 | −0,5 / −0,2 / 0,0 | — |
| **H+NOKIN** | +10,2 ± 0,4 | +8,5 ± 1,0 | +9,4 ± 2,6 | +9,4 ± 0,5 / +5,4 ± 1,3 / +6,0 ± 1,6 | +8,3 ± 0,5 / +7,3 ± 1,0 / +4,1 ± 1,1 |
| H+NOKIN_pid | +11,0 ± 0,5 | +8,2 ± 0,9 | +8,3 ± 2,7 | +9,9 / +5,1 / +6,7 | — |
| H+pidbb | +0,2 ± 0,2 | −0,4 ± 0,5 | +2,2 ± 1,7 | 0,0 / −0,3 / +1,5 | — |
| H+kin2 | +1,0 ± 0,3 | +0,1 ± 0,6 | +2,4 ± 1,7 | +1,0 / 0,0 / +0,4 | — |
| H+kinall | +3,2 ± 0,3 | +0,7 ± 0,6 | −0,1 ± 1,9 | +7,2 / +3,4 / −3,3 | +37,0 ± 0,9 / +21,6 ± 1,4 / +5,1 ± 1,0 |
| H+ALLnew | +12,2 ± 0,5 | +8,9 ± 1,0 | +6,5 ± 2,3 | +13,3 ± 0,6 / +6,7 ± 1,4 / +3,1 ± 1,6 | +45,9 ± 1,0 / +23,6 ± 1,5 / +8,3 ± 1,2 |

Błędy transferu: bootstrap ε_B. Z_H (MC, WP70): 3,53 ± 0,02 / 3,92 ± 0,05 / 3,66 ± 0,13.

Wynik: cały zysk niesie grupa **ip** (istotności d0, χ²(d0), Δz0). Pseudo-czas, cosθ*, asymetria pędów i PID dE/dx
nie dają nic ponad szum ziarna (≤ +1 % w 3–5 GeV). Próg 10 % jest osiągnięty tylko na MC w 3–5 GeV (+10,2 ± 0,4 %);
po przeniesieniu na dane: +9,4 / +5,4 / +6,0 % — **poniżej 10 %**.

## 3. Rzeźbienie masy (każdy zestaw)

| zestaw | R WP70, MC truth ≥ 8 (SR) | R WP70, dane (pasmo) | nachylenie ε_B [/100 MeV] | fałszywy sygnał N_fake/σ (3–30) | R² regresji masy z cech |
|---|---|---|---|---|---|
| H | 0,953 ± 0,008 | 0,963 ± 0,012 | −0,26 | −1,4σ | 0,278 |
| H+ip | 0,973 ± 0,009 | 0,978 ± 0,013 | −0,17 | −1,5σ | 0,277 |
| H+NOKIN | 0,965 ± 0,009 | 0,979 ± 0,013 | −0,17 | −1,9σ | 0,277 |
| H+pidbb | 0,955 ± 0,008 | 0,966 ± 0,012 | −0,26 | −2,4σ | 0,278 |
| H+kin2 | 0,961 ± 0,009 | 0,964 ± 0,012 | −0,24 | −1,2σ | 0,326 |
| H+kinall | 1,056 ± 0,010 | 1,017 ± 0,014 | −0,16 | −1,0σ | **0,800** |
| H+ALLnew | 1,026 ± 0,010 | 1,019 ± 0,014 | −0,10 | −1,4σ | **0,800** |

R² masy z samych nowych grup: ip 0,013, time 0,001, pidbb 0,013, kin2 0,017, **kinall 0,655**.

- **H+ip / H+NOKIN: brak rzeźbienia** — R zgodne z H w 2σ i bliżej 1, nachylenie płytsze niż w H, R² masy bez zmian,
  fałszywy sygnał w fitach na tle MC zgodny z H (−1,9σ vs −1,4σ; ujemny znak jest już w samym H).
- **kinall / ALLnew: odrzucone.** Min/max pT i stosunki pT córek razem z cechami H pozwalają odtworzyć masę (R² 0,80), R przeskakuje
  powyżej 1 (1,056 ± 0,010 vs 0,953 ± 0,008), a rozjazd „MC +3 % vs trening na SB danych +37 %” to sygnatura uczenia się masy
  (tło treningowe pochodzi z SB, sygnał z SR). Fit fałszywego sygnału tego nie wyłapał — sam ten test nie wystarcza.
- kin2 (cosθ*, asymetria): R² rośnie 0,278 → 0,326 przy zysku ~1 % — nie warto.

## 4. Modelowanie (C2ST SB danych vs MC mb, 3 ziarna)

H 0,5811; H+ip 0,5865 (Δ +0,0054); H+NOKIN 0,5866 (**Δ +0,0055**, poniżej progu 0,01); H+pidbb +0,0014; H+kinall +0,0016; H+ALLnew +0,0070.
C2ST 1D nowych cech ≤ 0,54.

Dwie uwagi techniczne:
- Przesunięcie dE/dx MC o −0,023 (jak w CharmNet/FEAT) zastosowane do cechy obok surowych dE/dx z H daje C2ST 0,9986
  (kontrola `pidshift_ctrl`) — to artefakt (deterministyczny znacznik MC), nie złe modelowanie. Cechy PID liczone są z surowego dE/dx.
- Pulle Bethe-Blocha i liczba punktów mają duże KS (0,36 i 0,61) przy małym C2ST 1D (0,53 / 0,50) — różnice w wartościach
  dyskretnych/wypełniających; bez znaczenia, bo PID nic nie wnosi.

## 5. Odpowiedź na pytanie zadania

- Czy jest ≥ 10 % dla D⁰ bez rzeźbienia? **Nie.** Bezpieczne maksimum: H+NOKIN (≈ H+ip), +10,2 ± 0,4 / +8,5 ± 1,0 / +9,4 ± 2,6 % na MC,
  +9,4 ± 0,5 / +5,4 ± 1,3 / +6,0 ± 1,6 % po przeniesieniu na dane, +8,3 / +7,3 / +4,1 % przy treningu na SB danych.
- Zproxy to nie σ(N_S)/N_S: przy S/B ≪ 1 zysk błędu ≈ zysk Zproxy, więc realnie oczekiwać ~5–9 % (hipoteza do sprawdzenia fitem
  piku D⁰ na danych — nie robione tutaj).
- Lekcja z D⁺ (`../DLZ_dplus/DATA_2026-11-17/RAPORT.md`): H+DLZ przeszedł testy MC/SB, a padł na fitach danych (closure, SF).
  H+ip dla D⁰ wymaga więc tej samej procedury na danych (fit przy kilku WP, ΔSF_rel, closure na MC) przed użyciem;
  σ(d0) jest tu przybliżone wzorem — docelowo błąd z dopasowania toru.

## 6. Ograniczenia

- mb_nd 40/80 plików, dane 6/39 plików v1.2; bin 8–30 GeV ma błędy 2–3 %.
- DCA torów i zmienne z drzewa Tracks nietestowane (brak w h5).
- Wariant „trening na SB danych” liczony tylko dla H_seed2, NOKIN, kinall, ALLnew.

## 7. Pliki

- Skrypty: `scripts/{d0_feats,d0_run,d0_sculpt,report_d0}.py`, `scripts/run_d0.sh` (kopie w `Q_lxplus/dlz_data/d0/`).
- Pozostałości: `frames/` (~1 GB, lokalne ramki z małej próbki), `screen/`, `logs/`, puste `out/` (lokalny przesiew przerwany z powodu obciążenia).

## Wynik (17.11.2026, etap `dlz_data`, węzeł d0) — nowe cechy dla D⁰
Tabele: `tables.md`. ΔZproxy/Z_H przy WP70, pT 3–5 / 5–8 / 8–30 GeV.
- **H+ip (parametry zderzenia córek)**: MC +10,6 ± 0,5 / +8,1 ± 1,0 / +8,0 ± 2,5 %; przeniesienie na dane (ε_B z dalekich SB 428770) ×1,095 ± 0,006 / 1,053 ± 0,013 / 1,057 ± 0,017; ΔC2ST tła +0,005; rzeźbienie nie większe niż H (R 0,973 ± 0,009 vs 0,953 ± 0,008; R² regresji masy 0,277 vs 0,278).
- H+kinall (kinematyka rozpadu): +3,2 % na MC, ale cechy odtwarzają masę (R² 0,80 vs 0,28; R 1,056 ± 0,010) → odrzucone. Trening na SB danych z kinall daje +37 % w 3–5 GeV — to rzeźbienie, nie zysk.
- time, pidbb, kin2: zero (|Δ| ≤ 1 % w 3–8 GeV). d0prod: +3,5 / +1,7 / +4,8 %.
- Wszystko bez kinematyki (H+NOKIN): +10,2 / +8,5 / +9,4 % MC, ×1,094 / 1,054 / 1,060 na danych, test z tłem z SB danych +8,3 / +7,3 / +4,1 %.
- Szum treningu (H_seed2): ≤ 0,7 % w 3–8 GeV.
Wniosek: jedyna grupa warta dodania dla D⁰ to parametry zderzenia (ok. +8–10 %). SF z tagu D* dla H+ip jeszcze niezmierzony — bez tego nie przyjmuję.

## Uzupełnienie 17.11.2026 18:00 — SF z tagu D* dla drzewa H + parametry zderzenia (z etapu `dnet1`)
Zestaw HX z `DNET` to H + parametry zderzenia + cτ (= „NOKIN”), drzewo trenowane na SB danych jak nominalne H (`D0.xgb.HX.sb`). Źródło: `DNET/RAPORT.md` §10.
- Zysk Zproxy względem H (WP70, 3–5 / 5–8 / 8–30 GeV): MC 1,101 ± 0,002 / 1,074 ± 0,006 / 1,041 ± 0,013; dane (ε_B z SB 428770) 1,081 ± 0,006 / 1,058 ± 0,012 / 1,055 ± 0,015; R_sculpt 0,959 ± 0,007 (H: 0,942 ± 0,006).
- SF (WP70, 5–8 / 8–30 GeV), HX wobec H: 427929 1,06 ± 0,13 / 0,98 ± 0,09 wobec 1,01 ± 0,13 / 0,97 ± 0,10; 428770 all 0,87 ± 0,12 / 0,97 ± 0,04 wobec 0,88 ± 0,12 / 0,97 ± 0,04; J12only 0,94 ± 0,18 / 0,99 ± 0,05 wobec 0,93 ± 0,18 / 0,99 ± 0,05. Bin 3–5 GeV niemierzalny.
Wniosek: dodanie parametrów zderzenia nie zmienia SF (różnice ≤ 0,05 przy błędach 0,04–0,13, te same zdarzenia) i daje +6–8 % na danych. Zysk jest poniżej progu 10 %, więc to bezpieczne, ale małe ulepszenie; decyzja o zmianie nominału pozostaje otwarta.
