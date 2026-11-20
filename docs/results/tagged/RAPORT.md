# TAGGED — mezony D w zdarzeniach z tagiem protonu AFP: pierwszy odczyt (projekt B), 19.11.2026

Fity pików D i uzyski w regionie z tagiem dopuszczone 19.11.2026 (wieczór). Kryteria zapisane **przed** ekstrakcją i odczytem: `KRYTERIA.md` (23:15 CET).
Kolejność: kryteria → ekstrakcja (lxplus) → przygotowanie tylko na pasmach (`out/prep.json`, 23:50) → **jeden odczyt** 23:51:09 CET (`out/UNBLINDED`, `out/unblind.json`).
Po odczycie nie zmieniono metody, cięć, binów ani modelu; wszystko, co nie wynika z kryteriów, jest w §7 oznaczone „post hoc”.
Wyniki lokalne; nota i podsumowania nietknięte. Pełne tabele: `tabele/wyniki.md`. Rysunki: `fig/`.

## 0. Wynik w skrócie

1. **Wynik główny — D⁰, run 427929, region MBTS, pT 3–30 GeV, WP70:** N_tag = **54,7 −18,7 +19,4 (stat) ± 12,6 (kształt)**; istotność lokalna Z = 3,05 (iloraz wiarygodności),
   p z 2000 zabawek = 0,0015 (2,97σ), bez systematyki. Wg kryterium (Z ≥ 3 z zabawek) to **„wskazanie” na granicy 3σ, nie obserwacja**.
   Ułamek tagowany: f_raw = 0,88 %, po odjęciu tagów przypadkowych (0,20 ± 0,10 zdarzenia) **f_true = 0,88 −0,30 +0,31 (stat) ± 0,22 (syst) %**; przedział 68 % [0,51; 1,26] %, 95 % [0,16; 1,63] %.
2. **Porównanie z przewidywaniami (D⁰):** z pasm 0,365 % (0,18–0,73) — **zgodny**; z Pythii 0,194 % (0,077–0,77) — **zgodny** (przedział 68 % zachodzi na oba zakresy).
   Pomiar leży powyżej obu wartości centralnych (+1,4σ i +1,9σ) i żadnej nie zawiera w przedziale 68 %; przewidywań nie rozróżnia.
   Odpowiednik w liczbie zdarzeń przy WP70: zmierzone 55 wobec 23 (11–46) z pasm i 12 (5–48) z Pythii.
3. **D⁺ (kontrolka K1), 427929 MBTS:** N_tag = 29,1 −10,8 +11,5 ± 8,8 (Z = 2,9); f_true = **0,50 −0,19 +0,20 ± 0,15 %**; zgodny z oboma przewidywaniami (0,237 % i 0,194 %; +1,1σ, +1,3σ).
   f(D⁺)/f(D⁰) = 0,56 ± 0,37 (−1,2σ od 1) — zgodne. Ułamek wspólny D⁰ + D⁺ (opisowo): 0,61 ± 0,21 %.
4. **428770:** R1 (`noalg_L1AFP_A_OR_C`, 0,469 nb⁻¹): D⁰ 0,9 −16,0 +16,7 ± 11,1 (Z = 0,1), D⁺ 24,8 −11,0 +11,6 ± 7,1 (Z = 2,4). R2 (`j20_L1AFP ∧ ¬R1`, trigger dżetowy): D⁰ 278 ± 49 ± 48 (Z = 6,0), D⁺ 185 ± 30 ± 24 (Z = 6,5).
   R3: D⁰ 19,9 −12,1 +12,8 ± 5,0, D⁺ −1,3 −6,1 +6,9 ± 4,7. Klasa J12 (ułamek): D⁰ 0,45 ± 0,17 ± 0,08 %, D⁺ 0,25 −0,13 +0,14 ± 0,05 %.
5. **Kontrolki:** stabilność WP — D⁰ spełniona (0,0 / 1,3 / 1,3 σ_Δ), **D⁺ WP50 niespełniona (3,6 σ_Δ)**; tag losowy z konstrukcji (K4) — D⁰ OK, **D⁺ ⟨pull⟩ = −0,29 (poza ±0,2)**; region przypadkowy z konstrukcji (K5) — brak nadwyżki;
   spójność runów (K6) — w 2σ; D* (K2) — bez mocy. **D⁰ przy pT ≥ 5 GeV: fit zdegenerowany (N_B → 0), wynik niemiarodajny** (§5, §7).
6. **Główne zastrzeżenie (model, nie statystyka):** w modelu D⁰ (M8) z sygnałem związane są szerokie składowe — odbicie K↔π (0,74 N_S) i DMM + CM (1,3–1,7 N_S); razem 3,4 N_S w zakresie fitu.
   N_tag(D⁰) zależy więc od założenia, że proporcje te są w zdarzeniach z tagiem takie jak inkluzywnie. Rachunek z samych zapisanych liczb (okno minus tło ekstrapolowane z pasm, §7) daje ~30 zamiast 55.

## 1. Definicje (bez zmian względem TAGREG / noty §6.3, §6.11)

- Tag: dokładnie jeden proton po stronie z 0,030 ≤ ξ ≤ 0,100; T2 = GRL ∧ (A ∨ C) ∧ `nPileupVertices == 0` ∧ ξ_cal,s ≤ ξ_AFP,s + 0,02 na każdej tagowanej stronie.
- 427929: bit `HLT_mb_mbts_L1MBTS_1`. 428770: R1, R2, R3 jak w nocie; klasa J12 = bit `HLT_noalg_L1J12` (dla ułamka).
- Kandydaci i klasyfikator: wiersze app etapu 4, drzewo H, WP70 (progi per bin pT z fd_nd, jak WPSCAN); pT 3–30 GeV jednym binem; pomocniczo 5–30 GeV oraz WP50 / 80 / 90.
- Fit: Poisson na binach 5 MeV, D⁰ 1,70–2,05, D⁺ 1,70–2,20; kształt sygnału stały z fitu inkluzywnego tego samego runu i klasy z wetem npu (D⁰ M8, D⁺ M0), suma binów pT z wagami uzysków inkluzywnych;
  wolne N_S, N_B, tło exp(Czebyszew n), n z samych pasm (n = 2 tylko D⁺ R1, R2, D⁺ MBAFP, D⁰ WP80, D⁰ R2 8–30). Przedziały z profilu, Z = √q₀.

## 2. Dane i złączenie

- Pliki score (`charm_q/dnet4/app/d150`) mają klucz `event_uid`, tabele tagu numer zdarzenia → nowa ekstrakcja: pakiet `Q_lxplus/tagged` (mapa uid → zdarzenie z `evtable` h5 v2.0).
  Wynik: `runs/lxplus/tg_cand_<run>.npz` (kandydaci w zdarzeniach GRL ∧ tag), `tg_incl_<run>.npz` (histogramy inkluzywne z wetem npu i bez).
- Kontrola złączenia (przed odczytem): histogramy inkluzywne bez weta npu = `wps_hist.npz` **co do zliczenia** w obu runach (Σ|Δ| = 0); bity `cls` zgodne z evtable (0 niezgodności);
  100 % zdarzeń tagowanych 428770 odnalezione w `bits_428770.parquet` (bit `j20_L1AFP_A_OR_C`). Zdarzenia T2 z kandydatem D⁰ poza oknem: 10 090 (427929), 74 018 (428770).
- Uzyski inkluzywne z wetem npu, WP70 (HESSE × 1,3): 427929 MBTS D⁰ 6208 ± 499 (4013 / 1529 / 666), D⁺ 5845 ± 321; 428770 „all” D⁰ 16 868 ± 922, D⁺ 13 237 ± 573; klasa J12 D⁰ 13 266 ± 796, D⁺ 10 589 ± 437.

## 3. Zabawki przed odczytem (tło z pasm tagowanych; `tabele/wyniki.md` T5)

| pozycja | N wstrz. | obciążenie | rms N̂ | rms pull | pokrycie 68 % | mediana Z |
|---|---|---|---|---|---|---|
| D⁰ 427929 MBTS | 0 / 10 / 20 / 40 (po 2000) | +1,0 / −0,3 / −0,6 / 0,0 (± 0,4) | 17,0 / 17,6 / 18,5 / 19,8 | 0,99–1,01 | 0,68 / 0,70 / 0,68 / 0,69 | 0,0 / 0,6 / 1,1 / 2,1 |
| D⁺ 427929 MBTS | 0 / 19 (po 300) | −1,1 ± 0,5 / +0,2 ± 0,6 | 9,2 / 10,6 | 1,04 / 1,03 | 0,68 / 0,67 | 0,0 / 2,0 |

- Obciążenie ≤ 0,06 σ, pokrycie 64–74 % we wszystkich pozycjach T2 poza: D⁺ 428770 J12 (0,61) i D⁺ 427929 MBAFP (0,62) → wg reguły błąd × 1,08 (zastosowane). W regionach K5 (kilka kandydatów) zabawki nie mają sensu (pokrycie 0–0,67).
- Przed odczytem było więc wiadomo: przy przewidywaniach 12–23 D⁰ oczekiwana istotność ~0,6–1,2σ; zmierzony uzysk jest większy.
- Uwaga (z góry znana wada zabawek): tło zabawek pochodzi z fitu samych pasm z N_S = 0, więc szerokie składowe związane z sygnałem są w nim wliczone w tło — zabawki nie testują degeneracji z §0 pkt 6.

## 4. Odczyt: uzyski i ułamki

| gat. | run / region | L [nb⁻¹] | kand. w zakresie / w oknie | N_tag | stat (68 %) | syst kształtu | Z lok. | f_raw [%] | f_true [%] (stat; syst) |
|---|---|---|---|---|---|---|---|---|---|
| **D⁰** | **427929 MBTS** | 0,4648 | 790 / 274 | **54,7** | −18,7 +19,4 | 12,6 | **3,05** (p_zab. 0,0015) | 0,881 | **0,878** (−0,301 +0,313; 0,216) |
| D⁺ | 427929 MBTS | 0,4648 | 366 / 100 | 29,1 | −10,8 +11,5 | 8,8 | 2,88 | 0,498 | 0,495 (−0,185 +0,197; 0,153) |
| D⁰ | 428770 R1 | 0,4687 | 627 / 192 | 0,9 | −16,0 +16,7 | 11,1 | 0,06 | – | – |
| D⁺ | 428770 R1 | 0,4687 | 331 / 94 | 24,8 | −11,0 +11,6 | 7,1 | 2,35 | – | – |
| D⁰ | 428770 R2 | 27,65 | 5077 / 1684 | 277,5 | −48,3 +49,1 | 48,1 | 5,96 | – | – |
| D⁺ | 428770 R2 | 27,65 | 2745 / 686 | 185,1 | −30,1 +30,8 | 23,8 | 6,50 | – | – |
| D⁰ | 428770 R3 | 29,41 | 340 / 113 | 19,9 | −12,1 +12,8 | 5,0 | 1,68 | – | – |
| D⁺ | 428770 R3 | 29,41 | 182 / 39 | −1,3 | −6,1 +6,9 | 4,7 | 0 | – | – |
| D⁰ | 428770 klasa J12 | – | 1125 / 365 | 60,6 | −22,5 +23,3 | 9,7 | 2,77 | 0,457 | 0,446 (−0,170 +0,175; 0,078) |
| D⁺ | 428770 klasa J12 | – | 637 / 155 | 27,4 | −13,8 +14,8 (× 1,08) | 5,4 | 2,20 | 0,259 | 0,249 (−0,130 +0,140; 0,052) |

- Tagi przypadkowe (wzór z kryteriów, N_incl · f_acc · P_tag,SB): 427929 D⁰ 0,20 ± 0,10, D⁺ 0,17 ± 0,08 — pomijalne; J12: 1,4 ± 0,4 i 1,1 ± 0,3. R1, R2: N_true = N_tag · (1 − 0,05 ± 0,03), orientacyjnie.
- Ułamek kandydatów pasm bocznych w T2 (P_tag,SB, WP70): D⁰ 0,372 %, D⁺ 0,328 % (427929); klasa J12: 0,223 % i 0,208 %.
- Budżet syst kształtu wyniku głównego (ΔN): rząd tła + 1: 9,2; zakres 1,75–2,00: 7,8; sc / dmu: 2,7; odbicie / κ: 2,3; wagi pT: 0,8. Do f dochodzi N_incl (8 %) i f_acc.
- **Co się nie skraca w f (niepoprawiane):** wydajność MBTS_1 / J12 dla zdarzeń z protonem względem inkluzywnych, akceptancja i wydajność AFP (założone 0,55; 0,35–0,80), wymóg dokładnie jednego protonu,
  wydajność cięcia ξ_cal (1,000 +0/−0,005 w 427929), różnica widma pT i krotności D tagowanych, ewentualna różnica wydajności klasyfikatora w zdarzeniach dyfrakcyjnych.
- **Przekroju tagowanego nie podaję** (akceptancja AFP niezmierzona). Wielkość widzialna N_tag/L, nieskorygowana: D⁰ 427929 118 ± 41 nb⁻¹, D⁺ 63 ± 24 nb⁻¹; R2: 10,0 ± 1,8 i 6,7 ± 1,1 nb⁻¹.
- R2 per bin pT (opisowo, tylko gdzie ≥ 2000 kandydatów w pasmach): D⁰ 3–5 GeV 90,8 −39,1 +39,8 (Z = 2,4); pozostałe biny i D⁺ poniżej progu — nie fitowane.

### Porównanie z przewidywaniami (kryterium: część wspólna przedziału 68 % / 95 % z zakresem założeń)

| | f_true [%], przedział 68 % | z pasm (S) | werdykt | z Pythii (M) | werdykt |
|---|---|---|---|---|---|
| D⁰ 427929 | 0,88 [0,51; 1,26] | 0,365 (0,18–0,73) | zgodny (+1,4σ od centralnej) | 0,194 (0,077–0,77) | zgodny (+1,9σ) |
| D⁺ 427929 | 0,50 [0,26; 0,75] | 0,237 (0,12–0,47) | zgodny (+1,1σ) | 0,194 (0,077–0,77) | zgodny (+1,3σ) |

Oba gatunki leżą powyżej wartości centralnych obu scenariuszy; górna połowa obwiedni Pythii (modele typu ABMST, R do 1,5 %) i górny brzeg scenariusza pasm są preferowane, ale bez mocy rozstrzygającej.
Pomiar nie wyklucza żadnego scenariusza na poziomie 95 %.

## 5. Kontrolki

| | kryterium | wynik | status |
|---|---|---|---|
| K1 D⁺ / D⁰ | f(D⁺)/f(D⁰) = 1 w 2σ | 0,56 ± 0,37 (−1,2σ) | spełnione |
| K2 D* (427929 MBTS; 428770 R1) | zgodność z f(D⁰) w 2σ | N_tag = 29 ± 18 (Z = 1,7) przy N_incl = 1990 ± 600 → f_raw = 1,5 ± 1,0 %; R1: 46 ± 18 (Z = 2,6) | bez mocy; zgodne w 2σ; patrz niżej |
| K3 stabilność WP, D⁰ | \|Δf\| ≤ 2 σ_Δ | WP50 0,88 %, WP80 0,51 %, WP90 1,20 % wobec 0,88 %: 0,0 / 1,3 / 1,3 σ_Δ | spełnione |
| K3 stabilność WP, D⁺ | \|Δf\| ≤ 2 σ_Δ | WP50 0,77 %, WP80 0,41 %, WP90 0,62 % wobec 0,50 %: **3,6** / 0,8 / 0,5 σ_Δ | **WP50 niespełnione** — zastrzeżenie, nominał zostaje |
| K4 tag losowy (200 podpróbek inkluzywnych o liczności próbki tagowanej) | \|⟨pull⟩\| ≤ 0,2; pokrycie ≥ 0,60 | D⁰: +0,16; 0,65 (obciążenie +3,3 ± 1,4 przy 23,8). D⁺: **−0,29**; 0,64 (−2,7 ± 0,8 przy 20,0) | D⁰ spełnione; **D⁺ niespełnione** (tło rzędu 1 za sztywne dla kształtu inkluzywnego na 1,70–2,20) |
| K5 tag przypadkowy z konstrukcji | brak nadwyżki ponad N_incl · P_CR,SB | 427929: D⁰ 0,1 [−1,6; 2,6] (oczek. 0,3), D⁺ 0 kandydatów; 428770 J12: D⁰ −2,7 [−6,2; 1,3] (oczek. 1,1), D⁺ 1,8 [−0,5; 4,9] (0,8) | brak nadwyżki; moc mała także w 428770 (22–36 kandydatów) |
| K6 spójność runów (N_tag/L, bit `noalg_L1AFP`) | zgodność w 2σ | D⁰: 427929 71 ± 34 wobec R1 2 ± 35 nb⁻¹ (+1,4σ); D⁺: 17 ± 24 wobec 53 ± 24 (−1,1σ) | spełnione; D⁰ w R1 nie widać |
| K7 strony (D⁰ 427929) | opisowo | tylko A 17,5 −13,1 +13,9; tylko C 36,3 −13,0 +13,8 | – |
| pT ≥ 5 GeV (427929) | opisowo | D⁰ 34,4 −8,0 +5,1 (Z = 5,0; f = 1,57 %) — **fit na granicy N_B → 0** (N_B = 5 przy 114 kandydatach); D⁺ 5,6 −4,0 +4,8 (f = 0,24 %) | D⁰ niemiarodajne; D⁰ i D⁺ rozbieżne |

- **D\*:** kandydaci slim nie mają klasyfikatora; S/B inkluzywnie ~0,3 %, kształt sygnału z fitu inkluzywnego źle określony (parametry na ograniczeniach, N_incl ± 30 %). Kontrolka nie ma mocy — zgodnie z oczekiwaniem z kryteriów, ale słabiej niż zakładałem.
- Stabilność w 428770 i K4 dla regionów 428770 nie były w kryteriach — nie liczone.

## 6. Co wysłano na lxplus, stan pętli

- Pakiet `charm_work/Q_lxplus/tagged` (`tg_extract.py`, `tg_job.sh`, `tg.sub`, `submit_tg.sh`, `tg_common.sh`), etap `ORCH/stages/tagged.sh`, wpis `tagged` na końcu `pipeline.yaml` (kopia `pipeline.yaml.bak_tagged`).
  Pętla ORCH zrestartowana 23:20 (stara 155762 → nowa 191419, `MAX_PARALLEL=8`); jedna pętla, `dnet6` pozostał `running`, jego joby nietknięte.
- Job `charm_tagged` (klaster 12793757, 1 rdzeń, 4 GB) wysłany 23:20, czekał w kolejce ~20 min; w tym czasie ten sam kod uruchomiony interaktywnie (`/tmp/klyczak_tagged`, 4 min 54 s + 1 min 51 s) → `runs/interactive/`.
  Job zakończył się 23:43, etap `tagged` = `done`, wynik w `runs/lxplus/` — **tablice identyczne** z przebiegiem interaktywnym. Analiza liczona na `runs/interactive/`.
- Wyniki na EOS: `/eos/user/k/klyczak/charm_q/tagged/`. Na AFS tylko logi (`Q_lxplus/logs/tagged`).
- **Pozostałości (nie usuwałem):** lxplus `/tmp/klyczak_tagged/` (skrypt, `o/` z testu 3 plików, `full/` ~21 MB); lokalnie `runs/interactive/` (duplikat `runs/lxplus/`), `ORCH/pipeline.yaml.bak_tagged`, katalog testowy w katalogu tymczasowym.

## 7. Uwagi po odczycie (post hoc — metody nie zmieniałem)

1. **Degeneracja szerokich składowych D⁰ z tłem.** W M8 z sygnałem związane jest 0,74 N_S odbicia i κ = 1,26–1,66 N_S składowej DMM + CM, obie szerokie na cały zakres. W próbce tagowanej (790 kandydatów)
   fit może przerzucać zdarzenia między tłem wykładniczym a tymi składowymi. Arytmetyka z liczb zapisanych przed i przy odczycie (bez nowego fitu): okno 274 − tło ekstrapolowane z samych pasm 227,9 = 46 zdarzeń
   sygnałopodobnych w oknie, czyli N_S ≈ 46 / 1,54 = **30** (f ≈ 0,48 %) wobec 54,7 z fitu, który obniża tło w oknie do 188. Dla D⁺ (bez szerokich składowych) ten sam rachunek daje 100 − 71,4 = 28,6 → N_S ≈ 32, zgodnie z fitem 29,1.
   Rachunek z pasm jest z kolei zaniżony (traktuje szerokie składowe w pasmach jako tło). Wniosek: **N_tag(D⁰) leży między ~30 a ~55 zależnie od traktowania szerokich składowych; przypisana systematyka 12,6 tego nie pokrywa w pełni.**
   Istotność liczona z okna i pasm byłaby zbliżona (~2,5–3σ), ułamek D⁰ zbliżyłby się do D⁺ (0,5 %).
   *Propozycja (niewykonana):* fit z κ i f_refl jako wolną wspólną normalizacją szerokich składowych albo fit jednoczesny próbki tagowanej i inkluzywnej; wybór przed kolejnym odczytem, na zabawkach z szerokimi składowymi w prawdzie.
2. **D⁰ pT ≥ 5 GeV:** N_B → 5 przy 114 kandydatach, wszystko poza pikiem opisane składowymi związanymi z sygnałem; przedział −8,0 +5,1 i Z = 5,0 są artefaktem granicy. Liczby nie należy cytować. Wąski pik w oknie jest widoczny (rys. 1), ale jego uzysk wymaga modelu z pkt 1.
3. **D⁺ WP50 (K3):** N_tag nie spada przy zaostrzeniu (30,1 przy WP50 wobec 29,1 przy WP70; inkluzywnie 3928 wobec 5845). Możliwe: fluktuacja tła przy WP70 pod pikiem, większa wydajność klasyfikatora dla D w zdarzeniach z tagiem
   (mniejsza krotność) albo niedoszacowanie przy WP70 związane z K4 (−14 %). σ_Δ dla próbek zagnieżdżonych jest przybliżeniem. Nierozstrzygalne tą statystyką; f(D⁺) przy WP50 = 0,76 ± 0,21 % jest bliżej D⁰.
4. **K4 D⁺:** tło exp rzędu 1 na 1,70–2,20 zaniża sygnał o 2,7 ± 0,8 zdarzenia, gdy tło ma kształt inkluzywny. Tło tagowane może mieć inny kształt (wybór rzędu z pasm dał n = 1, dq = 1,6). *Propozycja:* dla D⁺ zakres 1,75–2,10 jako nominał albo rząd 2 — do ustalenia z góry.
5. **R1 wobec 427929:** D⁰ w R1 = 0,9 ± 16 przy tej samej świetlności co 427929 (tam z bitem `noalg_L1AFP` 33 ± 15), D⁺ odwrotnie (25 ± 11 wobec 8 ± 11). Każde w 1,1–1,4σ; razem nie wskazują na problem, ale R1 nie potwierdza samodzielnie sygnału D⁰.
   W drzewie 428770 jest 44 % zdarzeń `noalg_L1AFP` (LUMI §4) — dla zdarzeń z D niezweryfikowane.
6. **Strony:** C ma więcej niż A (36 wobec 18), choć wydajność L1 AFP jest wyższa po stronie A — tag offline nie wymaga triggera AFP, więc nie musi to być sprzeczne; 1σ.
7. **Istotność:** Z = 3,05 nie zawiera systematyki modelu z pkt 1 ani poprawki na 2 gatunki i ~30 fitów pomocniczych; liczy się tylko wynik główny. D⁰ i D⁺ razem (niezależne próbki, Z = 3,0 i 2,9) to wyraźny sygnał D w zdarzeniach z tagiem;
   kombinacji formalnej nie było w kryteriach — nie liczę.
8. **Odstępstwa od kryteriów (wszystkie przed odczytem):** (a) D*: bin Δm 0,5 MeV w 0,13945–0,15995 zamiast 0,25 MeV w 0,139–0,160 — Δm w tabelach `afp_ext` jest skwantowane co 0,1 MeV (aliasing), ostatni bin niepełny; okno 0,14345–0,14745;
   tło inkluzywne D* z dodatkowym czynnikiem exp(bx + cx²). (b) MBAFP: kształt i N_incl z klasy MB (klasa inkluzywna z bitem AFP ma za mało zdarzeń). (c) Zabawki wyniku głównego: 2000 na punkt zamiast 1000.
   (d) K5 w 428770 liczone w klasie J12 (kryteria nie podawały klasy). (e) Przebieg interaktywny obok etapu ORCH (identyczny wynik). (f) Fity per bin pT w R2 tylko D⁰ 3–5 GeV (próg 2000 kandydatów pasm).
   (g) Systematyka kształtu policzona dla pozycji głównych, regionów i K6; dla WP pomocniczych, K5, K7 nie.

## 8. Czego nie dało się zrobić

- Przekrój tagowany bezwzględny (brak akceptancji AFP), ułamek w R1 / R2 (trigger wymaga protonu — brak mianownika inkluzywnego), poprawka na wydajność triggera dla zdarzeń dyfrakcyjnych.
- Miarodajna kontrolka D* (brak score dla D* w plikach app; bez klasyfikatora S/B ~0,3 %).
- Rozstrzygnięcie między przewidywaniem z pasm a Pythią; pomiar różniczkowy.
- f_acc po cięciu klasyfikatora (nota §6.3.3: „must be recomputed”) — użyto wartości bez score; przy 0,2 zdarzenia w 427929 bez znaczenia, w 428770 orientacyjne.

## 9. Pliki (`~/Documents/charm/charm_work/TAGGED_2026-11-19/`)

- `KRYTERIA.md` (przed odczytem), `RAPORT.md`.
- `skrypty/`: `tg_lib.py` (selekcje, kształt, fit, profil, zabawki), `t1_prep.py` (tylko pasma; `out/prep.json`, `out/prep.log`, `out/toys_null_q0.npz`), `t1b_dst.py` (D*: `prep` / `unblind`),
  `t2_unblind.py` (jeden odczyt; `out/unblind.json`, `out/unblind.log`, `out/hists.npz`, `out/UNBLINDED`), `t3_fig.py`, `t4_tables.py`. Ekstrakcja: `../Q_lxplus/tagged/`.
- `tabele/wyniki.md` (T1 uzyski, T2 ułamki, T3 przewidywania, T4 systematyka, T5 zabawki i K4, T6 D*, T7 kontrolki).
- `fig/fig1_427929_MBTS.png` (wynik główny, D⁺, warianty 5–30), `fig2_428770_regiony_D0.png`, `fig3_428770_regiony_Dplus.png`, `fig4_kontrolki_CR_MBAFP.png`, `fig5_Dstar.png`.
- `runs/lxplus/` (wynik etapu ORCH), `runs/interactive/` (identyczny).
