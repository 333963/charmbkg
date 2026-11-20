# R3 — kalibracja kształtu score H w binach (D⁰←D*, fit Δm jednoczesny w 6 kategoriach) — 18.11.2026, projekt B

Kryteria i projekt zapisane przed fitami danych: `KRYTERIA.md` (§1–4 przed wszystkim; §5 uzupełnienia po teście silnika na pseudo-danych,
nadal przed pierwszym fitem z sygnałem na danych). Pełne tabele: `tables.md`; liczby: `verdict.json`; rysunki Δm: `fig_dm_428770.png`, `fig_dm_427929.png`.

## 1. Co było lokalnie
Wszystko potrzebne — nic nie pobierano, nic na lxplus. `B_trig/runs/lxplus/bt_{data_428770,data_427929,fd_nd,mb_nd}.npz` (h5 v2.0, L_xy ≥ 0,25 mm):
pary D⁰←D* per kandydat z `score` (H nominalny), `dm`, `m`, `pt`, `run`, `event`, bity triggera/GRL, prawda w MC. Pary: 1 233 580 (428770),
402 470 (427929), 178 721 (fd_nd). Progi WP i kształty sygnału: `B_trig/runs/lxplus/out/btrig_sf.json`. `FIT_scorebins/out`, `DLZ_d0`, `DNET/runs`
nie były potrzebne (FIT_scorebins = D⁰ inkluzywne bez tagu, tylko pasma boczne).

## 2. Metoda (skrypty)
| plik | rola |
|---|---|
| `r3_common.py` | silnik: fit Poissona Δm (0,1 MeV, 140,5–160) w kategoriach score; sygnał 2 Gaussy z MC (wspólne μ, skala s), tło u^a(1+c1t+c2t²) per kategoria, gradient analityczny; wydajki S_c → f_c (Σ = 1), E(WP), SF, kowariancja |
| `r3_mc_width.py` → `mc_width.json` | reguła szerokości per bin score na MC (wymagana tylko w 3–5 GeV) |
| `r3_run.py toys/data` → `out/` | test wstrzyknięcia (300 toyów × 12), fit nominalny + 200 replik bootstrapu sparowanego po zdarzeniach, pass/fail tym samym silnikiem przy WP95/80/70/50, wariant równych ułamków |
| `r3_report.py` → `tables.md`, `verdict.json`, rysunki | tabele i kryteria |
| `r3_posthoc_count.py` → `posthoc_count.json` | kontrola post hoc wydajków zliczeniem w oknie |

Kategorie: 0 = odrzucone (< WP95), potem WP95–80, 80–70, 70–50, 50–25, > WP25 (ułamki MC 5/16/10/20/24/25 %) — granice = oficjalne WP, żeby
porównanie z pass/fail było przy tym samym cięciu (odstępstwo od „równej liczby sygnału” zapisane z góry; wariant równy policzony jako kontrola, tab. D6).
Kontrola silnika: pass/fail z 2 kategoriami odtwarza oficjalny `fit_passfail` we wszystkich 18 punktach (|Δ| ≤ 0,08σ; tab. D3).

## 3. Test wstrzyknięcia (kryterium a) — SF(WP70), 300 toyów, tło z pasm bocznych Δm danych per kategoria, N_S oficjalne
| run | pT | wstrzyknięte | odzyskane ⟨SF⟩ (obciążenie) | σ toya | pull HESSE (odporny, post hoc) | pass/fail: obciążenie | (a) |
|---|---|---|---|---|---|---|---|
| 428770 | 5–8 | 0,843 / 1,000 | 0,872 (+3,4 ± 0,8 %) / 1,023 (+2,3 ± 0,7 %) | 0,113 / 0,123 | 1,05 / 1,02 | +3,4 / +2,1 % | **NIE** (3,4 > 3 %) / tak |
| 428770 | 8–30 | 0,845 / 1,000 | 0,849 (+0,4 %) / 1,012 (+1,2 %) | 0,041 / 0,045 | 1,12 / 1,09 | +0,2 / +1,1 % | tak / tak |
| 427929 | 5–8 | 0,843 / 1,000 | 0,867 (+2,8 %) / 1,008 (+0,8 %) | 0,101 / 0,122 | 1,02 / **2,01** (0,98) | +3,2 / +0,7 % | tak / **NIE** (pull) |
| 427929 | 8–30 | 0,845 / 1,000 | 0,851 (+0,6 %) / 1,007 (+0,7 %) | 0,074 / 0,085 | 1,06 / **2,53** (1,15) | +0,4 / +0,7 % | tak / **NIE** (pull) |
| 427929 | 3–5 | 0,843 / 1,000 | 1,01 (+20 ± 9 %) / 1,15 (+15 ± 5 %) | 1,30 / 0,77 | 1,09 / 1,11 | +12 / +11 % | **NIE** |
| 428770 | 3–5 | — | σ toya 2,8–3,1 (niemierzalny) | | | | pominięty |

- Fity dobre (valid ∧ kowariancja dodatnio określona): 97–100 %. Obciążenie w σ: ≤ 0,27σ wszędzie.
- Obciążenie +2…+3,4 % w 5–8 GeV jest **takie samo w pass/fail** — pochodzi z modelu fitu (kształt piku 2 Gaussy vs rzeczywiste Δm z MC przy S/B ~ 2 %),
  nie z binowania; rośnie ku ciasnym WP (WP50: +3,5…+5,1 %, WP25: +4,8…+6,5 % w 428770 5–8).
- Szerokość pulla 2,0 / 2,5 w 427929 przy c = 0 to 1–10 toyów z 300 z zaniżonym błędem HESSE (min 0,004 przy medianie 0,08); odporna szerokość
  (połowa przedziału 16–84 %) 0,98–1,15 — **post hoc**, kryterium w zapisanej postaci niespełnione. Wniosek praktyczny: błąd brać z bootstrapu, nie z HESSE.

## 4. Dane: SF kształt vs pass/fail (błąd = bootstrap sparowany po zdarzeniach, 200 replik; statystyka MC ≤ 0,011, wspólna)
| run | pT | WP | SF kształt | SF pass/fail (ten silnik) | oficjalny | σ_k/σ_pf | kształt − pf |
|---|---|---|---|---|---|---|---|
| 428770 | 5–8 | 95 | 0,836 ± 0,099 | 0,839 ± 0,097 | 0,839 ± 0,099 | 1,02 | −0,003 ± 0,017 |
| 428770 | 5–8 | 80 | 0,864 ± 0,124 | 0,865 ± 0,123 | 0,865 ± 0,115 | 1,01 | −0,001 ± 0,014 |
| 428770 | 5–8 | **70** | **0,909 ± 0,133** | **0,913 ± 0,134** | – | **1,00** | −0,004 ± 0,014 (0,3σ) |
| 428770 | 5–8 | 50 | 0,954 ± 0,139 | 0,957 ± 0,139 | 0,957 ± 0,135 | 1,00 | −0,003 ± 0,023 |
| 428770 | 8–30 | 95 | 1,008 ± 0,035 | 1,006 ± 0,035 | 1,006 ± 0,036 | 1,00 | +0,002 ± 0,002 |
| 428770 | 8–30 | 80 | 0,981 ± 0,040 | 0,974 ± 0,041 | 0,974 ± 0,041 | 0,98 | +0,007 ± 0,004 |
| 428770 | 8–30 | **70** | **0,991 ± 0,041** | **0,986 ± 0,043** | – | **0,97** | +0,006 ± 0,006 (1,0σ) |
| 428770 | 8–30 | 50 | 0,947 ± 0,045 | 0,942 ± 0,047 | 0,938 ± 0,048 | 0,96 | +0,006 ± 0,011 |
| 427929 | 5–8 | 95 | 1,045 ± 0,111 | 1,053 ± 0,112 | 1,053 ± 0,106 | 0,99 | −0,009 ± 0,022 |
| 427929 | 5–8 | 80 | 1,032 ± 0,123 | 1,032 ± 0,128 | 1,031 ± 0,122 | 0,97 | 0,000 ± 0,022 |
| 427929 | 5–8 | **70** | **0,995 ± 0,123** | **0,995 ± 0,131** | – | **0,94** | 0,000 ± 0,029 (0,0σ) |
| 427929 | 5–8 | 50 | 1,003 ± 0,125 | 0,985 ± 0,133 | 0,987 ± 0,134 | 0,94 | +0,018 ± 0,049 |
| 427929 | 8–30 | 95 | 0,913 ± 0,069 | 0,913 ± 0,070 | 0,913 ± 0,066 | 0,99 | +0,001 ± 0,006 |
| 427929 | 8–30 | 80 | 0,945 ± 0,087 | 0,953 ± 0,098 | 0,951 ± 0,087 | 0,89 | −0,008 ± 0,030 |
| 427929 | 8–30 | **70** | **0,983 ± 0,097** | **1,002 ± 0,106** | – | **0,91** | −0,019 ± 0,035 (0,6σ) |
| 427929 | 8–30 | 50 | 0,965 ± 0,110 | 1,009 ± 0,133 | 1,010 ± 0,127 | 0,83 | −0,044 ± 0,066 |

- Nowe: SF(WP25) = 0,99 ± 0,16 / **0,72 ± 0,05** (428770 5–8 / 8–30), 0,99 ± 0,15 / **0,74 ± 0,13** (427929).
- Korelacja SF kształt–pass/fail na replikach: 0,95–0,99 (WP ≤ 80) — obie metody mierzą to samo z tych samych zdarzeń.
- **3–5 GeV**: 428770 N_tot = 823 ± 423 (1,9σ) — pik widać tylko w dwóch górnych kategoriach (rys.), poniżej WP50 wydajki zgodne z zerem lub ujemne
  (S_0 = −296 ± 359) → niemierzalny, pominięty. 427929: N_tot = 788 ± 253 (3,1σ, formalnie „mierzalny” wg reguły), ale SF70 = 1,37 ± 0,77 (bootstrap;
  HESSE 0,42 — zaniżony ×1,8), f_0 = −31 %, toy obciążony +15…+20 % → bez wartości pomiarowej. W obu runach mianownik (sygnał bez cięcia) nie jest wyznaczony.
- Błąd bootstrapowy jest większy od HESSE o 0–15 % (5–8 GeV: 0,133 vs 0,116), także dla pass/fail.

### Macierz korelacji SF między WP (brakujący element; bootstrap ⊕ MC; pełne kowariancje w `verdict.json` → `SF_cov_tot`)
| | 428770 5–8 | 428770 8–30 | 427929 5–8 | 427929 8–30 |
|---|---|---|---|---|
| ρ(95, 80) | 0,78 | 0,81 | 0,87 | 0,81 |
| ρ(95, 70) | 0,75 | 0,76 | 0,84 | 0,77 |
| ρ(95, 50) | 0,72 | 0,68 | 0,78 | 0,70 |
| ρ(80, 70) | 0,97 | 0,94 | 0,96 | 0,95 |
| ρ(80, 50) | 0,93 | 0,82 | 0,90 | 0,84 |
| ρ(70, 50) | 0,97 | 0,86 | 0,93 | 0,87 |
| ρ(50, 25) | 0,94 | 0,72 | 0,86 | 0,68 |

## 5. Werdykt wg kryteriów
| run | pT | (a) odzysk | (b) σ_k/σ_pf ≤ 0,8 | (c) zgodność ≤ 2σ | werdykt |
|---|---|---|---|---|---|
| 428770 | 5–8 | NIE (c = 1: +3,4 %) | NIE (1,00) | tak | NIEUŻYTECZNA |
| 428770 | 8–30 | tak | NIE (0,97) | tak | UŻYTECZNA CZĘŚCIOWO |
| 427929 | 5–8 | NIE (pull HESSE przy c = 0) | NIE (0,94) | tak | NIEUŻYTECZNA |
| 427929 | 8–30 | NIE (pull HESSE przy c = 0) | NIE (0,91) | tak | NIEUŻYTECZNA |
| 427929 | 3–5 | NIE | NIE (1,25) | tak | NIEUŻYTECZNA |
| 428770 | 3–5 | – | – | – | pominięty |

**Werdykt łączny wg zapisanej reguły: NIEUŻYTECZNA** jako dokładniejsza miara SF. Sedno: (b) nie jest spełnione nigdzie i nie o włos —
fit w 6 kategoriach daje ten sam błąd co pass/fail (0,91–1,02×), bo ogranicza go znajomość liczby sygnału w brudnych kategoriach (< WP80),
a ta jest wspólna dla obu metod. Niepowodzenia (a) są graniczne (3,4 % vs 3 %, identyczne w pass/fail) lub techniczne (ogony błędu HESSE);
metoda nie jest *gorsza* od pass/fail, tylko nie lepsza. Co z niej zostaje: (i) pełna kowariancja SF między WP z bootstrapu, (ii) SF przy
dowolnej granicy z jednego fitu, (iii) lokalizacja niezgodności w score (niżej), (iv) lekki silnik (fit ~1 s), zgodny z `fit_passfail` co do 0,08σ.
Dokładność 10–13 % w 5–8 GeV to granica statystyczna próbki z tagiem D*, nie metody.

## 6. Deficyt 5–8 GeV w 428770 — **post hoc** (wg KRYTERIA §4 tylko przy zdanym (a); tu (a) niezdane o 0,4 pkt proc., obciążenie +3 % ≪ efekt 17 %)
R_c = f_c(dane)/f_c(MC), błąd bootstrap ⊕ MC:

| run, pT | < WP95 | WP95–80 | WP80–70 | WP70–50 | WP50–25 | > WP25 | χ²/5 (p) |
|---|---|---|---|---|---|---|---|
| 428770, 5–8 | 3,9 ± 1,8 | 0,70 ± 0,40 | 0,56 ± 0,24 | 0,80 ± 0,16 | 0,91 ± 0,14 | 0,99 ± 0,16 | 7,2 (0,21) |
| 427929, 5–8 | 0,2 ± 2,0 | 1,11 ± 0,32 | 1,28 ± 0,28 | 0,98 ± 0,18 | 1,02 ± 0,15 | 0,99 ± 0,14 | 1,4 (0,92) |
| 428770, 8–30 | 0,85 ± 0,63 | 1,14 ± 0,13 | 0,91 ± 0,11 | 1,10 ± 0,08 | 1,17 ± 0,07 | **0,72 ± 0,05** | 44,9 (< 10⁻⁷) |
| 427929, 8–30 | 2,6 ± 1,3 | 0,76 ± 0,26 | 0,68 ± 0,23 | 1,03 ± 0,17 | 1,19 ± 0,17 | 0,74 ± 0,12 | 9,9 (0,08) |

- **5–8 GeV, 428770: deficyt nie jest istotny w kształcie** (p = 0,21; żadna kategoria > 2σ). Wzór jest jednak spójny: górne 50 % sygnału zgodne z MC
  (0,91–0,99), braki w kategoriach pośrednich WP95–50 (0,56–0,80), nadmiar w odrzuconej (256 ± 137 zamiast ~65 par). Czyli „0,83 przy WP95” to
  sygnał widziany **poniżej WP95**, w kategorii z S/B ≲ 1,5 % w oknie ±2,5 MeV — tam fit tła decyduje o wyniku. Wariant równych ułamków: to samo (p = 0,14).
  427929 nie pokazuje tego wzoru (p = 0,92), ale przy błędach ±2,0 w kat. 0 nie wyklucza go.
- **Znalezisko uboczne, 8–30 GeV:** w obu runach za mało sygnału przy najwyższym score (> WP25: 0,72 ± 0,05 i 0,74 ± 0,12), nadmiar w WP50–25.
  Potwierdzone zliczeniem w oknie ±3 MeV bez kształtu sygnału (`posthoc_count.json`: 0,72) i w wariancie równych ułamków (0,67 ± 0,05).
  Pass/fail przy WP95/80/70 tego nie widzi (SF ≈ 1); przy WP50 było 0,94 ± 0,05. Przyczyna niezbadana (rozkład pT wewnątrz 8–30, modelowanie
  cech przy score → 1); ma znaczenie tylko dla punktów pracy ciaśniejszych niż WP50.

## 7. Zastrzeżenia
- Generator tła w toyach jest z tej samej rodziny funkcji co fit → test nie mierzy systematyki modelu tła, która w kategorii odrzuconej dominuje.
- Kształt sygnału: 2 Gaussy z MC ze wspólną skalą; obciążenie +2…+3 % (5–8 GeV) z tego modelu nie jest korygowane w tabelach.
- MC bez wag emulacji triggera (dla klasy „all” oficjalne ε_MC z emulacją i bez są równe); klasy triggera (J12only, MBTSbased) nie liczone.
- 427929: passGRL; kompletność kontenera niepotwierdzona (jak w B_trig).
- Uzupełnienia §5 KRYTERIA (złożoność tła wg statystyki kategorii, obsługa parametrów na granicy) wprowadzono po teście na pseudo-danych,
  przed fitami danych — nie post hoc względem danych, ale po pierwotnym zapisie.
- Post hoc: odporna szerokość pulla (§3), cała sekcja 6 dla 428770 5–8, kontrola zliczeniem.
- 300 toyów (plan 200), 200 replik bootstrapu; liczenie lokalne ~24 min zegarowe (8 rdzeni). Błąd błędu bootstrapowego ~5 %, więc stosunki σ_k/σ_pf
  znane do ok. ±0,03 (szacunek, nie liczone).
- Liczby z literatury (LHCb, CMS, ATLAS w RECON_lit_B3) nie były używane.
