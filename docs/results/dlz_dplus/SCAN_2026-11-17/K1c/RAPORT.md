# K1c — zysk H+DLZ dla D⁺ na kandydatach: sparowany bootstrap, skan WP 50–90, pseudo-dane z kandydatami MC (17/18.11.2026, projekt B)

**Werdykt wg `../KRYTERIA_K1c.md` (zapisane przed liczeniem): NIEROZSTRZYGNIĘTY.**
427929: średni zysk dla rzędu z F-testu +8,7 ± 1,9 % (< 10 %), ale przy wspólnym rzędzie pol3/pol3 +10,7 ± 2,0 % (≥ 10 %) — wniosek zależy od rzędu tła (warunek pomocniczy 4 niespełniony).
428770: zysk +17,1 ± 1,9 % (PASS, także pol3/pol3 i pol4/pol4), ale obciążenie z pseudo-danych poza 5 % w 3–5 GeV → wg kryteriów NIEPOTWIERDZONY.
Pull width w normie wszędzie. Test obciążenia okazał się nierozstrzygalny tą metodą (§4) — to jest główne ograniczenie, nie statystyka.

Tabele: `tables.md`; post hoc `fakegap.md`; maszynowo `verdict.json`, `out/k1c.json`, `out/prep.json`; logi `logs/`.

## 1. Warunek wstępny i wykonanie
- Selekcja i progi jak `analyze.py`. Odtworzenie `res_RUN`: liczby kandydatów i zdarzeń identyczne (2 006 962 / 1 037 458 i 7 167 801 / 2 417 158), progi = `thr.json` (Δ = 0),
  histogramy identyczne, N_S i σ(N) przy WP70 identyczne (ΔN = 0,000 %) w 24 fitach (H i H+DLZ, 3 biny, 2 runy, fit swobodny i z kształtem ustalonym) — `logs/prep.log`.
- `fastfit.py` (ten sam NLL, szybsza implementacja) vs `fitlib.fit`: |ΔN| ≤ 0,03 %, |Δσ| ≤ 0,4 %.
- 200 sparowanych replik po zdarzeniach (wspólne wagi dla H, H+DLZ i wszystkich WP), 200 pseudo-eksperymentów; czas 11 min zegarowych (limit 45).
- Kształt sygnału: 3–5 GeV ustalony na MC, 5–8 i 8–30 swobodny (żaden wybrany fit swobodny nie stanął na granicy).

## 2. Zysk przy WP70, definicja K1 (1 − rel_DLZ/rel_H), %; błąd ze sparowanego bootstrapu
| run | pT | pol2/pol2 | F-test (rząd H/DLZ) | pol3/pol3 | pol4/pol4 | N_DLZ/N_H (pol3/pol3) |
|---|---|---|---|---|---|---|
| 427929 | 3–5 | +20,6 ± 3,3 | +12,7 ± 3,9 (2/3) | +19,2 ± 3,6 | +17,6 ± 4,1 | 1,03 ± 0,05 |
| 427929 | 5–8 | +9,1 ± 2,9 | +9,1 ± 2,9 (2/2) | +8,3 ± 3,3 | +6,0 ± 4,1 | 1,08 ± 0,05 |
| 427929 | 8–30 | +4,2 ± 2,8 | +4,2 ± 2,8 (2/2) | +4,4 ± 3,6 | +3,0 ± 4,3 | 0,96 ± 0,04 |
| 427929 | średnio | +11,3 ± 1,7 | **+8,7 ± 1,9** | **+10,7 ± 2,0** | +8,9 ± 2,4 | |
| 428770 | 3–5 | +32,6 ± 3,4 | +23,5 ± 4,1 (2/3) | +29,7 ± 4,0 | +26,1 ± 4,6 | 1,15 ± 0,07 |
| 428770 | 5–8 | +10,9 ± 3,1 | +22,3 ± 3,2 (4/2) | +9,3 ± 3,5 | +7,2 ± 4,3 | 1,09 ± 0,05 |
| 428770 | 8–30 | +6,4 ± 1,4 | +5,4 ± 2,6 (4/4) | +5,8 ± 1,7 | +5,4 ± 2,6 | 0,98 ± 0,02 |
| 428770 | średnio | +16,7 ± 1,6 | **+17,1 ± 1,9** | **+14,9 ± 1,9** | +12,9 ± 2,3 | |

- pol2/pol2 odtwarza K1 (+20,7/+9,1/+4,2 i +32,9/+10,9/+6,4; różnice ≤ 0,3 pkt z innego ziarna replik). Błędy sparowane 2–4 pkt (w K1b, bez parowania, 5–10).
- Rząd z F-testu jest asymetryczny w 3–5 GeV (H pol2, H+DLZ pol3) i w 428770 5–8 (H pol4, H+DLZ pol2): w pierwszym przypadku zaniża zysk H+DLZ o ~7 pkt, w drugim zawyża o ~12 pkt.
  Przy wspólnym rzędzie zysk maleje z rzędem łagodnie i tak samo w obu runach: 427929 11,3 → 10,7 → 8,9 %; 428770 16,7 → 14,9 → 12,9 %.
- Zysk jest istotny tylko w 3–5 GeV (18–20 % i 26–33 %, 4–9σ); 5–8 GeV: 6–11 % (1,5–3,5σ); 8–30 GeV: 3–6 % (≤ 1,5σ w 427929, 2–4,6σ w 428770).
- σ z bootstrapu zdarzeniowego / σ hesjanu = 0,89–1,08 → błąd fitu jest wiarygodny w danych dla pol2–pol4.

## 3. Skan WP 50/60/70/80/90 i stosunek R = (N/ε)_DLZ/(N/ε)_H
Surowy R (pol3/pol3; pełne tabele z F-testem i pol4/pol4: `tables.md` §3–4):
| run | pT | WP50 | WP60 | WP70 | WP80 | WP90 | ⟨R⟩ | ⟨R_corr⟩ (V1) |
|---|---|---|---|---|---|---|---|---|
| 427929 | 3–5 | 0,98 ± 0,03 | 0,99 ± 0,04 | 1,03 ± 0,05 | 1,18 ± 0,08 | 1,21 ± 0,10 | 1,016 ± 0,026 | 0,98 ± 0,04 |
| 427929 | 5–8 | 0,92 ± 0,04 | 1,04 ± 0,04 | 1,08 ± 0,05 | 1,00 ± 0,05 | 1,20 ± 0,08 | 0,999 ± 0,026 | 0,97 ± 0,04 |
| 427929 | 8–30 | 1,02 ± 0,05 | 1,01 ± 0,04 | 0,96 ± 0,04 | 1,01 ± 0,05 | 1,06 ± 0,09 | 0,992 ± 0,024 | 1,02 ± 0,03 |
| 428770 | 3–5 | 1,03 ± 0,04 | 1,09 ± 0,05 | 1,15 ± 0,07 | 1,38 ± 0,11 | 1,84 ± 0,26 | 1,067 ± 0,034 | 1,12 ± 0,06 |
| 428770 | 5–8 | 1,10 ± 0,05 | 1,02 ± 0,05 | 1,09 ± 0,06 | 1,01 ± 0,06 | 1,32 ± 0,17 | 1,065 ± 0,032 | 1,05 ± 0,05 |
| 428770 | 8–30 | 1,01 ± 0,02 | 1,01 ± 0,02 | 0,98 ± 0,02 | 1,01 ± 0,03 | 1,10 ± 0,04 | 1,021 ± 0,013 | 0,99 ± 0,02 |

- 427929 (czysty minbias): z tłem pol3/pol4 R jest zgodny z 1 przy WP ≤ 70 we wszystkich binach (już bez korekty), średnio 0,99–1,02 ± 0,025. Nadmiar zostaje przy WP80–90 w 3–5 GeV (1,18–1,21, ~2,3σ) i WP90 w 5–8.
- 428770: w 3–5 GeV R rośnie z WP także dla pol4/pol4 (1,04 → 1,63) i po korekcie (χ²(R_corr = 1) = 11–20 na 5 punktów) — nierozwiązane; 5–8 GeV R ≈ 1,05–1,09 (1–2σ); 8–30 zgodne z 1 poza WP90.
- Płaskość N/ε: dla H+DLZ przejście pol2 → F-test/pol4 zmniejsza χ²_flat/4 z 45 do 15 (427929 3–5) i z 93 do 11–13 (428770 3–5); w 428770 8–30 ani H (χ² 14–18), ani H+DLZ (22–47) nie są płaskie.

## 4. Pseudo-dane (kandydaci MC z prawdą + tło n_all(m)·ε_B(m) z danych)
Obciążenie b = ⟨N_fit⟩/N_true i szerokość pulla przy WP70, rząd z F-testu:
| run | pT | b H / H+DLZ, tło V1 (pasma > 45 MeV) | b H / H+DLZ, tło V2 (> 70 MeV) | pull width H / H+DLZ (V1) | zysk skor. V1 / V2 [%] |
|---|---|---|---|---|---|
| 427929 | 3–5 | 1,08 / 1,08 | 1,30 / 1,09 | 1,06 / 1,00 | +12,8 / +26,9 |
| 427929 | 5–8 | 0,88 / 1,01 | 0,91 / 0,94 | 1,06 / 1,09 | −3,8 / +6,1 |
| 427929 | 8–30 | 0,92 / 0,95 | 1,05 / 1,13 | 1,03 / 1,19 | +1,2 / −3,5 |
| 428770 | 3–5 | 1,22 / 1,17 | 1,48 / 1,38 | 1,05 / 0,99 | +26,5 / +28,8 |
| 428770 | 5–8 | 0,98 / 1,02 | 1,01 / 1,02 | 1,02 / 0,99 | +19,2 / +21,7 |
| 428770 | 8–30 | 0,98 / 0,99 | 1,13 / 1,15 | 0,94 / 0,97 | +4,0 / +3,5 |

- **Pull width 0,89–1,19 we wszystkich komórkach, rzędach i obu wariantach tła (WP70)** → σ(N) z fitu pokrywa rozrzut; kryterium 3 PASS.
- **Kryterium 2 (|b − 1| < 5 %) nie przechodzi**: V1 — wszystkie 3 biny 427929 i 3–5 GeV w 428770; V2 — jeszcze gorzej. Ale b zmienia się między V1 i V2 o 0,03–0,26, czyli wynik zależy od
  wariantu interpolacji bardziej niż od fitu. Błędy statystyczne b: 0,002–0,006.
- Średni zysk skorygowany (F-test): 427929 +3,4 % (V1) / +9,9 % (V2); 428770 +16,6 % / +18,0 %. Dla pol3/pol3: 427929 +8,2 / +12,5 %; 428770 +13,5 / +16,3 %.
- Korekta V1 psuje płaskość H (427929 3–5: χ²_flat/4 1,5 → 16,9; 428770 3–5: 2,2 → 55), czego dane nie potwierdzają — kolejny znak, że b z tej konstrukcji nie jest obciążeniem fitu.

**POST HOC (po obejrzeniu wyników; nie wchodzi do werdyktu) — domknięcie samej interpolacji** (`fakegap.md`, `scripts/fakegap.py`): tę samą interpolację ε_B zastosowałem do „fałszywej przerwy”
bez sygnału (m_c = M0 + 130 MeV) i porównałem przewidziane tło z zaobserwowanym. Przy WP70 błąd interpolacji tła to od −2,5 % do +8,6 % (V1) i od −15 % do +15 % (V2),
co przeliczone na N_S pod pikiem daje pozorne „obciążenie” od −3 % do +37 % (V1) i od −40 % do +47 % (V2); sama precyzja statystyczna tego domknięcia to 1,5–7 % N_S.
Wniosek: przy S/B = 0,06–1 tło pod pikiem trzeba znać na 0,3–3 %, żeby wyznaczyć obciążenie N_S na 5 %; interpolacja z pasm bocznych tego nie daje w żadnym binie.
Kryterium 2 w obecnej postaci jest więc nietestowalne — FAIL w tabeli oznacza „nie wykazano”, nie „fit jest obciążony”.

## 5. Kryteria
| | 427929 | 428770 |
|---|---|---|
| 1. średni zysk (def. K1), rząd F-test ≥ 10 % | **FAIL**: +8,7 ± 1,9 % (−0,7σ od progu) | PASS: +17,1 ± 1,9 % (+3,7σ) |
| 2. \|b − 1\| < 5 %, V1 | **FAIL** (3–5: 1,08; 5–8: H 0,88; 8–30: H 0,92) | **FAIL** (3–5: 1,22 / 1,17) |
| 3. pull width 0,8–1,2 | PASS | PASS |
| 4. pol3/pol3 ten sam znak wniosku co F-test | **NIE** (+10,7 ± 2,0 % ≥ 10 %) | tak (+14,9 ± 1,9 %) |
| wynik | NIEROZSTRZYGNIĘTY | NIEPOTWIERDZONY (przez kryt. 2) |

## 6. Wnioski
1. Zysk statystyczny H+DLZ jest realny i nie znika przy elastyczniejszym tle, ale w czystym minbiasie wynosi 9–11 % (średnio po binach) — dokładnie na progu 10 %; cały pochodzi z 3–5 GeV (+18–20 %).
   W 428770 13–17 %.
2. Z tłem pol3 N_DLZ/N_H przy WP70 w 427929 jest zgodne z 1 (1,03 / 1,08 / 0,96, błędy 0,04–0,05); deficyt K3 był artefaktem exp(pol2). Otwarte: wzrost R z WP w 3–5 GeV w 428770.
3. Obciążenia fitu na poziomie 5 % nie da się wykazać ani wykluczyć pseudo-danymi opartymi o pasma boczne (żadną z trzech prób: pol4, pol6, ε_B·n_all). Potrzebne inne narzędzie
   (np. MC minbias z pełnym tłem przy S/B danych, albo kontrola na D⁰/D*-tag, gdzie tło pod pikiem jest mierzalne).
4. Rząd tła trzeba ustalić wspólnie dla obu klasyfikatorów (pol3), bo wybór osobny przesuwa zysk o ±7–12 pkt.

## 7. Zastrzeżenia
- Kryterium 1 nie było ślepe: wartości centralne dla F-testu znałem z K1b (zapisane w KRYTERIA_K1c §0).
- Pseudo-dane nie zawierają składowych piko-podobnych tła; ε_B pod pikiem = gładka interpolacja.
- Plik kandydatów zawiera tylko kandydatów „luźnych” (WP95 dowolnego tagu) — wystarcza dla WP ≤ 90.
- Rząd „F-test” w replikach ustalony na danych nominalnych (niepewność wyboru rzędu nie wchodzi do błędu).
- 428770 to próbka jetowa przy MC minbias (ε^MC i kształty z minbias).

