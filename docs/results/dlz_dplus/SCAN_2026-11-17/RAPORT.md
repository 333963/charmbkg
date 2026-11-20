# H+DLZ dla D⁺→Kππ: skan cięcia (R1) i obciążenie fitu z pseudo-danych (R2) — raport, 17.11.2026 (projekt B)

**Werdykt wg `KRYTERIA.md` (zapisane 22:25, przed liczeniem): 427929 — NIEROZSTRZYGNIĘTE, 428770 — FAIL (płaskość w 8–30 GeV). H+DLZ nie spełnia warunku przyjęcia; H zostaje nominałem.**
Wynik fizyczny jest jednak inny niż w K3: po korekcie obciążenia stosunek (N/ε)_DLZ/(N/ε)_H przy WP70 jest zgodny z 1 we wszystkich 6 komórkach
(średnio 0,990 ± 0,015 i 1,010 ± 0,012), czyli deficyt K3 to w większości artefakt fitu, a nie złe ε_S. Nie przechodzi płaskość w WP (3–5 GeV w obu runach, 8–30 w 428770).

Tabele: `tables.md` (kryteria z góry), `tables_posthoc.md` (wariant post hoc, §5); maszynowo `verdict.json`, `verdict_posthoc.json`, `out/scan*.json`.

## 1. Zakres faktycznie wykonany (odstępstwa od planu, stwierdzone przed liczeniem — KRYTERIA §1)
- Danych per kandydat lokalnie nie ma (slim npz ze score, OOF MC i modele analizy są tylko na EOS; lokalne `oof/base.npz` uszkodzone, `h5/` bez score).
  Niczego nie pobierałem (reguła „nic na lxplus”). Liczone z `DATA_2026-11-17/runs/lxplus/ana/{res,rep}_RUN.*`: histogramy masy po cięciu, kształty i ε^MC, 200 sparowanych replik bootstrapu.
- **Brak WP60** (nie ma histogramu) → skan WP 50/70/80/90 (3 st. swobody), WP95 opisowo.
- **R2: sygnał parametryczny** (2 Gaussy dopasowane do MC z prawdą po cięciu, N_true = S0·ε^MC, S0 = N_H(WP70)/0,7), nie prawdziwi kandydaci MC → obciążenie od ogonów sygnału spoza modelu nie jest mierzone.
- **R2: tło** = exp(pol4) dopasowane do pasm bocznych danych po cięciu (A: |Δm| > 45 MeV, C: > 70 MeV), interpolowane pod pik; S/B jak w danych danego runu, binu i WP.
  Pominięte: wszystko, co ma strukturę w ±45 MeV (refleksje D_s⁺→KKπ, D⁺→πππ/KKπ, resztkowy D*⁺, lokalna zmiana ε_B(m) pod pikiem). R2 jest więc dolnym ograniczeniem obciążenia.
- 100 pseudo-eksperymentów na komórkę, fit identyczny (`fitlib`). Fit: swobodny w 5–8 i 8–30; w 3–5 GeV kształt ustalony na MC (reguła granicy; fit swobodny H+DLZ stoi na granicy przy WP80–95 w obu runach, 0–7 % replik udanych — pokazane w `tables.md` §3, nie maskowane).
- **Czas CPU przekroczony**: bieg nominalny 7 min zegarowych × 6 procesów, bieg post hoc 21 min × 6 (obciążony komputer) — łącznie ponad zadane ~40 min CPU. Liczby toyów nie zmniejszałem.

## 2. R1 — surowy skan (dane), wariant fitu komórki
| run | pT | χ²_flat/3 H | χ²_flat/3 DLZ | nachylenie H [1/ε] | nachylenie DLZ | R = (N/ε)_DLZ/(N/ε)_H: WP50 / 70 / 80 / 90 |
|---|---|---|---|---|---|---|
| 427929 | 3–5 | 0,4 | 36,8 | +0,06 ± 0,19 | +0,73 ± 0,14 | 1,01 ± 0,03 / 1,08 ± 0,04 / 1,28 ± 0,07 / 1,48 ± 0,14 |
| 427929 | 5–8 | 7,0 | 11,1 | −0,12 ± 0,18 | +0,53 ± 0,16 | 0,94 ± 0,03 / 1,10 ± 0,04 / 1,04 ± 0,05 / 1,25 ± 0,09 |
| 427929 | 8–30 | 4,6 | 6,0 | +0,17 ± 0,20 | +0,16 ± 0,19 | 1,05 ± 0,05 / 1,00 ± 0,04 / 1,05 ± 0,05 / 1,16 ± 0,07 |
| 428770 | 3–5 | 0,9 | 108,0 | −0,13 ± 0,26 | +1,15 ± 0,12 | 1,07 ± 0,04 / 1,25 ± 0,07 / 1,61 ± 0,12 / 2,35 ± 0,27 |
| 428770 | 5–8 | 4,5 | 19,1 | +0,26 ± 0,21 | +0,52 ± 0,17 | 1,14 ± 0,05 / 1,17 ± 0,05 / 1,14 ± 0,07 / 1,59 ± 0,16 |
| 428770 | 8–30 | 15,7 | 35,9 | +0,35 ± 0,09 | +0,47 ± 0,09 | 1,02 ± 0,02 / 1,03 ± 0,02 / 1,06 ± 0,03 / 1,17 ± 0,04 |

H jest płaskie (poza 428770 8–30, gdzie i H rośnie z WP: χ² 15,7/3 — w tym runie 82 % L1_J12, MC = minbias). H+DLZ rośnie ku luźnym WP wszędzie, najmocniej w 3–5 GeV.
Błąd różnicy selekcji zagnieżdżonych: bootstrap zgadza się z √|σ²_a − σ²_b| w ~5 % (np. 427929 5–8 H: 150 vs 140; nieskorelowane dałoby 181) — `tables.md` §4.

## 3. R2 — obciążenie b = ⟨N_fit/N_true⟩ (tło A; w nawiasie tło C)
| run | pT | H: WP50 / 70 / 80 / 90 | H+DLZ: WP50 / 70 / 80 / 90 |
|---|---|---|---|
| 427929 | 3–5 (fix) | 1,02 / 1,01 (1,01) / 0,99 / 0,97 | 1,05 / **1,07** (1,05) / 1,10 / **1,35** |
| 427929 | 5–8 | 0,99 / 0,99 (0,97) / 0,98 / 1,02 | 1,00 / 1,03 (0,99) / 1,04 / 1,10 |
| 427929 | 8–30 | 0,98 / 1,00 (0,97) / 1,00 / 0,98 | 1,06 / 1,05 (1,01) / 1,06 / 1,15 |
| 428770 | 3–5 (fix) | 0,99 / 1,01 (1,00) / 1,01 / 1,09 | 1,02 / **1,14** (1,10) / 1,29 / **1,73** |
| 428770 | 5–8 | 0,95 / 0,93 (0,92) / 0,98 / 1,01 | 1,00 / 1,07 (1,04) / 1,13 / 1,38 |
| 428770 | 8–30 | 1,00 / 1,00 (0,98) / 0,99 / 0,98 | 1,03 / 1,05 (1,04) / 1,04 / 1,07 |

Błędy statystyczne b: 0,3–1,5 %. Kontrola maszynerii (tło prawdy = exp(pol2) z fitu danych): b = 0,99–1,02. Rozrzut pojedynczego pseudo-eksperymentu przy WP70: 3–8 %.
Odczyt: tło po cięciu H+DLZ ma w pasmach bocznych krzywiznę, której exp(pol2) nie opisuje (χ²/ndf SB do 2–8 w 3–5 GeV przy luźnych WP; exp(pol4) daje 1,0–1,7),
i fit oddaje ją sygnałowi: zawyżenie N rośnie z WP dokładnie tak jak surowy R1. Dla H obciążenie jest ≤ 3 % (wyjątek 428770 5–8: −7 %).
Niezależne potwierdzenie kierunku: istniejące domknięcie na MC mb z prawdziwymi kandydatami (`tables.md` §7): H+DLZ N_fit/N_truth przy WP70 = 1,10 ± 0,06 / 1,13 ± 0,05 / 0,98 ± 0,06, przy WP90 1,87 / 1,41 / 1,06.

## 4. R1 skorygowane o R2 i kryteria (zapisane z góry)
| run | pT | (c) czuła | R_corr(WP70) | (a1) | ⟨R_corr⟩ po WP | Δχ²_flat (DLZ − H) surowe → skor. | (b) |
|---|---|---|---|---|---|---|---|
| 427929 | 3–5 | NIE (b_DLZ = 1,067) | 1,023 ± 0,044 | zgodne | 1,016 ± 0,027 | +36,4 → +11,4 | FAIL (komórka nieczuła) |
| 427929 | 5–8 | tak | 1,056 ± 0,044 | PASS | 0,993 ± 0,026 | +4,0 → −3,8 | PASS |
| 427929 | 8–30 | NIE (b_DLZ = 1,053) | 0,947 ± 0,038 | zgodne | 0,970 ± 0,024 | +1,3 → −0,8 | PASS |
| 428770 | 3–5 | NIE (1,137) | 1,109 ± 0,066 | zgodne (1,7σ) | 1,066 ± 0,037 | +107,1 → +7,4 | FAIL (nieczuła) |
| 428770 | 5–8 | NIE (H 0,926; DLZ 1,072) | 1,012 ± 0,047 | zgodne | 1,040 ± 0,033 | +14,6 → −0,4 | PASS |
| 428770 | 8–30 | tak | 0,978 ± 0,021 | PASS | 0,998 ± 0,013 | +20,3 → +5,7 | **FAIL** |

- (a1): |R_corr(WP70) − 1| < 2σ we wszystkich 6 komórkach (największe +1,7σ w 428770 3–5); błędy 4–7 % (2 % w 428770 8–30).
- (a2): ⟨R_corr⟩ po WP i binach: 427929 **0,990 ± 0,015**, 428770 **1,010 ± 0,012** → PASS w obu.
- (b): po korekcie nie przechodzi 3–5 GeV w obu runach (resztkowy wzrost z WP: R_corr(WP80) = 1,15 ± 0,07 i 1,26 ± 0,11; WP90 w 428770 1,48 ± 0,18) oraz 428770 8–30 (+5,7; tam samo H nie jest płaskie, χ² 18,6/3 w wariancie fix, 13,5 w free).
- (c): obciążenie H+DLZ przy WP70 ≥ 5 % w 4 z 6 komórek → test czuły tylko w 427929 5–8 i 428770 8–30.
- Werdykt: 427929 NIEROZSTRZYGNIĘTE (brak FAIL w komórce czułej, ale 2 z 3 nieczułe); 428770 FAIL ((b) w 8–30).

## 5. POST HOC (po obejrzeniu wyników; nie zmienia werdyktu z §4)
Po zobaczeniu tabel zauważyłem wadę konstrukcji tła prawdy: pasma boczne |Δm| > 45 MeV zawierają ~7 % sygnału (szeroki Gauss), więc exp(pol4) dopasowane do nich
wciąga ogony sygnału do „tła” i zawyża tło pod pikiem — efekt silniejszy przy małym tle, czyli różny dla H i H+DLZ. Wariant poprawiony: exp(pol4) dopasowane do SB
z odjętym oczekiwanym sygnałem S0·ε·kształt MC (`scan.py --tailsub`, `tables_posthoc.md`). Skutek: wszystkie b spadają o 3–5 % (H: 0,90–1,00; H+DLZ przy WP70: 1,04 / 0,99 / 1,00 w 427929, 1,13 / 1,03 / 1,00 w 428770),
ale **R_corr się nie zmienia**: WP70 = 1,025 ± 0,046 / 1,057 ± 0,040 / 0,953 ± 0,038 (427929), 1,104 ± 0,075 / 1,025 ± 0,042 / 0,978 ± 0,018 (428770); średnie 0,996 ± 0,015 i 1,011 ± 0,011.
Zmienia się natomiast klasyfikacja (c): w 427929 wszystkie 3 biny stają się czułe → werdykt post hoc 427929 = FAIL przez (b) w 3–5 GeV (Δχ² +13,8); 428770 = nierozstrzygnięte (H wychodzi poza 5 %).
Wniosek z porównania: kryterium (c) na bezwzględnym b jest kruche (zależy od konstrukcji tła na poziomie 3–5 %), stosunek DLZ/H jest stabilny.

## 6. Wnioski
1. Deficyt K3 (ΔSF_rel 0,53–0,88) i nadmiar N_DLZ/N_H są odtwarzane przez samo niedopasowanie exp(pol2) do gładkiego tła po cięciu H+DLZ. Po korekcie ε_S(H+DLZ)/ε_S(H) przy WP70 = 1 w granicach 4–7 % per komórka i 1,0–1,5 % średnio. Brak przesłanek, że cechy DLZ są źle modelowane w sygnale.
2. H+DLZ z obecnym fitem (exp(pol2)) nie nadaje się do pomiaru: zawyża N_S o 7–14 % w 3–5 GeV i 3–7 % wyżej już przy WP70, a przy luźnych WP w 3–5 GeV fit nie działa (granica, χ²/ndf do 6,3, obciążenie ×1,3–1,7, resztkowa niepłaskość po korekcie). Część zysku K1 (+11/+17 %) to większy dopasowany N_S, nie mniejszy błąd.
3. Kryteria przyjęcia nie są spełnione. Warunek powrotu: fit z tłem exp(pol3/pol4) (lub szablonem z SB) dla H i H+DLZ, powtórzone K1 i ten skan; w 3–5 GeV nie luźniej niż WP70.

## 7. Zastrzeżenia
- R2 nie widzi składowych piko-podobnych ani ogonów sygnału spoza 2 Gaussów; zgodność R_corr z 1 zakłada, że jedyne obciążenie to gładkie tło. Resztkowy wzrost z WP w 3–5 GeV może być właśnie taką składową.
- S0 wzięte z fitu H przy WP70 (obciążenie H ≤ 3 %, w 428770 5–8 −7 %) — wpływa na S/B prawdy w tym stopniu.
- Korekta jest duża względem progu: w 428770 3–5 R surowe 1,25 → 1,11 opiera się na obciążeniu 1,14 z modelu tła, którego pod pikiem nie da się sprawdzić.
- 428770 to próbka jetowa przy MC minbias; niepłaskość samego H w 8–30 GeV wskazuje na efekt ε(WP) niezwiązany z DLZ.
- Kowariancje z 200 replik; w komórkach „free” użyto replik udanych dla wszystkich 8 punktów (121–200).
- Pełna wersja testu (WP60, wstrzykiwanie prawdziwych kandydatów MC, bootstrap po kandydatach) wymaga pobrania `slim/data_*.npz` i OOF MC z EOS `charm_q/dlz_data` — do decyzji.

## 8. Pliki (`DLZ_dplus/SCAN_2026-11-17/`)
`KRYTERIA.md`, `RAPORT.md`, `tables.md`, `tables_posthoc.md`, `verdict.json`, `verdict_posthoc.json`, `out/scan.json`, `out/scan_posthoc.json`, `scripts/scan.py`, `scripts/tables.py`, `logs/`.
W `DATA_2026-11-17/` niczego nie zmieniono.
