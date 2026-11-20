# K1b — zysk K1 dla D⁺ z elastyczniejszym tłem (17.11.2026, projekt B)

**Werdykt wg `../KRYTERIA_K1b.md` (zapisane 23:05, przed liczeniem): zysk NIEPOTWIERDZONY.**
427929: średni zysk (b) +9,7 ± 4,5 % (< 10 %) i obciążenie poza 5 % w 5–8 i 8–30 GeV; 428770: średni zysk (b) +21,9 ± 4,9 % (PASS), ale obciążenie H+DLZ 1,12 w 5–8 GeV.
Szerokości pulli w normie we wszystkich komórkach. Zysk utrzymuje się wyraźnie tylko w 3–5 GeV.

Tabele: `tables.md`; maszynowo `verdict.json`, `out/k1b.json`; skrypty `scripts/{k1b,report}.py`.

## 1. Wynik przy WP70 (g = [σ/N]_H/[σ/N]_DLZ − 1, %)
| run | pT | (a) pol2/pol2 | rząd F-test H / DLZ | (b) wybrany | (c) pol4/pol4 | obciążenie H / DLZ (b) | pull w H / DLZ | (b) skor. o obciążenie | N_DLZ/N_H (b) |
|---|---|---|---|---|---|---|---|---|---|
| 427929 | 3–5 | +26,0 ± 4,6 | pol2 / pol3 | +14,6 ± 8,2 | +21,3 ± 9,9 | 0,97 / 1,02 | 1,01 / 0,94 | +8,8 | 1,02 ± 0,08 |
| 427929 | 5–8 | +10,0 ± 3,6 | pol2 / pol2 | +10,0 ± 7,1 | +6,4 ± 9,2 | **0,87** / 0,97 | 0,95 / 1,02 | −1,6 | 1,10 ± 0,08 |
| 427929 | 8–30 | +4,4 ± 3,2 | pol2 / pol2 | +4,4 ± 8,1 | +4,9 ± 10,1 | **0,90 / 0,91** | 1,18 / 1,08 | +3,3 | 1,00 ± 0,09 |
| 427929 | średnio | +13,5 ± 2,2 | | **+9,7 ± 4,5** | +10,9 ± 5,6 | | | +3,5 | |
| 428770 | 3–5 | +48,9 ± 7,8 | pol2 / pol3 | +31,1 ± 10,3 | +35,2 ± 12,9 | 0,96 / 1,04 | 0,96 / 0,95 | +21,3 | 1,14 ± 0,09 |
| 428770 | 5–8 | +12,3 ± 3,9 | pol4 / pol2 | +28,9 ± 8,9 | +7,9 ± 8,2 | 1,03 / **1,12** | 1,00 / 1,04 | +19,4 | 1,08 ± 0,10 |
| 428770 | 8–30 | +6,9 ± 1,6 | pol4 / pol4 | +5,8 ± 5,1 | +5,8 ± 5,1 | 1,03 / 1,02 | 0,94 / 0,87 | +6,2 | 0,98 ± 0,05 |
| 428770 | średnio | +22,7 ± 3,0 | | **+21,9 ± 4,9** | +16,3 ± 5,4 | | | +15,6 | |

- **(a) odtwarza K1**: w definicji K1 (1 − rel_DLZ/rel_H) dostaję dokładnie +20,7 / +9,1 / +4,2 % (427929) i +32,9 / +10,9 / +6,4 % (428770). Liczby (a) w tabeli są wyższe,
  bo definicja z planu (rel_H/rel_DLZ − 1) jest inna: g = g_K1/(1 − g_K1).
- **F-test**: H wybiera pol2 poza 428770 5–8 i 8–30 (pol4); H+DLZ wybiera pol3 w 3–5 GeV w obu runach (p = 4·10⁻⁵ i 8·10⁻¹¹), pol2 w 5–8, pol4 w 428770 8–30.
  Przejście pol2 → pol3 obniża N_S(H+DLZ) w 3–5 GeV o 5,5 % (427929: 3309 → 3126) i 8,6 % (428770: 5231 → 4780) przy niezmienionym H → N_DLZ/N_H spada z 1,08 / 1,25 do 1,02 / 1,14.
- **Zysk (b)** w 3–5 GeV spada mniej więcej o połowę względem (a) (26 → 15 %, 49 → 31 %), w 427929 5–8 i 8–30 bez zmian (ten sam rząd).
  W 428770 5–8 zysk (b) +28,9 % jest sztuczny: H dostał pol4 (p(3→4) = 0,039, p(2→4) = 0,030), H+DLZ pol2; przy równych rzędach jest +12,3 % (pol2) lub +7,9 % (pol4).
- **Pull**: szerokości 0,87–1,18 (błąd ±0,05) — błąd hesjanu pokrywa rozrzut pseudo-eksperymentów także dla pol3/pol4; σ z bootstrapu binowego danych / σ hesjanu = 0,98–1,06.
- **Obciążenie**: kryterium < 5 % niespełnione w 427929 5–8 (H 0,87) i 8–30 (H 0,90, H+DLZ 0,91) oraz 428770 5–8 (H+DLZ 1,12). Uwaga na wiarygodność — §3.
- WP50 i WP80 (kontrola, `tables.md` §2): ten sam obraz; w 3–5 GeV przy WP80 zysk (b) +31 % i +41 %, N_DLZ/N_H 1,18 i 1,31 (tam fit nadal oddaje tło sygnałowi).

## 2. Kryteria
| | 427929 | 428770 |
|---|---|---|
| 1. średni zysk (b) ≥ 10 % | **FAIL** (+9,7 ± 4,5 %) | PASS (+21,9 ± 4,9 %; bez artefaktu rzędów w 5–8: ~+16 %) |
| 2. \|obciążenie − 1\| < 5 % (H i H+DLZ, każdy bin) | **FAIL** (5–8, 8–30) | **FAIL** (5–8) |
| 3. pull width 0,8–1,2 | PASS | PASS |
| opisowo: średni zysk (b) skorygowany o obciążenie | +3,5 % | +15,6 % |

## 3. Zastrzeżenia
- **Czas**: bieg trwał 39 min zegarowych (limit ~30 min) — przekroczone; liczby replik/toyów (200/200) nie zmniejszałem, bo przekroczenie wyszło dopiero na końcu.
- **Błędy (b), (c) są zawyżone ~2×**: brak kandydatów lokalnie → bootstrap binowy niezależny dla H i H+DLZ zamiast sparowanego po zdarzeniach
  (dla pol2: ±8,7 zamiast ±4,6 pkt proc. w 427929 3–5). Stąd „+9,7 ± 4,5 %” nie odróżnia się od progu 10 %.
- **Obciążenie jest słabo wyznaczone**: prawda tła exp(pol6) z samych pasm bocznych ma swobodę pod pikiem (przerwa 90 MeV); wariant C (pasma > 70 MeV) różni się od A
  o 3–10 % (np. 428770 8–30 H: 0,97 vs 1,07; 427929 5–8 H: 0,87 vs 0,82). W skanie R2 (prawda pol4) to samo H pol2 miało obciążenie 0,90–1,01 (zależnie od wariantu). Wynik kryterium 2 zależy więc
  od konstrukcji prawdy co najmniej tak samo jak od fitu; „zysk skorygowany” ma tę samą niepewność (±5–10 pkt proc.), a w dodatku koryguje H i H+DLZ w przeciwne strony.
- Rząd wybierany osobno per klasyfikator daje niesymetryczne porównania (428770 5–8). Porównanie przy równym rzędzie (c) jest stabilniejsze: +10,9 ± 5,6 % i +16,3 ± 5,4 %.
- Sygnał w pseudo-danych parametryczny (2 Gaussy MC), bez składowych piko-podobnych — jak w skanie.
- 3–5 GeV: kształt sygnału ustalony na MC (jak w K1).

## 4. Wniosek
Z tłem dobranym do danych zysk H+DLZ nie znika, ale kurczy się do 3–5 GeV: tam +15 % (427929) i +31 % (428770) surowo, +9 % i +21 % po korekcie obciążenia; w 5–8 GeV ~+10 % surowo
i nieodróżnialny od zera po korekcie w 427929; w 8–30 GeV +4–6 %. Średnia w czystym minbiasie (+9,7 ± 4,5 %) nie spełnia progu 10 %, a test obciążenia nie przechodzi.
Do rozstrzygnięcia potrzebny bootstrap sparowany po zdarzeniach i prawdziwi kandydaci MC (slim npz i OOF z EOS `charm_q/dlz_data`).

