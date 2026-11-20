# Prawda swapu K↔π z v1.2 (`truth_swap_source: pdgid`) vs heurystyka `closest_mass`

Stan: 16.11.2026. Wejście: wyciąg `hyp_extract` z h5 v1.2 (`charm_work/HYP_d0/runs/lxplus/fd_nd.parquet`, MC fd_nd, pełna próbka,
D⁰ 3–30 GeV, m 1,75–2,00 GeV). Skrypt porównania: `mark_swaps` z `closest_mass` vs `Dmeson_hypo_swap`; WP z `fd_nd.meta.json`.

## Po co
Szablon refleksji D⁰ w ficie masy (`signal.reflection_values`, truth += 100) i kształt sygnału MC (`truth_values`) zależą od tego,
który wiersz pary hipotez uznamy za zamieniony. W v1.1 nie było prawdy hipotezy → heurystyka `closest_mass` (bliższy m0 = poprawny,
działa tylko gdy obie hipotezy są truth-matched). v1.2 ma `Dmeson_hypo_swap` z truthPdgId córek.

## Wynik (fd_nd, prawdziwe D⁰ truth 0/1: 321 378 wierszy; hypo_swap = 1: 106 390, −1: 4)

`closest_mass` vs prawda (wiersze):

| | prawda: poprawna | prawda: zamieniona |
|---|---|---|
| closest_mass: sygnał | 193 118 | **23 054** |
| closest_mass: swap | **21 866** | 83 336 |

Heurystyka myli ~20 % wierszy w każdą stronę; pomyłki prawie się znoszą w liczbach, ale nie w kształtach:

| WP | pT | r = N_refl/N_sig (cm → pdgid) | RMS sygnału w ±100 MeV [MeV] | ogon sygnału \|m−m0\|>55 MeV |
|---|---|---|---|---|
| bez cięcia | 3–5 | 0,493 → 0,501 | 28,0 → 29,9 | 7,5 % → 9,1 % |
| bez cięcia | 5–8 | 0,483 → 0,492 | 27,7 → 29,6 | 7,2 % → 8,8 % |
| bez cięcia | 8–30 | 0,467 → 0,475 | 28,2 → 30,1 | 7,7 % → 9,3 % |
| WP70 | 3–5 | 0,600 → 0,609 | 26,9 → 29,2 | 6,5 % → 8,3 % |
| WP70 | 5–8 | 0,603 → 0,613 | 27,0 → 29,4 | 6,6 % → 8,6 % |
| WP70 | 8–30 | 0,567 → 0,576 | 27,9 → 30,2 | 7,3 % → 9,3 % |

Szablon refleksji: średnia i RMS prawie bez zmian (m − m0 ≈ +7…+9 MeV; RMS 68–69 → 66–67 MeV).

## Wnioski
1. **r refleksji**: +1,5 % względnie — poniżej syst. szablonu (±5 %); normalizacja refleksji praktycznie bez zmian.
2. **Kształt sygnału MC jest zawężony przez heurystykę**: RMS +8–9 %, ogon poza ±55 MeV ×1,25–1,3. `closest_mass` przenosi
   prawdziwe D⁰ z ogonów masy do „refleksji”, a swapy blisko m0 do „sygnału”.
3. Dotyczy anomalii **EV** (5–8 GeV, WP95: ε_dane/ε_MC w ±55 MeV ≈ 0,78): część niedoszacowania ogona w MC może pochodzić
   z heurystyki, nie z danych. Skala z tabeli (+~2 pp ogona) nie tłumaczy całego efektu (~20 %). **Sprawdzone z pdgid
   (sekcja „Anomalia EV” niżej): nie tłumaczy, anomalia zostaje bez zmian.**
4. Zalecenie: `signal.truth_swap_source: pdgid` dla v1.2 (sygnał i szablon refleksji z prawdy). Analogiczna zmiana dotyczy
   fitu D⁰ w charmtag (`charmtag.samples.mark_swaps`, kod identyczny w chwili podziału) — osobna decyzja w projekcie A.

## Anomalia EV z prawdą pdgid (16.11.2026)

Bieg: `charm_work/EV_dstar` — `relabel_pdgid.py` przepisuje prawdę swapu w `ev_fd.npz` (v1.1, closest_mass) na `Dmeson_hypo_swap`
z wyciągu HYP v1.2 (dopasowano 178 637 / 178 747 wierszy; sygnał→swap 1 733, swap→sygnał 6 955), potem `ev_analysis.py` na tych
samych danych (D⁰←D*, 1,23 M par) dla obu wariantów. Wyniki: `EV_dstar/lxplus/ev_analysis_cm.md` (closest_mass),
`EV_dstar/pdgid/ev_analysis.md` (pdgid).

| bin 5–8 GeV, WP95 | ε_dane cm → pdgid | ε_MC | SF95 cm → pdgid |
|---|---|---|---|
| core (<36 MeV) | 0,81±0,08 → 0,81±0,07 | 0,95 | 0,85 → 0,85 |
| ring (36–55 MeV) | 0,72±0,31 → 0,71±0,21 | 0,93 → 0,94 | 0,77 → 0,75 |
| nom (<55 MeV) | 0,79±0,08 → 0,79±0,07 | 0,95 | 0,83 → 0,83 |

Udział sygnału Δm w pierścieniu względem rdzenia (ogony masy D⁰):

| bin | dane cm → pdgid | MC prawdziwe D⁰ cm → pdgid |
|---|---|---|
| 3–5 | 0,00±0,00 → 0,20±0,22 | 0,125 → 0,143 |
| 5–8 | 0,49±0,22 → 0,54±0,17 | 0,121 → 0,137 |
| 8–30 | 0,114±0,026 → 0,113±0,027 | 0,126 → 0,144 |

Bin 8–30: SF95 nom 1,00 w obu wariantach; 3–5: bez zmian w granicach (dużych) błędów.

**Wniosek:** prawda z pdgid poszerza ogon MC (ring/core +~15 % względnie, zgodnie z tabelą wyżej), ale anomalia 5–8 GeV
(ε_dane 0,81 vs ε_MC 0,95 w rdzeniu; ring/core w danych ~4× MC, ~2,4σ) zostaje bez zmian. Heurystyka `closest_mass` jej
nie tłumaczy — źródło leży w danych/modelowaniu detektora (np. rozdzielczość masy D⁰ w 5–8 GeV, tło Δm w pierścieniu).
Tabele cech pierścienia w obu raportach mają `nan` po stronie danych (S_est po odjęciu SB < 0) — nieużyteczne w obu wariantach.

## Konfiguracja
`signal.truth_swap_source: pdgid | closest_mass` (domyślnie brak wpisu = stare zachowanie: arbiter → pdgid, closest_mass →
closest_mass, keep_both → bez oznaczania). Niezależne od `swap_resolution` (wyboru wiersza). Bez kolumny `hypo_swap` (v1.1) →
closest_mass z ostrzeżeniem. Test: `tests/test_hyp.py::test_truth_swap_source_pdgid_independent_of_choice`.
