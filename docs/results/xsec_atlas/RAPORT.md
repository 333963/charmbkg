# XSEC_atlas — uzupełnienie pomiaru inkluzywnego 427929 o elementy wzorca ATLAS (20.11.2026, projekt B)

**Status: PRELIMINARY, POST HOC względem `XSEC_2026-11-19`.** Kryteria zapisane przed fitami danych: `KRYTERIA.md`. Pełne tabele: `tabele/wyniki.md` (generowane przez `skrypty/xa_tab.py` z `out/*.json`). Nic nie zmieniono w `XSEC_2026-11-19`, nocie ani ORCH. Bez AFP.

Dane: run 427929, GRL ∧ MBTS_1, L = 0,4648 ± 0,0103 nb⁻¹. MC: mb_nd (19 220 986 zdarzeń), szablony fd_nd. Wszystkie σ: D + anty-D, prompt + nieprompt, |η(D)| < 2,2.

## 0. Definicja wydajności i niezależny sprawdzian poprawki `decayType`

C_ekskl(bin) = N_reco(selekcja, kategoria prawdy 0/1, bin pT reco, zakres mas fitu) / N_gen^ekskl(bin pT prawdy, |η| < 2,2); σ = N_fit / (C_ekskl · BR_PDG · L); BR: 3,936 % (D⁰→K⁻π⁺, bez DCS), 9,31 % (D⁺), 2,665 % (łańcuch D*).
N_gen^ekskl z listy córek `TruthDmesons_daughterPdgId` (własna ekstrakcja `skrypty/lxplus/xa_truth.py`, 81 plików ROOT mb_nd, suma zdarzeń zgodna co do sztuki): dokładnie K⁻π⁺ / K⁻π⁺π⁺ / D⁰(→K⁻π⁺)π⁺ i brak brakującej masy (lub tylko radiacyjna < 100 MeV). Tabela: `tabele/ngen_prawda.md`.

| | 3–5 | 5–8 | 8–30 | wniosek |
|---|---|---|---|---|
| D⁰: N_ekskl/N_decayType = σ_popr/σ_XSEC | 0,7337 ± 0,0030 | 0,7344 ± 0,0047 | 0,7367 ± 0,0082 | potwierdza 0,734 |
| D⁺ | 0,9386 ± 0,0016 | 0,9387 ± 0,0026 | 0,9386 ± 0,0044 | potwierdza 0,939 |
| D*⁺ (`decayType` 1) | 0,742 ± 0,006 | 0,731 ± 0,009 | 0,724 ± 0,015 | **też topologiczny** (D⁰ z niewidocznym K⁰/η/…) |

- N_gen po `decayType` odtwarza mianowniki XSEC co do sztuki (22 306 / 8706 / 2909; 21 458 / 8440 / 2913; D* 5905 / 2484 / 862). BR ekskluzywne w MC: D⁰ 3,94–3,96 %, D⁺ 9,24–9,46 %, łańcuch D* 2,64–2,76 % — zgodne z PDG. „Łańcuch D* 3,70 %” to ten sam artefakt znacznika.
- Lista córek D* w `TruthDmesons` zawiera D⁰ **i** jego córki (kod „D⁰ K⁻ π⁺ π⁺”) — pułapka przy klasyfikacji.
- Poniżej używam wyłącznie C_ekskl i poprawionych σ(D⁰): 276 ± 23 ± 42 ± 6 µb (5–8), 97 ± 8 ± 18 ± 2 µb (8–30) (`CVER_2026-11-20`).

## A. σ(D*± → D⁰π± → K∓π±π±)

Selekcja bez klasyfikatora: tory D⁰ pT > 1 GeV, |m(Kπ) − m_D⁰| < 40 MeV, L_xy(D⁰) ≥ 0,25 mm; fit Δm `dstar_refit.fit_single` (2 Gaussy z MC, skala wolna, tło progowe u^a(1 + c₁t + c₂t²)), 140,47–156 MeV, bin 0,25 MeV. Poprawka na składową pikującą z innych rozpadów D⁰ (f_pk = 1–3 % z mb_nd).

| bin pT(D*) | N_fit | C_ekskl [%] | σ [µb] ± stat ± syst ± L | dane / ATLAS 7 TeV × K_FONLL | dane / Pythia A3 | σ(D*)/σ(D⁰_popr) | oczekiwane D*/D⁰ (prawda mb_nd) |
|---|---|---|---|---|---|---|---|
| 3,5–5 (orientacyjny) | 321 ± 72 | 10,09 ± 0,56 | 252 ± 56 ± 45 ± 6 | 0,68 ± 0,23 | 0,82 ± 0,24 | — | — |
| 5–8 | 337 ± 45 | 19,67 ± 0,93 | **134 ± 18 ± 23 ± 3** | 0,76 ± 0,19 | 0,67 ± 0,15 | 0,49 ± 0,08 (stat) | 0,416 |
| 8–30 | 233 ± 27 | 30,77 ± 1,85 | **60,6 ± 7,1 ± 6,8 ± 1,3** | 0,90 ± 0,18 | 0,87 ± 0,14 | 0,62 ± 0,09 (stat) | 0,437 |

- Jakość: χ²/ndf 59/55, 37/55, 33/55; domknięcie na mb_nd (cały MC jako pseudodane) 0,88 ± 0,12 / 1,09 ± 0,11 / 1,06 ± 0,10 — bez obciążenia w granicach ~10 % (dokładniej tą statystyką MC nie da się sprawdzić).
- Systematyka (tabela w `tabele/wyniki.md`): 5–8: 17 % (dominuje stabilność selekcji: wariant L_xy ≥ 0,10 daje 173 ± 26 µb, 2,1σ nieskorelowanego błędu → ½ różnicy = 14 %; bez tego wiersza 10 %); 8–30: 11 % (śledzenie 3 tory 7,7 %, stat MC 6 %, fit 4 %).
- Bin 3,5–5: N/σ_N = 4,5, ale skala szerokości w wariantach dochodzi do granicy 0,80, a warianty selekcji dają 181–346 µb → tylko orientacyjnie. Bin 3–5: 437 ± 111 µb, nie do użytku.
- ATLAS 7 TeV (arXiv:1512.02913, tab. 4; liczby odczytane z PDF pracy, `lit/1512.02913.pdf`): 217,5 / 96,3 / 32,2 µb w 3,5–5 / 5–8 / 8–30 GeV, |η| < 2,1. Skalowanie energii z FONLL (formularz sieciowy v1.3.4, NNPDF30 NLO; c→D* z f = 0,2429 + B→D*; `tabele/fonll_skalowanie.md`): K = 1,70 / 1,84 / 2,08 (wybór f bez znaczenia; min/min skali: 1,54 / 1,73 / 2,08). Suma 3,5–30 GeV: dane 447 ± 60 (stat) wobec 614 ± 60 → **0,73 ± 0,13**.
- FONLL bezwzględnie (×2 dla c + c̄; kontrola: 7 TeV 3,5–20 GeV daje ≈ 186 µb wobec 202 µb cytowanych w pracy ATLAS): 111 (69–174) µb w 5–8, 43 (32–58) µb w 8–30 → dane/FONLL 1,2 / 1,4, w paśmie.
- D*/D⁰: oczekiwanie z frakcji fragmentacji f(c→D*⁺)/f(c→D⁰) = 0,399 (e⁺e⁻, 1509.01061) / 0,408 (ALICE 13 TeV, 2308.04877); wprost z prawdy mb_nd w tych samych binach pT mezonu 0,416 / 0,437. Dane: 0,49 ± 0,08 (+0,9σ stat) i 0,62 ± 0,09 (+2,1σ stat; ≈ 1,4σ z systematyką nieskorelowaną). Równoważnie D* „przewiduje” σ(D⁰) = 323 ± 44 µb (5–8) i 139 ± 16 µb (8–30) wobec 276 ± 23 i 97 ± 8 (stat).

## B. σ(D⁰), σ(D⁺) z selekcją cięciową bez klasyfikatora

Cięcia z prac ATLAS (D⁺: 2302.00336 tab. 3 bez ΔR i izolacji; D⁰: tory > 1 GeV wg 1512.02913, wierzchołek wzorem D⁺, L_xy > 0,45 mm skalowane czasem życia) — pełna lista w `KRYTERIA.md` §B. Własna ekstrakcja (DAG `charm_xsatlas`, 17 węzłów). Ten sam kod fitu i te same tabele dla obu selekcji: D⁺ M0, D⁰ M8 (κ × 1; szablony po tych samych cięciach).

Kontrola: selekcja WP70 z moich tabel odtwarza XSEC — D⁺ N = 3542 / 1649 / 687 (XSEC 3547 / 1651 / 687), C(decayType) = 15,12 / 21,65 / 28,94 % (identycznie); D⁰ N = 4050 / 1507 / 663 przy κ × 1 (XSEC z κ × 1,364: 4105 / 1530 / 681; u mnie z × 1,364: 4066 / 1421 / 646 — różnica definicji κ, ≤ 7 %).

| | bin | σ_cięcia [µb] | σ_WP70 [µb] (te same tabele) | R = σ_c/σ_k ± stat(boot) ± stat MC | SF_klas = 1/R |
|---|---|---|---|---|---|
| D⁰ | 3–5 | 954 ± 127 | 1036 ± 95 | 0,92 ± 0,13 ± 0,02 | 1,09 ± 0,16 |
| D⁰ | 5–8 | **286 ± 28** | 288 ± 23 | 0,99 ± 0,09 ± 0,02 | 1,01 ± 0,10 |
| D⁰ | 8–30 | **107 ± 11** | 96 ± 8 | 1,11 ± 0,09 ± 0,02 | 0,90 ± 0,08 |
| D⁺ | 3–5 | 393 ± 45 | 508 ± 37 | 0,77 ± 0,11 ± 0,03 | 1,29 ± 0,19 |
| D⁺ | 5–8 | **158 ± 11** | 165 ± 8 | 0,96 ± 0,07 ± 0,03 | 1,04 ± 0,08 |
| D⁺ | 8–30 | **57,8 ± 3,9** | 51,5 ± 3,0 | 1,12 ± 0,08 ± 0,03 | 0,89 ± 0,07 |
| D⁰ | 5–30 | | | 1,05 ± 0,07 | **0,95 ± 0,06** |
| D⁺ | 5–30 | | | 1,04 ± 0,06 (χ² = 2,1/1) | **0,97 ± 0,05** |

- Błąd R: wspólny bootstrap po zdarzeniach (200 replik, te same wagi Poissona per zdarzenie dla obu selekcji; korelacja uzysków ρ = 0,2–0,6) ⊕ stat MC z uwzględnieniem części wspólnej. Błędy σ w tabeli: tylko stat (HESSE).
- Warianty zapisane z góry (L_xy luźniej/ciaśniej, bez σ_3D, pT torów): D⁺ R stabilne (5–8: 0,98–1,01; 8–30: 1,00–1,10; 3–5: 0,78–0,79); D⁰ rozrzut większy (5–8: 0,91–1,07; 8–30: 0,92–1,16), żaden wariant > 2σ. Przyjmuję niepewność modelowania cięć ≈ 5 % (D⁺) i ≈ 8 % (D⁰) — z rozrzutu, nie z pomiaru.
- Domknięcia fitów na mb_nd: wszystkie w ≤ 1,4σ od 1 poza D⁺ WP70 5–8 (1,11 ± 0,05) i D⁰ WP70 3–5 (0,78 ± 0,10; znane z XSEC). Poprawka domknięć przesuwa R do 0,87 (D⁰ 5–8) i 0,99 (D⁺ 5–8) — w granicach błędu.
- **Co to mówi o SF:** dla D⁺ (dotąd scenariusz SF = 1 ± 10 %) jest to pierwsze niezależne ograniczenie: SF⁺(WP70) = 1,04 ± 0,08 (5–8), 0,89 ± 0,07 (8–30), razem 0,97 ± 0,05 (stat) ± 0,05 (modelowanie cięć) — scenariusz potwierdzony, niepewność 10 % adekwatna. Dla D⁰: 0,95 ± 0,06 ± 0,08, zgodne z SF z tagu D* w XSEC (1,005 ± 12,5 %, 0,965 ± 10,1 %). Bin 3–5 D⁺: R = 0,77 ± 0,12 (1,9σ) — fit WP70 ma tam p = 0,04 i skalę szerokości 1,16; traktować jako ostrzeżenie o uzysku WP70 3–5, nie jako pomiar SF.
- **Poziom σ(D⁰):** selekcja bez klasyfikatora daje 286 ± 28 i 107 ± 11 µb — **potwierdza poziom 276 / 97 µb**. Klasyfikator (ani jego SF) nie jest przyczyną niskiego wyniku względem CMS.

## Wniosek o napięciu z CMS (D⁰ 5–8 GeV: dane/CMS = 0,53 ± 0,11)

1. Cięcia prostokątne potwierdzają 276 µb (R = 0,99 ± 0,10) → przyczyna nie leży w klasyfikatorze.
2. D* jest nieco wyżej: równoważnik σ(D⁰) 323 ± 44 µb (5–8, +1σ) i 139 ± 16 µb (8–30, +2σ), ale sam D* leży na 0,76 ± 0,19 i 0,90 ± 0,18 przeskalowanego ATLAS 7 TeV oraz 0,67 / 0,87 Pythii. D* **łagodzi, nie usuwa** napięcia: wszystkie trzy niezależne od klasyfikatora pomiary (D⁰ cięcia, D⁺ cięcia, D*) leżą na 0,6–0,9 odniesień.
3. Wspólne dla wszystkich: świetlność (kalibracja 2,2 % = założenie), wydajność śledzenia i rekonstrukcji wierzchołka dane/MC (nie zmierzona u nas, w budżecie z literatury), preselekcja ntupli, model fitu D⁰ (M8; w 8–30 p = 0,01 dla cięć). To są kandydaci, nie rozstrzygnięcie.

## C. Budżet (własna kopia; `tabele/wyniki.md` §C)

| | D⁰ 5–8 | D⁰ 8–30 | D⁺ 3–5 (scen.) | D⁺ 5–8 (scen.) | D⁺ 8–30 (scen.) |
|---|---|---|---|---|---|
| razem syst XSEC [%] | 15,3 | 18,1 | 13,5 | 16,4 | 14,2 |
| śledzenie (2,6 %/tor, skorelowane) | 5,2 | 5,2 | 7,8 | 7,8 | 7,8 |
| rozdzielczość d₀/z₀ + materiał (lit.) | 1,5 | 1,5 | 1,5 | 1,5 | 1,5 |
| rekonstrukcja wierzchołka (lit.) | — | — | 3,0 | 3,0 | 3,0 |
| dE/dx (∓0,023 na MC, próg WP70 stały) | 0,5 | 0,3 | 0,5 | 0,4 | 0,1 |
| normalizacja składowej pikującej ±5,8 % | 0,9 | 0,4 | — | — | — |
| model generatora (pT, η ⊕ krotność lit.) | 2,1 | 2,2 | 2,5 | 2,3 | 2,2 |
| **razem syst (nowy) [%]** | **16,4** | **19,0** | **16,2** | **18,6** | **16,7** |
| wariant: śledzenie Run 2 (1,35 %/tor) | 15,8 | 18,5 | 14,7 | 17,4 | 15,3 |
| wariant: pikująca ±20 % | 16,7 | 19,1 | — | — | — |

- σ(D⁰) z nowym budżetem: **276 ± 23 (stat) ± 45 (syst) ± 6 (L) µb** (5–8), **97 ± 8 ± 18 ± 2 µb** (8–30). D⁺ (scen. SF = 1 ± 10 %, poprawione): 509 ± 40 ± 82, 165 ± 8 ± 31, 51,5 ± 4,5 ± 8,6 µb.
- D⁺/D⁰ (bez SF): 18,3 → 18,8 % (5–8), 21,6 → 22,0 % (8–30); nowe wiersze nieskorelowane ≈ 4,5 % (jeden tor, wierzchołek 3 torów).
- dE/dx: wpływ ≤ 0,6 % — drzewo H słabo zależy od dE/dx przy pT ≥ 3 GeV. Wynik odtwarza zapisane score do 8·10⁻⁴ (D⁰) / 4·10⁻³ (D⁺).
- Składowa pikująca: N(κ × 1,364)/N(κ × 1) = 0,943 (5–8), 0,975 (8–30) → dN/N = −0,16 % / −0,07 % na 1 % κ. Uwaga osobna od budżetu: uzyski XSEC mają κ × 1,364 (artefakt); przy κ × 1 uzysk D⁰ rośnie o +6 % (5–8) / +3 % (8–30) w moim ficie.
- Model generatora: czułość C na widmo pT (±10 % na krańcach binu) i η (±10 % dla |η| > 1,1) ≤ 1,2 % dla WP70, do 3,8 % dla cięć D⁺ 3–5; krotność nie do policzenia z tabel (literatura 2,0–2,5 %).

## Zastrzeżenia

- Wszystko post hoc; cięcia nie były optymalizowane na danych, ale wybór „bez ΔR” i progu L_xy D⁰ jest mój (na MC).
- Wiersze śledzenia, IP, wierzchołka i krotności są z literatury ATLAS (7 TeV / Run 2), nie z pomiaru na tych danych — to największa niezmierzona część budżetu i główny kandydat na wspólne przesunięcie względem CMS.
- SF = 1/R zakłada, że MC opisuje wydajność cięć prostokątnych; obie selekcje dzielą śledzenie, wierzchołkowanie i preselekcję, więc R nie widzi błędów wspólnych.
- D* 5–8: niestabilność względem L_xy (2,1σ) wpisana jako 14 %; fit w 3,5–5 słaby. Domknięcia Δm sprawdzone tylko do ~10 %.
- FONLL: ×2 dla c + c̄ i f(c→D*) = 0,2429 przyjęte przeze mnie; K nie zależy od f, ale pasmo skal traktowane jako skorelowane między energiami.
- Porównanie z CMS nie było tu przeliczane (liczba 0,53 ± 0,11 z `CVER_2026-11-20`).
- Wyników `SFABS_wynik_2026-11-20` nie czytałem — pkt B jest niezależny.

## Pliki

`KRYTERIA.md`; `skrypty/` (`xa_lib.py`, `xa_truth_tab.py`, `xa_dstar.py`, `xa_fit.py`, `xa_dedx.py`, `xa_gen.py`, `xa_fonll.py`, `xa_tab.py`, `lxplus/{xa_extract.py,xa_truth.py,xa_job.sh,xa.sub,make_dag_xa.sh}`); `out/*.json`; `tabele/{wyniki,ngen_prawda,dedx,model_generatora,fonll_skalowanie}.md`; `in/{parts,truth}` (kopie z EOS); `fonll/` (surowe odpowiedzi formularza + `run_atlas.sh`); `lit/1512.02913.pdf`; `NOTA_tekst_EN.md`.
