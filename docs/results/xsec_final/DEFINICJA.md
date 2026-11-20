# XSEC_final — definicja końcowego rachunku σ(D⁰), σ(D⁺), run 427929 (projekt B)

20.11.2026, ~00:50 CET. **Zapisane przed fitami tego zadania** (fity §3 uruchomione po zapisie tego pliku). Całość POST HOC względem XSEC, FITBIAS, BKGCAT, CVER, XSEC_atlas, SFABS — niczego nowego nie wprowadza, łączy ustalone elementy.
Tylko piki inkluzywne, bez tagu protonu. Praca lokalna, bez lxplus. Nic nie zmieniono poza tym katalogiem i notatką w bazie wiedzy.
Wyniki znane przed zapisem (nie do uniknięcia): uzyski FITBIAS/BKGCAT/XSEC_atlas, σ XSEC/CVER, porównanie z CMS 0,53 ± 0,11.

## 1. Wielkość mierzona i wzór

σ(bin) = N / (C_ekskl × SF × BR_PDG × L)

| symbol | definicja | źródło |
|---|---|---|
| σ | przekrój na produkcję D⁰ + anty-D⁰ (suma) albo D⁺ + D⁻ (suma), **prompt + nieprompt**, wszystkie rozpady, pT(D) w binie 3–5 / 5–8 / 8–30 GeV, **\|η(D)\| < 2,2** (η, nie y); D⁰ z D*⁺ wliczone | `XSEC_2026-11-19/DEFINICJA.md` §1 (bez zmian) |
| N | uzysk z fitu piku przy WP70 drzewa H, dane 427929 (`427929\|MB427929` = HLT_mb_mbts_L1MBTS_1 ∧ GRL), histogramy `WPSCAN_2026-11-19/runs/lxplus/wps_hist.npz` (kandydaci bez cięcia η reco). Szablon sygnału = fd_nd, `truth_category` 0/1, bez odbicia → N to wyłącznie prawdziwe K⁻π⁺(γ) / K⁻π⁺π⁺(γ) | §3 |
| C_ekskl | N_reco / N_gen^ekskl. **Licznik:** kandydaci mb_nd, WP70, `truth_category` ∈ {0, 1}, D⁰: `swap ≠ 1`, masa w zakresie histogramu, bin pT **reco**, bez cięcia η (3653 / 1919 / 848; 3245 / 1827 / 843; `wps_fit.json` `N_truth`; 100 % trybu ekskluzywnego wg CVER §3). **Mianownik:** prawdziwe D w mb_nd w trybie **dokładnie K⁻π⁺ / K⁻π⁺π⁺ (+ fotony FSR, m_brak < 100 MeV)**, bin pT **prawdy**, \|η\| < 2,2, bez DCS: `XSEC_atlas_2026-11-20/out/truth_ngen.json` `N_excl`; sprawdzian: CVER tab. §2 (ekskl. bez γ + z γ) — mają się zgadzać co do sztuki | CVER, XSEC_atlas |
| BR_PDG | D⁰→K⁻π⁺ **3,936 ± 0,030 %** (bez DCS: mianownik nie zawiera DCS, a DCS K⁺π⁻ jest w MC w składowej pikującej DMM szablonu M8, nie w sygnale — BKGCAT §2); D⁺→K⁻π⁺π⁺ 9,31 ± 0,14 % (PDG 2026, z FSR) | `BR_frag_2026-11-18`, CVER §2 |
| SF | dane/MC wydajności klasyfikatora przy WP70 — §4 | DNET, SFABS |
| L | 0,4648 ± 0,0103 nb⁻¹ (2,2 %; kalibracja = założenie) | `LUMI_2026-11-19/lumi.json`, XSEC2 §8.1 |

Względem `XSEC_2026-11-19`: mianownik C po prawdziwym trybie zamiast znacznika topologicznego `decayType` (σ × 0,734 D⁰, × 0,939 D⁺), BR D⁰ bez DCS (+0,4 %), uzysk D⁰ z r = 1 i szablonem z drzewa H (§3), SF z SFABS (§4), nowe wiersze budżetu (§5). XSEC nie jest zmieniany.

## 2. Rozbieżność κ — co sprawdzam (przed fitami)

BKGCAT: zmiana uzysku D⁰ 427929 przy przejściu FITBIAS → (szablon wprost, r = 1) = −1,1 / +0,4 / +2,8 %. XSEC_atlas: „przy κ × 1 uzysk rośnie o +6 % (5–8) / +3 % (8–30)”.
Hipotezy do rozróżnienia (z lektury kodu, bez nowych fitów): (a) inna baza porównania — BKGCAT porównuje z zapisem FITBIAS (szablon T interpolowany, κ × 1,364), XSEC_atlas porównuje N(κ × 1)/N(κ × 1,364) wewnątrz własnego fitu;
(b) inna próbka — tabele XSEC_atlas mają cięcie \|η_reco\| < 2,2 na kandydatach (18 187 kandydatów w 5–8), histogramy WPSCAN nie; (c) inny szablon — XSEC_atlas buduje szablon i κ z własnych tabel (kategorie 2, 3, 6, 7), BKGCAT z `bkgcat.npz`;
(d) drugi minimum fitu (inny punkt startowy; `xa_fit.dofit` startuje od zera, bez `start`).
Sprawdzian: `skrypty/kappa_fit.py` — (i) odtworzenie liczb BKGCAT (ref/dir × r_FB/1) na histogramach WPSCAN; (ii) odtworzenie liczb XSEC_atlas na jego tabelach (`in/parts`) z κ × 1 i × 1,364, z zapisem dmu, sc, −2 ln L i skanem po κ (r = 1,0 … 1,364 co 0,091) oraz fit z różnych punktów startu;
(iii) ten sam fit na tabelach XSEC_atlas bez cięcia \|η\| < 2,2. Rozstrzygnięcie = wskazanie, który z (a)–(d) tłumaczy różnicę liczbowo.

## 3. Uzysk D⁰ — nominał i warianty (kryteria zapisane z góry)

- **Nominał:** model M8 z FITBIAS bez żadnego dostrajania (`fb_fit.fit`, `fb_var.fit_set`: ograniczenie Gaussa dmu 0,8 MeV / sc 10 %, odbicie o stałym f_refl, składowa pikująca DMM + CM o stałym κ, tło exp(Czebyszew 2), zakres 1,70–2,05 GeV, biny 5 MeV),
  **szablon pikujący i κ wprost z drzewa H** (`BKGCAT_2026-11-20/out/peak_templates_H.npz`, klucze `D0|<bin>|WP70|{h5,kap}`, wczytanie przez `bc_m8.patch(mode="dir")`), **r = 1**.
- **Wariant szablonu:** szablon T interpolowany (FITBIAS, `mode="ref"`), r = 1. Systematyka „szablon pikujący” = |N_ref − N_dir| / N_dir.
- **Normalizacja κ:** refit przy κ × (1 ± 0,058) (5 % BR ⊕ 3 % fragmentacja, wartość z XSEC_atlas KRYTERIA §C) → ½|N₊ − N₋|/N. Zastępuje wiersz „norm. składowej pikującej” XSEC_atlas (liczony na innej próbce). Informacyjnie ± 20 %.
- **Błąd statystyczny:** max(rms 30 replik bootstrapu po zdarzeniach z `wps_hist.npz`, 1,3 × HESSE) — FITBIAS: rms pull M8 = 1,28.
- **Model fitu:** max |ΔN/N| wobec M4 i M9 (ten sam szablon wprost, r = 1). Wariant „κ bez skali BR” z FITBIAS odpada (jest teraz nominałem).
- **Obciążenie fitu:** domknięcie na mb_nd (M8, szablon wprost, WP70): b = N_fit/N_prawda − 1; systematyka |b| ⊕ σ_b, **nominał bez poprawki** (jak XSEC); wariant N/(1 + b) tylko informacyjnie.
- Fit uznaję za ważny, jeśli `valid` i nie na granicy; p(χ²) raportowane, bez progu odrzucenia (FITBIAS: 8–30 GeV p = 0,04 było znane).
- Biny: 3–5, 5–8, 8–30 GeV. 428770 poza zakresem zadania.
- **D⁺:** bez zmian — M0, N = 3547 / 1651 / 687, błąd max(HESSE, bootstrap), model tła i obciążenie z `XSEC_2026-11-19/out/xsec_results.json` (nie przefitowuję). Informacyjnie podaję błąd stat. × 1,3 (rms pull M0 = 1,26), nie wchodzi do nominału.

## 4. SF

| komórka | SF nominalny | źródło |
|---|---|---|
| D⁰ 5–8 | 1,005 ± 0,126 | tag D* w 427929, WP70 (`DNET/runs/lxplus4/dnet4_results.json`; jak XSEC) |
| D⁰ 8–30 | 0,965 ± 0,097 | jw. |
| D⁰ 3–5 | 0,967 ± 0,033 | SFABS `out/final.json` `D` (łączny) |
| D⁺ 3–5 / 5–8 / 8–30 | 0,949 ± 0,048 / 0,993 ± 0,052 / 0,990 ± 0,051 | SFABS `D` (łączny; drabina + przeniesienie) |

Warianty D⁺ pokazywane obok nominału: (i) SF = 1 ± 0,10 (dawny scenariusz), (ii) wartość SFABS z błędem ± 0,08 (wariant ostrożny). Niezależny sprawdzian (nie wchodzi do średniej): selekcja cięciowa XSEC_atlas 1/R = 0,97 ± 0,05 ± 0,05 (D⁺, 5–30 GeV), 0,95 ± 0,06 ± 0,08 (D⁰).

## 5. Budżet systematyczny (wiersze; ⊕ kwadratowo)

| wiersz | wartość / sposób | źródło | uwagi o podwójnym liczeniu |
|---|---|---|---|
| fit: model | D⁰: §3 (M4, M9); D⁺: max z 4 wariantów tła | §3; XSEC `var` | |
| fit: obciążenie (domknięcie) | D⁰: §3; D⁺: 7,8 / 12,4 / 5,6 % | §3; XSEC | |
| szablon pikujący (kształt) | D⁰: dir vs ref | §3 | zastępuje „κ bez skali” w fit_model FITBIAS |
| normalizacja κ ± 5,8 % | D⁰: refit | §3 | zastępuje wiersz XSEC_atlas (0,9 / 0,4 %) |
| SF | §4 | | SF z danych zawiera różnicę dane/MC dE/dx i rozdzielczości **w cięciu klasyfikatora** → **wiersz dE/dx XSEC_atlas (0,1–0,5 %) pominięty** (byłby liczony podwójnie) |
| stat. MC (C) | dwumianowy na N_gen^ekskl | truth_ngen | |
| BR | 0,76 % / 1,5 % | PDG | |
| odbicie K↔π | ≤ 0,12 % | XSEC | |
| nieprompt | 0,4–3,1 % | XSEC | |
| wybór zdarzeń + tryger | 0,2 % | XSEC2 | |
| śledzenie | 2,6 % na tor, skorelowane między torami: D⁰ 5,2 %, D⁺ 7,8 % | XSEC_atlas (1512.02913) | nowy; nie ma go w SF (SF = tylko klasyfikator przy danej rekonstrukcji) |
| rozdzielczość d₀/z₀ + materiał | 1,5 % | XSEC_atlas (2302.00336) | nowy; dotyczy preselekcji/wierzchołka; część klasyfikatorowa jest w SF — zostawiam w całości (zachowawczo) |
| rekonstrukcja wierzchołka 3 torów | 3,0 %, tylko D⁺ | XSEC_atlas | nowy |
| model generatora (pT, η ⊕ krotność) | 2,1 / 2,1 / 2,2 % (D⁰), 2,5 / 2,3 / 2,2 % (D⁺); D⁰ 3–5: z `gen_model.json` (1,2 %) ⊕ 2,0 % | XSEC_atlas | nowy |
| świetlność | 2,2 % | osobno | |

Wariant informacyjny: śledzenie Run 2 (1,35 %/tor).

**Stosunek D⁺/D⁰** = [N⁺/(C⁺ SF⁺ BR⁺)] / [N⁰/(C⁰ SF⁰ BR⁰)]. Skracają się: L, wybór zdarzeń. Częściowo: śledzenie → 1 tor (2,6 %); d₀/materiał → 1,0 %; wierzchołek 3,0 % (zostaje); model generatora → 1,5 %; nieprompt — jak XSEC T2 (kwadratowo, bez skracania);
SF: 3–5 GeV — SF⁺/SF⁰ = 0,981 ± 0,042 (wspólny SF_ref się skraca; SFABS `ratio_3-5`); ≥ 5 GeV — SF⁰ z tagu D* 427929, SF⁺ z SFABS; wspólna część: komórki 427929 w SF_ref (waga ~0,06–0,07) → kowariancja w·σ²_komórki odjęta. Pozostałe wiersze nieskorelowane (fit, stat MC, BR z korelacją 0).
**Stosunek D*/D⁰:** σ(D*) z `XSEC_atlas` (bez zmian: 134,4 ± 18,1 ± 23,1 ± 3,0 i 60,6 ± 7,1 ± 6,8 ± 1,3 µb). Skracają się: L, 2 z 3 torów (zostaje 2,6 %), d₀/materiał; SF⁰ dotyczy tylko D⁰. Korelację statystyczną (kandydaci D* ⊂ próbka D⁰) pomijam (zachowawczo — błąd stosunku zawyżony).

## 6. Porównania (definicje zapisane z góry)

Istotność odchylenia: (dane − odniesienie)/√(σ²_dane + σ²_odn), σ_dane = stat ⊕ syst ⊕ L (dla stosunków bez L), plus niepewność przeliczenia tam, gdzie jest.
1. **CMS 13 TeV (arXiv:2107.01476):** dane → definicja CMS: × f_η (ułamek D z \|η\| < 2,1 wśród \|η\| < 2,2, policzony z prawdy mb_nd per bin; wariant 2,1/2,2) × (1 − f_np) (mb_nd per bin; niepewność ± 50 % f_np) ÷ K_E (FONLL 13,6/13 TeV: 1,0285 / 1,036 / 1,0455 — jak XSEC2). CMS: Σ (dσ/dpT × Δ) po binach 5–6, 6–7, 7–8 oraz 8–12 … 24–40.
   Przed porównaniem lista kontrolna z liczbą i źródłem dla każdego czynnika: η vs y i zakres; prompt vs prompt + nieprompt; suma vs średnia cząstka/antycząstka; szerokości binów (czy tabela to dσ/dpT); 30–40 GeV; energia; świetlność w systematyce CMS; BR użyte przez CMS.
2. **FONLL 13,6 TeV, \|η\| < 2,2:** prompt = dane × (1 − f_np) wobec 2 f(c→D) σ_FONLL, **pasmo (min–max)**; odchylenie liczone do najbliższej krawędzi pasma (0, jeśli w paśmie) i osobno stosunek do wartości centralnej. f = 0,609 / 0,240.
3. **Pythia 8 A3:** wszystkie D gatunku w obszarze (dowolny rozpad): N_gen_all × σ_ND / N_zdarzeń ÷ f_ND — niezależne od `decayType` (CVER §1: spójne). Liczby z XSEC `mc`.
4. **ATLAS 7 TeV (arXiv:1512.02913), D± i D*±:** definicja jak nasza (c + b, suma ładunków), \|η\| < 2,1; skalowanie K = FONLL(13,6; 2,2)/FONLL(7; 2,1) per bin; D*: z XSEC_atlas (1,84 / 2,08); D⁺: K z plików `XSEC_atlas/fonll/prompt_Dp_*` (sam prompt; udział b zmienia K < 1 %); niepewność K: różnica K_centr − K_min/min.
5. **D⁺/D⁰:** ALICE 13 TeV \|y\| < 0,5 prompt (0,458 / 0,470 / 0,430 ± ~6 %), LHCb 13 TeV (0,405 ± 0,017), e⁺e⁻ (0,395 ± 0,015), CMS 13 TeV (z tabel), Pythia A3 (0,521 / 0,528 / 0,534).
6. **D*/D⁰:** prawda mb_nd w binie pT mezonu (0,416 / 0,437), e⁺e⁻ 0,399, CMS (z tabel).

Niczego nie dostrajam do porównań; jeśli lista kontrolna z pkt 1 ujawni błąd definicji, poprawiam przeliczenie i raportuję obie wersje.
