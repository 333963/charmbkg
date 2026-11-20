# Wydajność weta wierzchołków `nPileupVertices == 0` — próbka inkluzywna i MC (18.11.2026)

Punkt otwarty głównego łańcucha analizy („ε weta wierzchołków — nie zaczęte”, `GITLAB_2026-11-17/RAPORT.md` §2(b), §2(f)).
Weto w kodzie charmpp: `nPileupVertices == 0`, stosowane w danych i w MC (`src/HistogramManager.cxx:577`), razem z tagiem AFP i ξ_cal.
Tutaj: **tylko próbka inkluzywna**, bez żadnej zmiennej protonowej (AFP, ξ). Nic nie uruchamiano na lxplus.

## 0. Wynik w skrócie

| run | próbka | ε_veto zdarzeń z D⁰ (tło) | ε_veto sygnału D* (pasma Δm odjęte) | MC: sygnał D* (prawda, fd_nd) | dane/MC |
|---|---|---|---|---|---|
| 428770 | cała | 94,08 ± 0,02 % | 93,95 ± 0,73 % | 99,25 ± 0,04 % | 0,949 |
| 428770 | J12only ∪ MBTSbased | 94,72 ± 0,02 % | 94,68 ± 0,75 % | 99,25 ± 0,04 % | 0,955 |
| 427929 | GRL (MBTS) | 98,48 ± 0,02 % | 98,60 ± 0,69 % | 99,25 ± 0,04 % | 0,993 |

- Strata ma dwa człony: **f** ≈ 1,0–1,3 % (wierzchołek „pile-up” bez pile-upu; jest też w MC: 0,8 %) oraz **człon pile-up** 1 − ⟨exp(−ε_v μ)⟩ =
  4,0 % w 428770 (⟨μ⟩ = 0,054) i 0,56 % w 427929 (⟨μ⟩ = 0,0085). MC ma μ = 0, więc członu pile-up nie zawiera.
- Model P(npu = 0 | μ) = (1 − f)·exp(−ε_v μ) opisuje dane (χ²/ndf 8–12/14); ε_v = 0,77 ± 0,01 (428770, J12 ∪ MBTS) i 0,66 ± 0,03 (427929).
- Sygnał D* i zdarzenia tła z kandydatem D⁰ mają to samo ε_veto: różnica w danych +0,01 ± 0,50 % (oba runy), w MC +0,06 ± 0,04 % przy tej samej krotności.
- Poprawka do normalizacji względem MC: **×0,955 ± 0,005 (428770, bez klasy „other”), ×0,993 ± 0,003 (427929)**. Pochodzenie niepewności w §4.
- W stosunku tagowane/inkluzywne weto **nie znosi się** bez dodatkowych założeń (§4); część zależy od regionu tagowanego i nie została policzona.

## 1. Dane wejściowe — dokładna zmienna, lokalnie

Lokalnie jest dokładnie ta zmienna: `EventInfo_nPileupVertices_NOSYS` z ROOT v2.0, przepisana przez konwersję h5 v2.0 i etap ORCH `btrig` (17.11) do
`B_trig/runs/lxplus/bt_{data_428770,data_427929,fd_nd,mb_nd}.npz` (kolumny `e_npu` na wierszu kandydata, `ev_npu` per zdarzenie) oraz
`trg_data_*.npz` (`npu` dla **wszystkich** zdarzeń drzewa, także bez kandydata). Obok: `actual_mu`, `lumiBlock`, `passGRL_MinBiasAFP`, `nSelBL`, bity 11 łańcuchów trygera.
Zastępniki nie były potrzebne. `tj2_*.npz` niosą tę samą kolumnę; `AFP_mix`, `DLZ_dplus/SCAN_2026-11-17/cand` nie były używane (zawierają zmienne AFP).

| próbka | zdarzenia w drzewie (GRL) | z kandydatem D⁰ | pary D⁰←D* po cięciach | uwagi |
|---|---|---|---|---|
| data_428770 | 7 988 801 | 1 297 217 | 56 913 | 97 plików; ⟨μ⟩ = 0,054 (1–99 %: 0,017–0,165) |
| data_427929 | 23 851 045 | 608 181 | 15 654 | 458 plików; tylko GRL (LB 119–324); ⟨μ⟩ = 0,0085 (0,002–0,028) |
| fd_nd | – | 700 176 | 65 044 (32 967 sygnału z prawdy) | μ = 0, LB = 1 |
| mb_nd | – | 506 644 | 12 376 (770 sygnału z prawdy) | μ = 0, LB = 1 |

Definicje: „zdarzenie z D⁰” = ≥ 1 kandydat D⁰→Kπ w oknie 1,75–2,00 GeV, pT ≥ 3 GeV, L_xy ≥ 0,25 mm (tabela `all_`, bez cięcia klasyfikatora) — to próbka tła.
„Sygnał D*” = D⁰←D* z |m(Kπ) − m_D⁰| < 40 MeV, H > 0,6, pT 3–30 GeV, okno Δm 144,4–146,6 MeV, tło odjęte z pasm Δm wagami `iso_common.sbs_w`
(fit progowy pasm w grupach pass/fail × pT). W MC sygnał z prawdy. Błędy: sygnał — bootstrap Poissona po zdarzeniach (300 replik, wagi SBS stałe);
ułamki zdarzeń — dwumianowe (sprawdzone: bootstrap po zdarzeniach daje to samo, 0,0516 vs 0,0519 %).
Klasy trygera jak w `btrig_common.classify`: 428770 — J12only, MBTSbased (podklasy j20MBTS, sptrk4), other; 427929 — MB427929 (= cała próbka GRL), sptrk4.

## 2. ε_veto w danych

### 2.1 Per run i klasa (pełne tabele: `out/tables.md`)

| run | klasa | ⟨μ⟩ | ⟨nSelBL⟩ | N_S | sygnał D* [%] | pasma Δm [%] | zdarzenia z D⁰ [%] | N zdarzeń |
|---|---|---|---|---|---|---|---|---|
| 428770 | all | 0,0541 | 43,6 | 2797 | 93,95 ± 0,73 | 92,89 ± 0,13 | 94,08 ± 0,02 | 1 297 217 |
| 428770 | J12only | 0,0539 | 44,7 | 2193 | 94,72 ± 0,82 | 93,42 ± 0,14 | 94,78 ± 0,02 | 951 433 |
| 428770 | MBTSbased | 0,0541 | 42,4 | 372 | 93,92 ± 2,13 | 93,19 ± 0,34 | 94,45 ± 0,05 | 216 517 |
| 428770 | j20MBTS | 0,0541 | 42,4 | 340 | 94,27 ± 2,08 | 93,26 ± 0,35 | 94,53 ± 0,05 | 198 739 |
| 428770 | sptrk4 | 0,0541 | 43,1 | 312 | 94,00 ± 2,20 | 93,24 ± 0,38 | 94,45 ± 0,06 | 165 562 |
| 428770 | other (tylko tło) | 0,0553 | 37,8 | – | – | – | 88,24 ± 0,09 | 129 267 |
| 427929 | all = MB427929 | 0,0085 | 33,3 | 954 | 98,60 ± 0,69 | 97,85 ± 0,14 | 98,48 ± 0,02 | 608 181 |
| 427929 | sptrk4 | 0,0085 | 37,5 | 508 | 98,20 ± 0,91 | 97,51 ± 0,23 | 98,27 ± 0,03 | 255 175 |

Klasa „other” (ani J12, ani żaden łańcuch MBTS; 10,0 % zdarzeń z D⁰; 24,6 % z nich ma jeden z dwóch zapisanych łańcuchów L1AFP) jest podana tylko jako ułamek
zdarzeń tła — bez uzysku D* i bez zmiennych protonowych. Jej ε_veto jest o 6,5 punktu niższe niż reszty próbki (omówienie w §2.3 i §4).

### 2.2 Zależność od μ, LB, krotności, pT, L_xy

- **μ:** 428770 — 96,90 % (μ = 0,021) → 94,40 % (0,050) → 91,12 % (0,090) → 81,3 % (0,21); 427929 — 98,86 % (0,0025) → 97,51 % (0,024). Zależność wykładnicza (§2.3).
  Sygnał D* w binach μ zgodny z modelem dopasowanym do zdarzeń z D⁰ bez swobodnych parametrów: χ²/n = 5,2/8 (428770), 2,0/7 (427929).
- **LB (okresy runu):** 428770 — 94,22 / 94,19 / 94,04 / 93,90 % w czterech okresach (⟨μ⟩ 0,0548 / 0,0523 / 0,0546 / 0,0551); 427929 — 98,51 → 98,45 % (⟨μ⟩ 0,0077 → 0,0091).
  Spadek o 0,3 punktu w 428770 nie wynika z samej zmiany ⟨μ⟩ (ta dałaby ≤ 0,2 punktu); przyczyny nie ustaliłem.
- **nSelBL:** silna — 428770: 98,0 % (10–20) → 94,6 % (40–50) → 82,4 % (> 80); 427929: 99,5 % → 97,9 % → 93,4 %. Dwa powody: f rośnie z krotnością (także w MC),
  a nSelBL liczy tory całego zdarzenia, więc zdarzenia z pile-upem mają większe nSelBL. Bin nSelBL nie jest zatem niezależny od weta
  (pozorne „ε_v” w binach nSelBL rośnie od 0,2 do 2,8). Do przenoszenia poprawki między próbkami nadaje się tylko f(nSelBL) przy μ → 0.
- **pT(D⁰):** zdarzenia z D⁰ — 428770: 93,96 / 94,23 / 94,64 % (3–5 / 5–8 / 8–30 GeV); 427929: 98,47 / 98,24 / 98,10 %. Sygnał D*: 91,2 ± 2,1 / 95,6 ± 1,2 / 95,1 ± 0,7 % (428770),
  98,8 ± 1,1 / 98,6 ± 0,8 / 97,6 ± 1,4 % (427929). Wzrost z pT w 428770 odpowiada zmianie składu trygerów (udział J12 rośnie z pT), nie kandydatowi.
- **L_xy(D⁰):** sygnał w danych — 428770: 92,2 ± 2,8 / 90,0 ± 2,0 / 95,5 ± 1,3 / 96,4 ± 1,0 / 95,0 ± 1,2 % (χ²/ndf wobec stałej 9,2/4, p ≈ 0,06); 427929: 2,9/4. Pasma Δm płaskie.

### 2.3 Oczekiwanie Poissona

Dla rzadkiego procesu liczba dodatkowych zderzeń w tym samym przecięciu ma rozkład Poissona(μ), więc P(npu = 0 | μ) = (1 − f)·exp(−ε_v μ),
gdzie ε_v to efektywna wydajność rekonstrukcji wierzchołka dodatkowego zderzenia, a f — ułamek zdarzeń z wierzchołkiem „pile-up” bez pile-upu.
Dopasowanie na zdarzeniach z D⁰ (biny kwantylowe μ, wiarygodność dwumianowa):

| próbka | f [%] | ε_v | χ²/ndf | strata pile-up 1 − ⟨exp(−ε_v μ)⟩ [%] |
|---|---|---|---|---|
| 428770 J12only | 1,24 ± 0,06 | 0,765 ± 0,012 | 8,1/14 | 4,03 ± 0,06 |
| 428770 MBTSbased | 1,44 ± 0,13 | 0,789 ± 0,025 | 20,8/14 | 4,16 ± 0,13 |
| 428770 J12only ∪ MBTSbased | 1,27 ± 0,05 | 0,770 ± 0,011 | 10,4/14 | 4,05 |
| 428770 all | 1,25 ± 0,05 | 0,899 ± 0,011 | 8,6/14 | 4,72 ± 0,05 |
| 428770 other | 1,62 ± 0,24 | 1,98 ± 0,05 | 16,2/14 | 10,3 ± 0,2 |
| 427929 GRL | 0,97 ± 0,02 | 0,656 ± 0,025 | 11,5/14 | 0,56 ± 0,02 |

- Dane spełniają model w każdej klasie. ε_v ≈ 0,66–0,79 to rozsądna wartość dla ułamka zderzeń nieelastycznych dających zrekonstruowany wierzchołek.
- ε_v w 428770 i 427929 różni się o 0,109 ± 0,028 (3,9σ). Możliwe przyczyny: skala `actual_mu` przy μ ~ 0,005 albo inne warunki wiązek; nie rozstrzygałem.
  Jednej uniwersalnej ε_v nie należy więc przenosić między runami (w 428770 zmieniłoby to stratę o 0,7 punktu) — ε_veto trzeba brać z danego runu.
- Klasa „other” ma nachylenie 1,98 > 1, czego czysty Poisson nie dopuszcza: wybór tych zdarzeń jest skorelowany z obecnością dodatkowego zderzenia.
  Przez nią „all” ma ε_v = 0,90 mimo 0,77 w obu głównych klasach.
- Poza GRL w 427929 (μ ≈ 22) ⟨npu⟩/⟨μ⟩ = 0,556 — niżej niż 0,66 przez łączenie wierzchołków przy dużym μ; tylko informacyjnie.

## 3. MC i porównanie dane/MC

**Profil μ w MC:** `actual_mu` = 0 i `lumiBlock` = 1 we wszystkich zdarzeniach fd_nd i mb_nd — MC nie ma pile-upu. Mimo to 0,7–1,0 % zdarzeń ma npu ≥ 1
(fd_nd: zdarzenia z D⁰ 99,02 ± 0,01 %, sygnał D* 99,25 ± 0,04 %; mb_nd: 99,19 ± 0,01 % i 99,48 ± 0,24 %). To jest człon f.

**Czy porównanie ma sens:** bezpośredni stosunek ε_veto dane/MC (0,949 / 0,955 / 0,993) nie testuje modelowania — mierzy głównie pile-up, którego w MC nie ma.
Sens ma porównanie członu f przy μ → 0 z MC, w binach krotności:

| nSelBL | mb_nd [%] | fd_nd [%] | fd_nd sygnał D* [%] | 428770 fit [%] | 428770/mb_nd | 427929 fit [%] | 427929/mb_nd |
|---|---|---|---|---|---|---|---|
| 10–20 | 0,25 ± 0,02 | 0,30 ± 0,02 | 0,29 ± 0,06 | 0,30 ± 0,12 | 1,21 ± 0,50 | 0,25 ± 0,03 | 1,01 ± 0,14 |
| 20–30 | 0,48 ± 0,02 | 0,56 ± 0,02 | 0,48 ± 0,08 | 0,51 ± 0,09 | 1,07 ± 0,20 | 0,50 ± 0,04 | 1,04 ± 0,09 |
| 30–40 | 0,83 ± 0,03 | 0,88 ± 0,02 | 0,85 ± 0,12 | 0,71 ± 0,10 | 0,86 ± 0,12 | 0,99 ± 0,05 | 1,19 ± 0,07 |
| 40–50 | 1,19 ± 0,04 | 1,25 ± 0,03 | 1,30 ± 0,17 | 1,17 ± 0,12 | 0,98 ± 0,10 | 1,39 ± 0,07 | 1,18 ± 0,07 |
| 50–65 | 1,59 ± 0,05 | 1,70 ± 0,04 | 1,41 ± 0,21 | 1,99 ± 0,14 | 1,25 ± 0,10 | 1,76 ± 0,09 | 1,11 ± 0,07 |
| 65–80 | 2,35 ± 0,14 | 2,17 ± 0,09 | 2,11 ± 0,53 | 2,61 ± 0,25 | 1,11 ± 0,13 | 2,86 ± 0,21 | 1,22 ± 0,12 |
| > 80 | 3,21 ± 0,42 | 2,80 ± 0,24 | 1,44 ± 1,01 | 4,06 ± 0,47 | 1,27 ± 0,22 | 3,59 ± 0,46 | 1,12 ± 0,20 |

MC opisuje f(nSelBL) w danych z dokładnością 0–25 % (dane wyżej średnio o ~10–15 %). Łącznie f: dane 1,0–1,3 %, MC 0,8 %. W 428770 część różnicy to większa krotność (⟨nSelBL⟩ 43,6 wobec 32,9 w mb_nd); w 427929 krotność jest taka sama (33,3), więc tam różnica to modelowanie.
Bin nSelBL < 10 w 427929 odstaje (0,90 ± 0,09 % wobec 0,10 % w MC, ujemne nachylenie w μ) — patrz niżej „tło wiązki”.

**Zależność od samego kandydata (punkt 3 zadania):**
- MC, sygnał D* z prawdy wobec zdarzeń z D⁰ o tej samej krotności: +0,06 ± 0,04 % (fd_nd), +0,22 ± 0,24 % (mb_nd); pasma Δm: −0,02 ± 0,09 %, −0,19 ± 0,12 %.
  Tory rozpadu D nie tworzą dodatkowych wierzchołków „pile-up” na poziomie > 0,1 %.
- MC, sygnał wobec L_xy: 99,43 ± 0,09 % (0,25–0,4 mm) → 99,03 ± 0,13 % (1–2 mm) → 98,73 ± 0,28 % (> 2 mm); χ²/ndf wobec stałej 10,4/4 (p ≈ 0,03).
  Słaby spadek z długością lotu, 0,7 punktu na całym zakresie; wobec pT płasko (99,31 / 99,15 / 99,29 %).
- Dane, sygnał D* − zdarzenia z D⁰: −0,13 ± 0,73 % (428770), +0,12 ± 0,69 % (427929), łącznie +0,01 ± 0,50 %.
- Dane, pasma Δm są niżej niż zdarzenia z D⁰ o 1,2 punktu (428770) i 0,6 punktu (427929). W 427929 znika to po wyrównaniu krotności (+0,08 ± 0,14 %),
  czyli to efekt krotności (kombinatoryka rośnie z liczbą torów). W 428770 wyrównanie nSelBL przestrzeliwuje (+1,41 ± 0,13 %), bo nSelBL zawiera tory pile-upu.
- Zdarzenia **bez żadnego kandydata** (drzewo): 428770 — 95,46 ± 0,01 % (⟨nSelBL⟩ = 13,7), fit f = 0,51 ± 0,03 %, ε_v = 0,765 ± 0,007 (J12only: 0,675 ± 0,013);
  427929 — 99,13 %, przy nSelBL ≥ 10: f = 0,34 %, ε_v = 0,633 ± 0,006. Nachylenie zgodne ze zdarzeniami z D⁰ (0,77 i 0,66), f niższe przez mniejszą krotność.
  Różnice między zdarzeniami z kandydatem i bez tłumaczy więc krotność, nie kandydat.

**Tło wiązki (uboczne):** w 427929 0,149 % zdarzeń drzewa w GRL ma npu ≥ 3 przy medianie nSelBL = 1; ich udział maleje jak 1/μ (0,37 % przy μ = 0,002, 0,025 % przy μ = 0,026),
czyli stała częstość na przecięcie — to nie pile-up, tylko zdarzenia niekolizyjne. Wśród zdarzeń z D⁰ jest ich 0,017 % (428770: 0,010 %), bez znaczenia dla pomiaru,
ale psują fit Poissona na zdarzeniach bez kandydata i w binie nSelBL < 10 (dlatego w tabelach wariant nSelBL ≥ 10).

## 4. Wniosek dla pomiaru

**Poprawka.** MC stosuje weto, ale zawiera tylko człon f. W danych dochodzi człon pile-up, niezależny od procesu twardego. Jeśli świetlność liczy wszystkie przecięcia,
uzysk trzeba podzielić przez ε_veto z danych, czyli pomnożyć wydajność MC przez SF = ε_veto(dane)/ε_veto(MC):

| run | próbka | SF (zdarzenia z D⁰, dane/mb_nd) | SF bezpośrednio z sygnału D* (dane SBS / fd_nd prawda) |
|---|---|---|---|
| 428770 | J12only | 0,9556 | 0,954 ± 0,008 |
| 428770 | MBTSbased | 0,9523 | 0,946 ± 0,021 |
| 428770 | J12only ∪ MBTSbased | 0,955 | 0,954 ± 0,008 |
| 428770 | all (z klasą „other”) | 0,9485 | 0,947 ± 0,007 |
| 427929 | GRL | 0,9928 | 0,993 ± 0,007 |

Proponowana wartość: SF ze zdarzeń z D⁰ (błąd statystyczny < 0,05 %), z niepewnością systematyczną:
- przeniesienie tło → sygnał: MC ≤ 0,1–0,2 %, dane potwierdzają tylko do ±0,5 % → **0,25 %**;
- modelowanie f w MC (dane/MC 1,0–1,25 przy f ≈ 0,8 %) → **0,2 %**;
- zależność od pT/L_xy, jeśli SF jest jedna na wszystkie biny (zakres 93,96–94,64 % w 428770; trend L_xy 0,7 punktu w MC) → **0,3 %** (428770), 0,2 % (427929);
- skład klas trygera w 428770 (J12only 94,78 %, MBTSbased 94,45 %) → **0,2 %**, pod warunkiem że regiony trygera są zdefiniowane; bez nich dochodzi
  0,65 punktu różnicy między „all” a „J12 ∪ MBTS”.
Razem: **428770: SF = 0,955 ± 0,005; 427929: SF = 0,993 ± 0,003.** To oszacowanie budżetu, nie wynik fitu; składniki są moimi przyjętymi wartościami na podstawie liczb powyżej.
Alternatywa bez założeń o przenoszeniu: pomiar wprost na sygnale D* z ±0,7 %.

**Czy znosi się w stosunku tagowane/inkluzywne.** Nie automatycznie. Znosi się człon pile-up prawdziwie tagowanych zdarzeń (proton i D z tego samego zderzenia),
jeżeli (a) dodatkowe zderzenie jest niezależne od procesu twardego, (b) obie próbki mają ten sam rozkład μ ważony świetlnością (te same LB i BCID, te same preskale w czasie)
i (c) wybór zdarzenia nie zależy od dodatkowego zderzenia. Nie znosi się:
1. człon f — zależy od krotności (0,25 % przy nSelBL 10–20, 1,2 % przy 40–50, 2–4 % powyżej 65). Jeżeli zdarzenia tagowane mają mniejszą krotność niż inkluzywne,
   ich ε_veto jest wyższe o ułamek punktu do ~1 punktu. Do policzenia z tabeli f(nSelBL), gdy znany będzie rozkład nSelBL w regionie tagowanym;
2. warunek (c) — dla tagu przypadkowego proton pochodzi z dodatkowego zderzenia, więc takie zdarzenia z definicji mają pile-up i weto usuwa je z inną (znacznie niższą)
   wydajnością niż sygnał. Klasa „other” w próbce inkluzywnej (nachylenie 1,98, ε_veto 88,2 %) pokazuje, że wybór skorelowany z dodatkowym zderzeniem zmienia ε_veto o kilka punktów.

**Czego nie da się ustalić bez regionu tagowanego (propozycje, niepoliczone):**
- ε_veto zdarzeń tagowanych w funkcji μ i nSelBL oraz jego rozkład na część prawdziwą i przypadkową;
- wydajność weta dla tagu przypadkowego (ile dodatkowych zderzeń dających proton w oknie ξ zostawia zrekonstruowany wierzchołek) — potrzebny region tagowany
  z odwróconym wetem (npu ≥ 1) albo MC dyfrakcyjne z nakładaniem;
- rozkład nSelBL zdarzeń tagowanych (do punktu 1 powyżej) i udział klas trygera wśród tagowanych kandydatów;
- czy ξ_cal i weto są skorelowane (oba reagują na dodatkową aktywność).

## 5. Zastrzeżenia

- Selekcja kandydatów jest nasza (D⁰→Kπ, pT ≥ 3 GeV, L_xy ≥ 0,25 mm, klasyfikator H), nie głównego łańcucha (CharmNet, D⁺, D* z pT od 2 GeV). ε_veto jest wielkością zdarzeniową
  i słabo zależy od kandydata, ale D⁺ i zakresu 2–3 GeV nie sprawdziłem.
- Sygnał D* w danych: N_S ≈ 2800 (428770) i ≈ 950 (427929); błąd nie zawiera niepewności skali pasm (wagi SBS stałe w bootstrapie). W małych binach (N_S < 100) błędy są zaniżone.
- 427929: 458 plików lokalnie; kompletność zbioru niepotwierdzona (znana sprawa). Na ε_veto nie wpływa, o ile brakujące pliki nie są wybrane po μ.
- `actual_mu` = 0 w 19 (428770) i 97 (427929) zdarzeniach — odrzucone.
- W `bt_fd_nd.npz` tabela `ev_` (235 408 wierszy) nie zgadza się kluczem (file_idx, event) z tabelami kandydatów (zgodność nSelBL 24 %; numery zdarzeń 1…10000 powtarzają się
  między plikami). Zmienne zdarzeniowe wziąłem z wiersza kandydata, a tabeli `ev_` dla fd_nd nie używam. Inne analizy łączące fd_nd po `ev_event` mogą być tym dotknięte — do sprawdzenia osobno.
- Klasa „other” jest zdefiniowana tylko przez 11 zapisanych łańcuchów; jej dokładny skład (inne łańcuchy z listy 106) wymagałby `tj2_*.npz` i nie był badany.
- Spadek ε_veto z LB w 428770 (0,3 punktu) i różnica ε_v między runami (3,9σ) są niewyjaśnione.

## 6. Pliki

- `v1_veto.py` — analiza główna (~75 s); `v1_veto_lib.py` — te same funkcje bez pętli (dla `v2_addendum.py`); `v2_addendum.py` — klasa „other”, J12 ∪ MBTS, zestawienia.
- `v0_explore.py`, `v0b_explore.py`, `v0c_explore.py`, `v0d_explore.py` (+ `v0_explore.txt`, `v0c_explore.txt`) — rozpoznanie plików, tło wiązki, klucz zdarzeń MC.
- `out/tables.md` — wszystkie tabele; `out/veto.json` — liczby; `out/addendum.{md,json}`; `out_v1.log`.
- Uruchomienie: `PYTHONDONTWRITEBYTECODE=1 ~/Documents/charm/.venv_h5/bin/python v1_veto.py` (z tego katalogu). Wejście tylko do odczytu: `B_trig/runs/lxplus/{bt,trg}_*.npz`,
  `B_trig/ISO_origin/iso_common.py`, `B_trig/btrig_common.py`.
