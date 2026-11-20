# MC_templates — test szablonów tła z MC w ficie TRExFitter głównego łańcucha (zadanie #7 planu prac z 17.11; projekt B)

Stan 17.11.2026 08:10 CET: **paczka gotowa i przetestowana, joby niewysłane, wyników fizycznych jeszcze nie ma.** Ten plik opisuje strukturę i kryteria;
liczby pojawią się w `runs/lxplus/ana/WYNIKI.md` (z condora) / `out/WYNIKI.md` (lokalnie) po etapie ORCH `mctmpl`.
Kryteria oceny: `KRYTERIA.md` — zapisane przed uruchomieniem analizy na jakichkolwiek histogramach (także testowych).

## 1. Pytanie
W ficie TRExFitter głównego łańcucha (`charmpp`, gałąź minimum-bias, `HistogramManager.cxx:555-556`, `minbias_trex_example.yaml:31-41`) swobodna jest tylko normalizacja `Comb`.
`Reflection` (D⁰), `CharmMismatch`, `DecayModeMismatch`, `PartComb` mają kształt i normalizację z MC, bez systematyk. Na ile te szablony zgadzają się z danymi
i jaką niepewność im przypisać (w % N_S)?

## 2. Co jest policzone inaczej niż w głównym łańcuchu (stałe zastrzeżenia)
- Klasyfikator: nasz H (XGB, 18 cech, bez izolacji), WP 95/80/50; nie CharmNet (ciaśniejsze WP, ptcone, dla D⁺ mKπ/m_D).
- Selekcja: D⁰ L_xy ≥ 0,25 mm (w głównym łańcuchu 0,10), D⁺ z wetem D* (Δm_min < 160 MeV; w głównym łańcuchu bez weta), pT 3–30 (bez binu 2–3 GeV), |η| < 2,2, okna masy i biny 10 MeV jak w głównym łańcuchu.
- Region: inkluzywny, GRL, **bez zdarzeń z tagiem AFP po cięciu ξ_cal i wecie wierzchołków** (odrzucane w ekstrakcji, nieliczone). Region sygnału głównego łańcucha nie jest oglądany;
  czy skład tła jest tam taki sam — nie wiemy.
- Normalizacja składowych: względem sygnału (κ_k = N_k/N_Sig z MC), nie absolutna σ·L — nie mamy L_eff per tryger. Równoważne konfiguracji głównego łańcucha tylko przy μ_sygnału = 1.
- MC: tylko ND (mb_nd 78 plików, fd_nd 24); bez SD/DD/CD. Bez podziału Sig/SigNonFid. Bez niepewności statystycznej MC w likelihoodzie.
- Kategorie: dokładnie mapowanie TruthMatchCategory → próbki z głównego łańcucha; Reflection wg definicji w charmpp (|pdgId hipotezy| ≠ |truthPdgId| którejś z dwóch córek).

## 3. Struktura wyniku (`WYNIKI.md`)
| § | treść | kryterium |
|---|---|---|
| 1 | Skład MC per kategoria (κ, podskład CharmMismatch/PartComb wg hadronu-rodzica) per pT i WP; stabilność κ i kształtu przy przeważeniu MC w (nSelBL, pT toru wiodącego) do pasm bocznych danych każdej klasy trygera | K1 |
| 2a | Refleksja K↔π: R/S i kształt z danych z tagu D* (bliźniak o zamienionej hipotezie) vs MC (ta sama procedura) i vs prawda | K2 |
| 2b | Rozpady częściowo zrekonstruowane / inne mody: μ̂ (dane/MC) z fitu z wolną normalizacją (nie jest to obszar niezależny) | K2 |
| 2c | D⁺: piki D* (Δm_min), φ (m_KK), D_s i D⁺→φπ (m_KKπ), Λ_c (m_pKπ): wielkość względem N_S i pozycja, dane vs MC | K2 |
| 2d | Udział z b: fit |impact| w oknie sygnału szablonami SigC/SigB vs prawda MC | K2 |
| 3 | Fity masy C0 („jak w głównym łańcuchu”) vs C0w (MC przeważone), C1–C3 (wolne μ składowych), C4 (refleksja z danych), C5, C6 (Comb exp(pol2)), C7 (wolny pik): ΔN_S, μ̂, σ; domknięcie na MC | K3, K4 |
| 4 | Proponowana systematyka per składowa i suma, w % N_S, per run/klasa trygera × gatunek × pT × WP | K3 |
| 5 | Co z tego wynika dla głównego łańcucha analizy | — |
Próbki danych: 427929 (wszystkie GRL; HLT_mb_mbts_L1MBTS_1), 428770 (wszystkie; J12only; MBTSbased; strumień AFP = któryś z łańcuchów L1AFP_A_OR_C).

## 4. Test (17.11, 1–4 pliki na próbkę; `runs/test/`)
Ekstrakcja, suma, DAG (35 + 1 węzłów), DRYRUN joba i analiza (lokalnie i w env LCG) przechodzą. Domknięcie kodu na MC: przy surowych szablonach z tego samego MC fit odtwarza N_Sig
dokładnie (pull 0,00 w komórkach bez dołączania rzadkich składowych); 18/18 komórek testowych w 1σ.
Zaobserwowane przy teście (własność metody, nie wynik fizyczny): wygładzanie szablonu PartComb funkcją exp(pol2) przesuwało N_S (PartComb ma część pikującą) → wygładzanie
PartComb usunięte przed jakimkolwiek spojrzeniem na fity danych; zostało tylko dla Comb z < 300 zdarzeń MC.
Liczby z plików testowych (1 plik danych, 2 pliki MC) są bez znaczenia statystycznego i nie są cytowane.

## 5. Znane ograniczenia metody (hipotezy do sprawdzenia na pełnej statystyce)
- 2a zakłada, że D⁰ z D* mają przy tym samym pT ten sam rozkład score bliźniaka co D⁰ inkluzywne; odjęcie tła Δm jedną funkcją progową i jednym α dla wszystkich binów masy.
- 2c: liniowe tło pod pikami; pik φ/D_s tylko dla części D_s→KKπ przez φ; r dotyczy sumy wszystkich źródeł piku, nie jednej kategorii prawdy.
- 2d: odjęcie pasm bocznych ze stałym współczynnikiem 2/3 (tło liniowe); cechy wierzchołkowe są w H.
- K1: przeważenie tylko w dwóch zmiennych zdarzeniowych; dla 428770 spodziewane małe n_eff (wtedy „nieokreślony”).
- Szablony po WP50 przy 8–30 GeV będą miały mało zdarzeń MC; fit nie zawiera niepewności statystycznej szablonów.
