# TRKJET2 — trzy testy z literatury deficytu izolacji D⁰ w trygerze dżetowym (projekt B)

Wygenerowane: `analyze.py`, wejście `runs/lxplus`, 2026-11-17 07:51. Liczby: `out/trkjet2.json`.

## 0. Werdykty

| test | werdykt |
|---|---|
| T1 stożki prostopadłe / losowe | **NIE piedestał — piedestał w J12 jest podniesiony, lecz nie tłumaczy deficytu; nadwyżka w stożku D⁰: udział piedestału nieokreślony** |
| T2 podpróbka nieobciążona trygerem | **NIE — deficyt zostaje także wtedy, gdy tryger tłumaczy dżet po przeciwnej stronie** |
| T3 zastępnik energii neutralnej (spekulatywne) | **f_n BEZ WPŁYWU — przeważenie w r_n nie zmienia SF w J12only** |
| łącznie: „deficyt = obciążenie otoczenia D⁰ usuwalne zmienną zdarzenia” | **ODRZUCONA: ani piedestał, ani podpróbka away-side, ani f_n nie zdejmują deficytu → ε z danych per strumień (T1 lit.) albo embedding (T7 lit.)** |

| kryterium | liczby | wynik |
|---|---|---|
| K1.1 SF_perp = ε(perp1 < 2 GeV) dane/MC | 9:MB: ε 0,904 / 0,889 → 1,017 ± 0,024 (fit); 8:J12only: ε 0,801 / 0,884 → 0,906 ± 0,021 (fit) | piedestał podniesiony w J12: TAK; kontrola 9:MB zgodna z 1: TAK |
| K1.2 udział piedestału F_ped (J12only) | nadwyżka dane − MC: stożek D⁰ 4,89 ± 1,79 GeV, stożek prostopadły 0,74 ± 0,15 GeV; F_ped = – | nieokreślone (mianownik < 3σ) |
| K1.3 domknięcie | odjęcie piedestału: „[odjęcie] pT”: 9:MB 1,025 ± 0,029, J12only 0,733 ± 0,029, χ²/ndf 87,0/4, cof. 0,00; najlepsze przeważenie: „pT×perp_avg×lead_pt”: 9:MB 0,974 ± 0,035, J12only 0,689 ± 0,036, χ²/ndf 95,9/4, cof. 0,08; odniesienie: „pT×lead_pt×nSelBL (ISO_origin)”: 9:MB 0,985 ± 0,034, J12only 0,761 ± 0,040, χ²/ndf 48,0/4, cof. 0,00 | NIE |
| K2.1 SF w podpróbce T(X, Y = 5; weto O) | punkt główny X = 15 GeV: N_sig 217 ± 22, n_MC 624, ε 0,135 ± 0,042 / 0,338 → SF_T 0,400 ± 0,125 (po pT: 0,604 ± 0,151); inkluzywny 0,455 ± 0,027; różnica -0,055 ± 0,128 | NIE |
| K2.2 trend SF_T(X) | punktów miarodajnych 6; ostatni − pierwszy: 0,7σ | rośnie: TAK (informacyjnie) |
| K2.3 bity ortogonalne | (i) 8:softMB ∧ TBP_J12=0: N_sig 28, SF 0,92 ± 0,24 (niemiarodajny); (ii) 8:softMB ∧ TBP_J12=1: N_sig 11, SF -0,30 ± 0,15 (niemiarodajny); 9:MB ∧ TBP_J12=1: N_sig 29, SF 0,33 ± 0,19 (niemiarodajny); ułamek J12only z TBP_MBTS_1: 0,996 ± 0,003 | (i) warunki runu w porządku: –; (ii) deficyt idzie za dżetem L1: – |
| K3.0 MC nadaje się do przeważania w r_n wprost | ⟨r_n⟩_MC / ⟨r_n⟩_9:MB = 21,86 | NIE |
| K3.1 ⟨r_n⟩ J12only vs 8:notJ12 | 0,1208 ± 0,0062 vs 0,0803 ± 0,0035; różnica 5,7σ | TAK |
| K3.2 SF po przeważeniu „pT×u (kwantyl r_n)” | J12only: 0,531 ± 0,028 (sam pT: 0,570 ± 0,029; ΔSF = -0,039); 9:MB: 0,975 ± 0,036 | zmienia: NIE; domyka: NIE |

## 1. Kryteria (zapisane przed obejrzeniem danych)

Zapisane 17.11.2026 przed uruchomieniem ekstrakcji `trkjet2` na pełnych danych (stałe w `analyze.py`, drukowane dosłownie). Wspólne:
cięcie izolacji ptcone40 < 2 GeV (gałąź ntupla), selekcja i klasy jak w TRKJET; ε z fitu Δm pass/fail przy N_sig(SBS) ≥ 150, inaczej SBS
(„słaby” N_sig < 150, „niemiarodajny” N_sig < 50 lub N_sig/σ < 3 — nie wchodzi do werdyktu). „Domknięcie” = |SF − 1| < 2σ JEDNOCZEŚNIE
w 8:J12only i 9:MB. Odniesienie: pT×lead_pt×nSelBL (TRKJET: J12only 0,761 ± 0,040). Schemat z udziałem cofnięć („cof.”) > 0,30 w J12only
nie liczy się jako domknięcie. Bieg na niepełnym wejściu albo N_sig < 150 w klasie kluczowej → „NIEROZSTRZYGNIĘTE”.

**T1 — stożki prostopadłe / losowe (piedestał czy aktywność skorelowana z D⁰).** Zmienne: `perp1` = Σ pT torów (pT > 0,2 GeV, ta sama definicja co
własny stożek D⁰) w JEDNYM stożku R = 0,4 o osi (η_D⁰, φ_D⁰ ± π/2), znak wg parzystości numeru zdarzenia; `perp_avg` = średnia z obu;
`rnd` = stożek o losowym φ (≥ 0,8 od D⁰ i od osi wiodącego dżetu torowego).
- K1.1 (czy tryger podnosi piedestał): SF_perp = ε(perp1 < 2 GeV; dane) / ε(MC przeważone w pT). TAK, gdy SF_perp(J12only) < 1 o > 3σ.
  Kontrola: |SF_perp(9:MB) − 1| < 2σ (inaczej MC źle opisuje zdarzenie podkładowe niezależnie od trygera — zaznaczyć).
- K1.2 (udział piedestału w nadwyżce): F_ped = [⟨perp_avg⟩_dane − ⟨perp_avg⟩_MC] / [⟨stożek D⁰⟩_dane − ⟨stożek D⁰⟩_MC] w J12only
  (średnie sygnału, MC przeważone w pT). F_ped ≥ 0,7 → piedestał; F_ped ≤ 0,3 → aktywność skorelowana z D⁰; pomiędzy → mieszane.
  Liczone tylko, gdy mianownik > 3σ.
- K1.3 (domknięcie): (a) SF po odjęciu piedestału — warunek (własny stożek − perp_avg) < 2 GeV w danych i MC (MC przeważone w pT);
  (b) SF nominalny po przeważeniu MC w pT×perp (schematy `T1_PERP`; zmienne nie używają torów stożka D⁰ → brak koła).
  TAK, gdy (a) lub najlepszy (najmniejsze χ² na 4 klasach rozłącznych) schemat (b) domyka. „Częściowo”, gdy reszta |SF − 1| w J12only
  spada co najmniej 2× względem odniesienia przy 9:MB nadal w 2σ.
- Werdykt T1: K1.3 TAK → „PIEDESTAŁ tłumaczy deficyt”; K1.3 częściowo → „CZĘŚCIOWO piedestał”; K1.3 NIE → „NIE piedestał” z dopiskiem
  wg K1.1/K1.2 (piedestał podniesiony, lecz nie tłumaczy / piedestał niepodniesiony → nadwyżka skorelowana z D⁰).

**T2 — podpróbka nieobciążona trygerem.** Topologia T(X, Y): max pT dżetu torowego (kol. B) w |Δφ(dżet, D⁰)| > 2π/3 większy od X
oraz brak dżetu pT > Y w półkuli D⁰ (|Δφ| < π/2). Weto w wersji **O** (dżety z torów poza stożkiem D⁰ — bez koła; wersja główna)
i **B** (z torami stożka — [koło], informacyjnie). MC z tym samym wymaganiem; n_MC raportowane; n_MC < 100 = „MC słabe”, n_MC < 25 → brak SF.
- Punkt główny: Y = 5 GeV (O), największe X ∈ {5, 8, 10, 12, 15} z N_sig(J12only) ≥ 150 i n_MC ≥ 100. Brak takiego X → „NIEROZSTRZYGNIĘTE”.
- K2.1: SF_T = ε_fit(dane, T) / ε_MC(T, bez wag) w punkcie głównym; porównanie z SF inkluzywnym J12only (ta sama metoda).
  TAK („otoczenie D⁰ nieobciążone”), gdy |SF_T − 1| < 2σ oraz SF_T − SF_inkl > 2σ. CZĘŚCIOWO, gdy SF_T − SF_inkl > 2σ, ale SF_T < 1 o > 2σ.
  NIE, gdy SF_T − SF_inkl ≤ 2σ (deficyt zostaje także wtedy, gdy tryger tłumaczy dżet po przeciwnej stronie).
- K2.2 (informacyjnie): SF_T w funkcji X przy Y = 5 — czy rośnie z X; ε(dane) w binach pT dżetu strony przeciwnej.
- K2.3 (bity ortogonalne; każda komórka liczy się przy N_sig ≥ 150): (i) 428770, zdarzenia z HLT minbias bez wymagania twardości
  (lista `SOFT_MB`) i trigTBP_L1_J12 = 0 — SF zgodny z 1 (|pull| < 2) → warunki runu nie są przyczyną; (ii) zdarzenia zapisane trygerem
  minbias (428770 SOFT_MB albo 9:MB), lecz z trigTBP_L1_J12 = 1 — SF < 1 o > 2σ → deficyt idzie za samą obecnością dżetu L1, nie za tym,
  który tryger zapisał zdarzenie. Uwaga zapisana z góry: trigTBP_L1_MBTS_* jest ustawiony w prawie każdym zdarzeniu z D* — sam w sobie
  NIE definiuje próbki nieobciążonej (zdarzenie zapisał J12); raportowany tylko ułamek.

**T3 — zastępnik energii neutralnej (SPEKULATYWNE).** r_n = ξ_trk/ξ_cal, ξ_trk = Σ pT·(e^{+η} + e^{−η})/√s po wszystkich torach drzewa,
ξ_cal = ξ_cal_A + ξ_cal_C (suma stron — niezależna od przypisania A/C); f_n = 1 − r_n.
- K3.0 (czy MC nadaje się do przeważania wprost): ⟨r_n⟩_MC / ⟨r_n⟩_9:MB ∈ [0,67; 1,5]. NIE → schemat „pT×r_n” jest nieważny, decyduje
  schemat kwantylowy „pT×u” (u = dystrybuanta r_n w klasie odniesienia tego samego runu: 9:MB dla 427929, 8:notJ12 dla 428770; MC: własna).
- K3.1: ⟨r_n⟩ sygnału w J12only różni się od 8:notJ12 (ten sam run, te same warunki pile-upu) o > 3σ.
- K3.2: ΔSF = SF(schemat decydujący) − SF(pT) w J12only. „Zmienia”, gdy |ΔSF| ≥ 0,05 i |ΔSF| > 2σ; „domyka” wg definicji wspólnej.
  Dodatkowo dane→dane w tym samym runie: ε(J12only) / ε przewidziane z 8:notJ12 w komórkach (pT, r_n) względem samych komórek pT.
- Werdykt T3: domyka → „f_n TŁUMACZY (spekulatywne — wymaga potwierdzenia)”; zmienia, nie domyka → „f_n MA WPŁYW, nie domyka”;
  inaczej → „f_n BEZ WPŁYWU”.

**Łącznie:** hipoteza „deficyt = obciążenie otoczenia D⁰ przez tryger dżetowy, usuwalne zmienną zdarzenia” jest potwierdzona, gdy domyka
T1 albo T3, i wsparta, gdy T2 = TAK. Jeżeli T1 = NIE, T2 = NIE i T3 bez wpływu: deficyt siedzi w samym stożku D⁰ przy dowolnie dobranym
otoczeniu → nie da się go zdjąć przeważeniem MC minbias; pozostaje ε z danych per strumień (RECON_lit_B2 T1) albo embedding (T7).

## 2. Wejście i metoda

| próbka | par | plików | komplet córek | stożek losowy znaleziony | ξ_cal brak |
|---|---|---|---|---|---|
| data_427929 | 402470 | 458 | 0,9902 | 1,0000 | 0,0000 |
| data_428770 | 1233580 | 97 | 0,9917 | 1,0000 | 0,0000 |
| fd_nd | 178721 | 24 | 0,9959 | 1,0000 | 0,0000 |
| mb_nd | 207021 | 78 | 0,9922 | 1,0000 | 0,0000 |

- Pary, selekcja, klasy, fit i SBS jak w TRKJET (`../TRKJET/RAPORT.md` §1). MC: fd_nd, truth ∈ {0,1}, N = 41236; „MC przeważone w pT” = wagi w 12 binach pT do widma sygnału danej klasy.
- Stożki prostopadłe i losowy: te same tory i ta sama definicja co własny stożek D⁰ (`v_s40_p2`: wszystkie tory drzewa `Tracks`, pT > 0,2 GeV, bez córek D⁰; zgodność własnego stożka z gałęzią `ptcone40`: ta sama decyzja < 2 GeV w ≥ 99 %, TRKJET §1). Stożek losowy: pierwsza z 16 liczb losowych φ spełniająca warunki; gdy żadna — użyty stożek prostopadły.
- Punkt wyjścia (odtworzenie):

| klasa | N_sig (SBS) | ε SBS | ε fit | ε_MC | SF | jakość |
|---|---|---|---|---|---|---|
| 9:MB | 957 ± 46 | 0,655 ± 0,024 | 0,649 ± 0,028 | 0,664 | 0,98 ± 0,04 (fit) | ok |
| 8:all | 2796 ± 83 | 0,324 ± 0,013 | 0,338 ± 0,016 | 0,664 | 0,51 ± 0,02 (fit) | ok |
| 8:J12only | 2194 ± 73 | 0,293 ± 0,015 | 0,302 ± 0,018 | 0,664 | 0,46 ± 0,03 (fit) | ok |
| 8:MBTSbased | 369 ± 30 | 0,404 ± 0,040 | 0,415 ± 0,045 | 0,664 | 0,63 ± 0,07 (fit) | ok |
| 8:other | 238 ± 23 | 0,480 ± 0,048 | 0,532 ± 0,060 | 0,664 | 0,80 ± 0,09 (fit) | ok |
| 8:notJ12 | 562 ± 36 | 0,444 ± 0,032 | 0,473 ± 0,038 | 0,664 | 0,71 ± 0,06 (fit) | ok |

## 3. T1 — stożki prostopadłe / losowe

**3.1 Średnia Σ pT sygnału [GeV]: stożek D⁰, stożki prostopadłe, losowy (MC przeważone w pT danej klasy)**

| klasa | stożek D⁰: dane | MC | perp_avg: dane | MC | perp +π/2 | perp −π/2 | losowy: dane | MC | nadwyżka stożek D⁰ | nadwyżka perp | F_ped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 9:MB | 1,99 ± 0,18 | 2,13 | 0,77 ± 0,06 | 0,78 | 0,75 | 0,79 | 0,90 ± 0,09 | 0,76 | -0,15 ± 0,18 | -0,01 ± 0,06 | – |
| 8:all | 7,54 ± 1,42 | 3,83 | 1,50 ± 0,12 | 0,82 | 1,57 | 1,44 | 1,26 ± 0,08 | 0,88 | 3,70 ± 1,43 | 0,69 ± 0,12 | – |
| 8:J12only | 9,01 ± 1,78 | 4,12 | 1,56 ± 0,15 | 0,82 | 1,59 | 1,53 | 1,30 ± 0,10 | 0,91 | 4,89 ± 1,79 | 0,74 ± 0,15 | – |
| 8:MBTSbased | 1,99 ± 1,61 | 2,96 | 1,35 ± 0,16 | 0,80 | 1,68 | 1,02 | 1,13 ± 0,17 | 0,81 | -0,97 ± 1,61 | 0,55 ± 0,16 | – |
| 8:other | 2,47 ± 0,57 | 2,55 | 1,19 ± 0,13 | 0,79 | 1,21 | 1,17 | 1,03 ± 0,17 | 0,78 | -0,08 ± 0,57 | 0,40 ± 0,14 | – |
| 8:notJ12 | 1,84 ± 1,04 | 2,66 | 1,15 ± 0,09 | 0,79 | 1,29 | 1,00 | 1,08 ± 0,11 | 0,79 | -0,81 ± 1,04 | 0,35 ± 0,09 | – |

**3.2 ε(perp1 < próg) — „izolacja” stożka prostopadłego, dane vs MC (przeważone w pT)**

| klasa | próg [GeV] | N_sig | ε dane | metoda | ε MC | SF_perp | jakość |
|---|---|---|---|---|---|---|---|
| 9:MB | 1 | 956 ± 46 | 0,725 ± 0,023 | SBS | 0,750 | 0,967 ± 0,030 | ok |
| 9:MB | 2 | 956 ± 46 | 0,904 ± 0,021 | fit | 0,889 | 1,017 ± 0,024 | ok |
| 8:all | 1 | 2797 ± 83 | 0,639 ± 0,015 | SBS | 0,748 | 0,854 ± 0,020 | ok |
| 8:all | 2 | 2798 ± 83 | 0,798 ± 0,015 | fit | 0,885 | 0,902 ± 0,018 | ok |
| 8:J12only | 1 | 2193 ± 73 | 0,642 ± 0,017 | SBS | 0,747 | 0,860 ± 0,023 | ok |
| 8:J12only | 2 | 2194 ± 73 | 0,801 ± 0,018 | fit | 0,884 | 0,906 ± 0,021 | ok |
| 8:MBTSbased | 1 | 369 ± 30 | 0,611 ± 0,040 | SBS | 0,750 | 0,815 ± 0,054 | ok |
| 8:MBTSbased | 2 | 369 ± 30 | 0,775 ± 0,040 | fit | 0,888 | 0,873 ± 0,045 | ok |
| 8:other | 1 | 238 ± 23 | 0,642 ± 0,047 | SBS | 0,751 | 0,855 ± 0,063 | ok |
| 8:other | 2 | 238 ± 23 | 0,796 ± 0,048 | fit | 0,890 | 0,895 ± 0,054 | ok |
| 8:notJ12 | 1 | 562 ± 36 | 0,624 ± 0,032 | SBS | 0,751 | 0,831 ± 0,042 | ok |
| 8:notJ12 | 2 | 562 ± 36 | 0,780 ± 0,032 | fit | 0,889 | 0,877 ± 0,036 | ok |

**3.3 SF izolacji: nominalna, własny stożek, po odjęciu piedestału (fit pass/fail; MC przeważone w pT)**

| definicja pass | klasa | N_sig | ε dane | metoda | ε MC (pT) | ε MC (bez wag) | SF (MC pT) | jakość |
|---|---|---|---|---|---|---|---|---|
| nominalna (gałąź) < 2 GeV | 9:MB | 957 ± 46 | 0,649 ± 0,028 | fit | 0,670 | 0,664 | 0,969 ± 0,043 | ok |
| własny stożek < 2 GeV | 9:MB | 957 ± 46 | 0,658 ± 0,028 | fit | 0,673 | 0,667 | 0,978 ± 0,042 | ok |
| własny − perp_avg < 2 GeV | 9:MB | 956 ± 46 | 0,750 ± 0,025 | fit | 0,756 | 0,751 | 0,992 ± 0,033 | ok |
| nominalna (gałąź) < 2 GeV | 8:all | 2796 ± 83 | 0,338 ± 0,016 | fit | 0,512 | 0,664 | 0,660 ± 0,033 | ok |
| własny stożek < 2 GeV | 8:all | 2796 ± 83 | 0,342 ± 0,017 | fit | 0,515 | 0,667 | 0,665 ± 0,033 | ok |
| własny − perp_avg < 2 GeV | 8:all | 2797 ± 83 | 0,492 ± 0,018 | fit | 0,604 | 0,751 | 0,814 ± 0,031 | ok |
| nominalna (gałąź) < 2 GeV | 8:J12only | 2194 ± 73 | 0,302 ± 0,018 | fit | 0,489 | 0,664 | 0,619 ± 0,038 | ok |
| własny stożek < 2 GeV | 8:J12only | 2194 ± 73 | 0,306 ± 0,018 | fit | 0,491 | 0,667 | 0,623 ± 0,038 | ok |
| własny − perp_avg < 2 GeV | 8:J12only | 2193 ± 73 | 0,444 ± 0,021 | fit | 0,582 | 0,751 | 0,763 ± 0,036 | ok |
| nominalna (gałąź) < 2 GeV | 8:MBTSbased | 369 ± 30 | 0,415 ± 0,045 | fit | 0,582 | 0,664 | 0,713 ± 0,077 | ok |
| własny stożek < 2 GeV | 8:MBTSbased | 369 ± 30 | 0,432 ± 0,045 | fit | 0,585 | 0,667 | 0,737 ± 0,078 | ok |
| własny − perp_avg < 2 GeV | 8:MBTSbased | 368 ± 30 | 0,638 ± 0,048 | fit | 0,674 | 0,751 | 0,947 ± 0,072 | ok |
| nominalna (gałąź) < 2 GeV | 8:other | 238 ± 23 | 0,532 ± 0,060 | fit | 0,620 | 0,664 | 0,859 ± 0,096 | ok |
| własny stożek < 2 GeV | 8:other | 238 ± 23 | 0,519 ± 0,059 | fit | 0,623 | 0,667 | 0,833 ± 0,095 | ok |
| własny − perp_avg < 2 GeV | 8:other | 237 ± 23 | 0,689 ± 0,055 | fit | 0,710 | 0,751 | 0,971 ± 0,078 | ok |
| nominalna (gałąź) < 2 GeV | 8:notJ12 | 562 ± 36 | 0,473 ± 0,038 | fit | 0,609 | 0,664 | 0,777 ± 0,063 | ok |
| własny stożek < 2 GeV | 8:notJ12 | 562 ± 36 | 0,480 ± 0,038 | fit | 0,612 | 0,667 | 0,783 ± 0,063 | ok |
| własny − perp_avg < 2 GeV | 8:notJ12 | 561 ± 36 | 0,659 ± 0,039 | fit | 0,700 | 0,751 | 0,941 ± 0,055 | ok |

**3.4 SF po przeważeniu MC (izolacja nominalna; „[odjęcie]” = warunek własny − perp_avg < 2 GeV)**

| schemat | 9:MB | 8:all | 8:J12only | 8:MBTSbased | 8:other | 8:notJ12 | χ²/ndf vs 1 (4 klasy rozłączne) | cof. J12only |
|---|---|---|---|---|---|---|---|---|
| brak | 0,986 ± 0,036 | 0,488 ± 0,020 | 0,441 ± 0,022 | 0,609 ± 0,060 | 0,724 ± 0,073 | 0,669 ± 0,048 | 703,5/4 | 0,00 |
| pT | 0,977 ± 0,036 | 0,632 ± 0,027 | 0,600 ± 0,031 | 0,695 ± 0,070 | 0,776 ± 0,078 | 0,729 ± 0,053 | 197,0/4 | 0,00 |
| pT×lead_pt×nSelBL (ISO_origin) | 0,985 ± 0,034 | 0,774 ± 0,033 | 0,761 ± 0,040 | 0,786 ± 0,077 | 0,837 ± 0,081 | 0,798 ± 0,056 | 48,0/4 | 0,00 |
| pT×perp_avg | 0,977 ± 0,035 | 0,623 ± 0,026 | 0,585 ± 0,030 | 0,703 ± 0,070 | 0,792 ± 0,079 | 0,743 ± 0,054 | 220,6/4 | 0,00 |
| pT×perp_max×perp_min | 0,978 ± 0,036 | 0,620 ± 0,026 | 0,583 ± 0,030 | 0,696 ± 0,069 | 0,788 ± 0,079 | 0,737 ± 0,053 | 226,3/4 | 0,03 |
| pT×perp_avg×nSelBL | 0,988 ± 0,034 | 0,691 ± 0,029 | 0,661 ± 0,033 | 0,750 ± 0,073 | 0,805 ± 0,078 | 0,776 ± 0,055 | 120,2/4 | 0,02 |
| pT×perp_avg×lead_pt | 0,974 ± 0,035 | 0,709 ± 0,030 | 0,689 ± 0,036 | 0,728 ± 0,072 | 0,809 ± 0,080 | 0,759 ± 0,055 | 95,9/4 | 0,08 |
| pT×stożek losowy | 0,983 ± 0,036 | 0,617 ± 0,026 | 0,580 ± 0,029 | 0,694 ± 0,070 | 0,782 ± 0,078 | 0,734 ± 0,053 | 232,3/4 | 0,00 |
| [odjęcie] pT | 1,025 ± 0,029 | 0,787 ± 0,025 | 0,733 ± 0,029 | 0,951 ± 0,064 | 0,938 ± 0,068 | 0,932 ± 0,047 | 87,0/4 | 0,00 |
| [odjęcie] pT×lead_pt×nSelBL | 1,025 ± 0,028 | 0,882 ± 0,029 | 0,838 ± 0,034 | 1,013 ± 0,067 | 0,977 ± 0,070 | 0,972 ± 0,049 | 23,4/4 | 0,00 |

\* N_sig < 150. Klasy rozłączne do χ²: 9:MB, 8:J12only, 8:MBTSbased, 8:other.

**3.5 Dane→dane: ε(klasa) / ε przewidziane z 9:MB w komórkach (bez MC)**

| komórki | klasa | ε dane | ε przewidziane z 9:MB | stosunek | cof. |
|---|---|---|---|---|---|
| pT | 8:all | 0,324 | 0,543 | 0,596 ± 0,037 | 0,00 |
| pT×perp_avg | 8:all | 0,324 | 0,539 | 0,600 ± 0,036 | 0,28 |
| pT×r_n | 8:all | 0,324 | 0,522 | 0,620 ± 0,037 | 0,17 |
| pT | 8:J12only | 0,293 | 0,529 | 0,554 ± 0,038 | 0,00 |
| pT×perp_avg | 8:J12only | 0,293 | 0,528 | 0,555 ± 0,038 | 0,30 |
| pT×r_n | 8:J12only | 0,293 | 0,506 | 0,579 ± 0,040 | 0,20 |
| pT | 8:MBTSbased | 0,404 | 0,582 | 0,695 ± 0,074 | 0,00 |
| pT×perp_avg | 8:MBTSbased | 0,404 | 0,569 | 0,711 ± 0,077 | 0,21 |
| pT×r_n | 8:MBTSbased | 0,404 | 0,567 | 0,714 ± 0,077 | 0,07 |
| pT | 8:other | 0,480 | 0,612 | 0,785 ± 0,084 | 0,00 |
| pT×perp_avg | 8:other | 0,480 | 0,594 | 0,809 ± 0,088 | 0,19 |
| pT×r_n | 8:other | 0,480 | 0,600 | 0,800 ± 0,085 | 0,08 |
| pT | 8:notJ12 | 0,444 | 0,601 | 0,739 ± 0,060 | 0,00 |
| pT×perp_avg | 8:notJ12 | 0,444 | 0,592 | 0,751 ± 0,063 | 0,20 |
| pT×r_n | 8:notJ12 | 0,444 | 0,589 | 0,754 ± 0,061 | 0,07 |

**3.6 ε_iso w binach perp_avg [GeV]**

| próbka | 0–0.3 | 0.3–0.6 | 0.6–1 | 1–1.5 | 1.5–2.2 | 2.2–3.2 | 3.2–5 | ≥ 5 |
|---|---|---|---|---|---|---|---|---|
| 9:MB | 0,73 ± 0,03 (378) | 0,70 ± 0,06 (162) | 0,58 ± 0,06 (171) | 0,58 ± 0,08 (108) | 0,79 ± 0,11 (72) | 0,30 ± 0,13 (48) | – | – |
| 8:J12only | 0,29 ± 0,03 (607) | 0,27 ± 0,03 (370) | 0,28 ± 0,04 (288) | 0,34 ± 0,04 (333) | 0,25 ± 0,05 (204) | 0,31 ± 0,06 (166) | 0,31 ± 0,07 (119) | 0,36 ± 0,07 (114) |
| 8:MBTSbased | 0,61 ± 0,09 (79) | 0,47 ± 0,09 (65) | 0,47 ± 0,09 (64) | 0,45 ± 0,09 (67) | 0,18 ± 0,14 (35) | 0,12 ± 0,15 (31) | 0,18 ± 0,13 (27) | – |
| 8:other | 0,59 ± 0,08 (75) | – | 0,49 ± 0,13 (34) | 0,41 ± 0,13 (35) | – | – | – | – |
| MC fd_nd | 0,71 ± 0,00 (16826) | 0,67 ± 0,01 (7279) | 0,65 ± 0,01 (6263) | 0,62 ± 0,01 (4535) | 0,62 ± 0,01 (3155) | 0,57 ± 0,01 (1900) | 0,56 ± 0,02 (918) | 0,57 ± 0,03 (360) |

χ²/ndf: 9:MB vs 8:J12only: 207,4/6 (p = 5.1e-42); 9:MB vs MC: 8,7/6 (p = 0.19); 8:J12only vs MC: 641,4/8 (p = 2.9e-133); 8:J12only vs 8:MBTSbased: 23,2/7 (p = 0.0015).

## 4. T2 — podpróbka nieobciążona trygerem

**4.1 Skan T(X, Y)** — X: min pT dżetu torowego w |Δφ| > 2π/3 (0 = bez wymagania); Y: weto na dżet w półkuli D⁰ („-” = bez weta). SF = ε(dane)/ε_MC(T) bez wag; „SF po pT” = MC przeważone w pT wewnątrz T.

| weto | X | Y | klasa | N_sig | ε dane | metoda | n_MC | ε MC | SF | SF po pT | uwagi |
|---|---|---|---|---|---|---|---|---|---|---|---|
| O | 0 | - | 8:J12only | 2189 ± 73 | 0,303 ± 0,018 | fit | 40926 | 0,662 | 0,457 ± 0,028 | 0,573 ± 0,029 |  |
| O | 0 | - | 8:all | 2791 ± 83 | 0,338 ± 0,016 | fit | 40926 | 0,662 | 0,511 ± 0,025 | 0,609 ± 0,026 |  |
| O | 0 | - | 9:MB | 951 ± 45 | 0,646 ± 0,028 | fit | 40926 | 0,662 | 0,975 ± 0,043 | 0,976 ± 0,036 |  |
| O | 0 | 5 | 8:J12only | 651 ± 36 | 0,211 ± 0,026 | fit | 29201 | 0,696 | 0,304 ± 0,037 | 0,450 ± 0,048 |  |
| O | 0 | 5 | 8:all | 951 ± 43 | 0,314 ± 0,024 | fit | 29201 | 0,696 | 0,452 ± 0,034 | 0,584 ± 0,040 |  |
| O | 0 | 5 | 9:MB | 684 ± 36 | 0,711 ± 0,029 | fit | 29201 | 0,696 | 1,022 ± 0,042 | 1,033 ± 0,037 |  |
| O | 0 | 3 | 8:J12only | 274 ± 21 | 0,197 ± 0,034 | fit | 17160 | 0,732 | 0,269 ± 0,047 | 0,391 ± 0,061 |  |
| O | 0 | 3 | 8:all | 407 ± 26 | 0,333 ± 0,033 | fit | 17160 | 0,732 | 0,455 ± 0,045 | 0,582 ± 0,054 |  |
| O | 0 | 3 | 9:MB | 378 ± 25 | 0,785 ± 0,030 | fit | 17160 | 0,732 | 1,072 ± 0,041 | 1,065 ± 0,040 |  |
| O | 5 | - | 8:J12only | 1762 ± 67 | 0,279 ± 0,020 | fit | 13772 | 0,534 | 0,523 ± 0,037 | 0,627 ± 0,039 |  |
| O | 5 | - | 8:all | 2114 ± 73 | 0,295 ± 0,018 | fit | 13772 | 0,534 | 0,553 ± 0,035 | 0,636 ± 0,035 |  |
| O | 5 | - | 9:MB | 352 ± 29 | 0,408 ± 0,050 | fit | 13772 | 0,534 | 0,763 ± 0,093 | 0,836 ± 0,077 |  |
| O | 5 | 5 | 8:J12only | 484 ± 31 | 0,173 ± 0,028 | fit | 7779 | 0,560 | 0,308 ± 0,051 | 0,483 ± 0,067 |  |
| O | 5 | 5 | 8:all | 644 ± 36 | 0,228 ± 0,026 | fit | 7779 | 0,560 | 0,407 ± 0,047 | 0,555 ± 0,058 |  |
| O | 5 | 5 | 9:MB | 206 ± 20 | 0,443 ± 0,053 | fit | 7779 | 0,560 | 0,791 ± 0,094 | 0,868 ± 0,087 |  |
| O | 5 | 3 | 8:J12only | 193 ± 18 | 0,168 ± 0,038 | fit | 3442 | 0,573 | 0,293 ± 0,067 | 0,472 ± 0,092 |  |
| O | 5 | 3 | 8:all | 263 ± 21 | 0,252 ± 0,037 | fit | 3442 | 0,573 | 0,441 ± 0,065 | 0,590 ± 0,082 |  |
| O | 5 | 3 | 9:MB | 85 ± 11 | 0,622 ± 0,067 | SBS | 3442 | 0,573 | 1,086 ± 0,117 | 1,120 ± 0,118 | słaby |
| O | 8 | - | 8:J12only | 1482 ± 61 | 0,271 ± 0,021 | fit | 6062 | 0,460 | 0,590 ± 0,047 | 0,662 ± 0,047 |  |
| O | 8 | - | 8:all | 1676 ± 65 | 0,290 ± 0,020 | fit | 6062 | 0,460 | 0,631 ± 0,045 | 0,684 ± 0,045 |  |
| O | 8 | - | 9:MB | 159 ± 19 | 0,310 ± 0,069 | fit | 6062 | 0,460 | 0,674 ± 0,149 | 0,665 ± 0,125 |  |
| O | 8 | 5 | 8:J12only | 418 ± 28 | 0,174 ± 0,029 | fit | 3021 | 0,485 | 0,360 ± 0,061 | 0,498 ± 0,075 |  |
| O | 8 | 5 | 8:all | 511 ± 31 | 0,230 ± 0,029 | fit | 3021 | 0,485 | 0,475 ± 0,060 | 0,607 ± 0,071 |  |
| O | 8 | 5 | 9:MB | 88 ± 12 | 0,292 ± 0,065 | SBS | 3021 | 0,485 | 0,603 ± 0,134 | 0,619 ± 0,136 | słaby |
| O | 8 | 3 | 8:J12only | 164 ± 17 | 0,194 ± 0,041 | fit | 1292 | 0,485 | 0,401 ± 0,085 | 0,549 ± 0,113 |  |
| O | 8 | 3 | 8:all | 206 ± 18 | 0,237 ± 0,040 | fit | 1292 | 0,485 | 0,490 ± 0,083 | 0,613 ± 0,101 |  |
| O | 8 | 3 | 9:MB | 32 ± 7 | 0,388 ± 0,103 | SBS | 1292 | 0,485 | 0,800 ± 0,213 | 0,845 ± 0,217 | niemiarodajny |
| O | 10 | - | 8:J12only | 1312 ± 57 | 0,253 ± 0,022 | fit | 3847 | 0,420 | 0,603 ± 0,054 | 0,656 ± 0,052 |  |
| O | 10 | - | 8:all | 1456 ± 61 | 0,266 ± 0,021 | fit | 3847 | 0,420 | 0,634 ± 0,052 | 0,667 ± 0,050 |  |
| O | 10 | - | 9:MB | 97 ± 15 | 0,327 ± 0,073 | SBS | 3847 | 0,420 | 0,778 ± 0,174 | 0,785 ± 0,174 | słaby |
| O | 10 | 5 | 8:J12only | 353 ± 26 | 0,163 ± 0,032 | fit | 1791 | 0,442 | 0,368 ± 0,073 | 0,496 ± 0,086 |  |
| O | 10 | 5 | 8:all | 419 ± 29 | 0,207 ± 0,031 | fit | 1791 | 0,442 | 0,469 ± 0,072 | 0,586 ± 0,083 |  |
| O | 10 | 5 | 9:MB | 52 ± 9 | 0,279 ± 0,083 | SBS | 1791 | 0,442 | 0,630 ± 0,189 | 0,632 ± 0,182 | słaby |
| O | 10 | 3 | 8:J12only | 138 ± 15 | 0,170 ± 0,038 | SBS | 776 | 0,441 | 0,385 ± 0,087 | 0,522 ± 0,127 | słaby |
| O | 10 | 3 | 8:all | 161 ± 17 | 0,207 ± 0,043 | fit | 776 | 0,441 | 0,471 ± 0,100 | 0,587 ± 0,121 |  |
| O | 10 | 3 | 9:MB | 19 ± 5 | 0,328 ± 0,137 | SBS | 776 | 0,441 | 0,745 ± 0,313 | 0,863 ± 0,339 | niemiarodajny |
| O | 12 | - | 8:J12only | 1145 ± 54 | 0,243 ± 0,024 | fit | 2568 | 0,372 | 0,655 ± 0,067 | 0,701 ± 0,062 |  |
| O | 12 | - | 8:all | 1246 ± 57 | 0,258 ± 0,023 | fit | 2568 | 0,372 | 0,694 ± 0,065 | 0,728 ± 0,061 |  |
| O | 12 | - | 9:MB | 71 ± 12 | 0,321 ± 0,081 | SBS | 2568 | 0,372 | 0,864 ± 0,218 | 0,853 ± 0,212 | słaby |
| O | 12 | 5 | 8:J12only | 292 ± 24 | 0,142 ± 0,035 | fit | 1093 | 0,391 | 0,364 ± 0,090 | 0,498 ± 0,106 |  |
| O | 12 | 5 | 8:all | 331 ± 26 | 0,163 ± 0,034 | fit | 1093 | 0,391 | 0,417 ± 0,088 | 0,554 ± 0,106 |  |
| O | 12 | 5 | 9:MB | 29 ± 7 | 0,214 ± 0,111 | SBS | 1093 | 0,391 | 0,548 ± 0,286 | 0,574 ± 0,279 | niemiarodajny |
| O | 12 | 3 | 8:J12only | 116 ± 14 | 0,141 ± 0,039 | SBS | 461 | 0,371 | 0,381 ± 0,107 | 0,474 ± 0,141 | słaby |
| O | 12 | 3 | 8:all | 133 ± 15 | 0,160 ± 0,039 | SBS | 461 | 0,371 | 0,431 ± 0,109 | 0,533 ± 0,141 | słaby |
| O | 12 | 3 | 9:MB | 12 ± 4 | 0,340 ± 0,166 | SBS | 461 | 0,371 | 0,916 ± 0,451 | 0,908 ± 0,401 | niemiarodajny |
| O | 15 | - | 8:J12only | 970 ± 50 | 0,225 ± 0,026 | fit | 1617 | 0,332 | 0,678 ± 0,081 | 0,734 ± 0,072 |  |
| O | 15 | - | 8:all | 1040 ± 51 | 0,236 ± 0,025 | fit | 1617 | 0,332 | 0,711 ± 0,080 | 0,744 ± 0,071 |  |
| O | 15 | - | 9:MB | 43 ± 10 | 0,418 ± 0,105 | SBS | 1617 | 0,332 | 1,260 ± 0,320 | 1,186 ± 0,304 | niemiarodajny |
| O | 15 | 5 | 8:J12only | 217 ± 22 | 0,135 ± 0,042 | fit | 624 | 0,338 | 0,400 ± 0,125 | 0,604 ± 0,151 |  |
| O | 15 | 5 | 8:all | 247 ± 23 | 0,137 ± 0,038 | fit | 624 | 0,338 | 0,406 ± 0,115 | 0,597 ± 0,135 |  |
| O | 15 | 5 | 9:MB | 18 ± 5 | 0,245 ± 0,128 | SBS | 624 | 0,338 | 0,725 ± 0,382 | 0,733 ± 0,348 | niemiarodajny |
| O | 15 | 3 | 8:J12only | 95 ± 13 | 0,144 ± 0,042 | SBS | 256 | 0,328 | 0,438 ± 0,134 | 0,495 ± 0,159 | słaby |
| O | 15 | 3 | 8:all | 103 ± 13 | 0,144 ± 0,041 | SBS | 256 | 0,328 | 0,440 ± 0,132 | 0,491 ± 0,152 | słaby |
| O | 15 | 3 | 9:MB | 10 ± 4 | 0,389 ± 0,169 | SBS | 256 | 0,328 | 1,186 ± 0,524 | 0,985 ± 0,385 | niemiarodajny |
| B [koło] | 0 | - | 8:J12only | 2189 ± 73 | 0,303 ± 0,018 | fit | 40926 | 0,662 | 0,457 ± 0,028 | 0,573 ± 0,029 |  |
| B [koło] | 0 | 5 | 8:J12only | 323 ± 22 | 0,435 ± 0,042 | fit | 27162 | 0,747 | 0,583 ± 0,056 | 0,714 ± 0,059 |  |
| B [koło] | 0 | 3 | 8:J12only | 97 ± 11 | 0,529 ± 0,059 | SBS | 15027 | 0,832 | 0,635 ± 0,070 | 0,801 ± 0,085 | słaby |
| B [koło] | 5 | - | 8:J12only | 1762 ± 67 | 0,279 ± 0,020 | fit | 13772 | 0,534 | 0,523 ± 0,037 | 0,627 ± 0,039 |  |
| B [koło] | 5 | 5 | 8:J12only | 215 ± 18 | 0,412 ± 0,055 | fit | 6688 | 0,650 | 0,634 ± 0,084 | 0,781 ± 0,084 |  |
| B [koło] | 5 | 3 | 8:J12only | 62 ± 9 | 0,559 ± 0,072 | SBS | 2599 | 0,753 | 0,742 ± 0,097 | 0,954 ± 0,122 | słaby |
| B [koło] | 8 | - | 8:J12only | 1482 ± 61 | 0,271 ± 0,021 | fit | 6062 | 0,460 | 0,590 ± 0,047 | 0,662 ± 0,047 |  |
| B [koło] | 8 | 5 | 8:J12only | 177 ± 16 | 0,428 ± 0,057 | fit | 2416 | 0,604 | 0,709 ± 0,095 | 0,805 ± 0,095 |  |
| B [koło] | 8 | 3 | 8:J12only | 53 ± 8 | 0,572 ± 0,075 | SBS | 878 | 0,707 | 0,808 ± 0,107 | 1,000 ± 0,130 | słaby |
| B [koło] | 10 | - | 8:J12only | 1312 ± 57 | 0,253 ± 0,022 | fit | 3847 | 0,420 | 0,603 ± 0,054 | 0,656 ± 0,052 |  |
| B [koło] | 10 | 5 | 8:J12only | 145 ± 15 | 0,373 ± 0,050 | SBS | 1363 | 0,580 | 0,644 ± 0,087 | 0,782 ± 0,107 | słaby |
| B [koło] | 10 | 3 | 8:J12only | 40 ± 7 | 0,586 ± 0,087 | SBS | 487 | 0,700 | 0,837 ± 0,127 | 1,016 ± 0,154 | niemiarodajny |
| B [koło] | 12 | - | 8:J12only | 1145 ± 54 | 0,243 ± 0,024 | fit | 2568 | 0,372 | 0,655 ± 0,067 | 0,701 ± 0,062 |  |
| B [koło] | 12 | 5 | 8:J12only | 117 ± 13 | 0,342 ± 0,056 | SBS | 778 | 0,546 | 0,625 ± 0,105 | 0,745 ± 0,127 | słaby |
| B [koło] | 12 | 3 | 8:J12only | 33 ± 6 | 0,497 ± 0,098 | SBS | 257 | 0,661 | 0,752 ± 0,152 | 0,893 ± 0,182 | niemiarodajny |
| B [koło] | 15 | - | 8:J12only | 970 ± 50 | 0,225 ± 0,026 | fit | 1617 | 0,332 | 0,678 ± 0,081 | 0,734 ± 0,072 |  |
| B [koło] | 15 | 5 | 8:J12only | 89 ± 11 | 0,359 ± 0,064 | SBS | 413 | 0,508 | 0,707 ± 0,130 | 0,822 ± 0,154 | słaby |
| B [koło] | 15 | 3 | 8:J12only | 27 ± 6 | 0,502 ± 0,105 | SBS | 129 | 0,643 | 0,780 ± 0,171 | 0,874 ± 0,187 | niemiarodajny |

**4.2 ε_iso w binach pT dżetu strony przeciwnej [GeV] (weto O, Y = 5 GeV)**

| próbka | 0–3 | 3–5 | 5–8 | 8–12 | 12–20 | ≥ 20 |
|---|---|---|---|---|---|---|
| 8:J12only | 0,29 ± 0,07 (76) | 0,27 ± 0,06 (93) | 0,21 ± 0,08 (64) | 0,22 ± 0,04 (128) | 0,12 ± 0,04 (160) | 0,15 ± 0,04 (131) |
| 8:all | 0,49 ± 0,05 (146) | 0,43 ± 0,05 (163) | 0,22 ± 0,06 (131) | 0,32 ± 0,04 (183) | 0,15 ± 0,04 (185) | 0,15 ± 0,04 (146) |
| 9:MB | 0,81 ± 0,03 (334) | 0,87 ± 0,07 (151) | 0,62 ± 0,07 (118) | 0,33 ± 0,08 (60) | – | – |
| MC fd_nd | 0,78 ± 0,00 (13501) | 0,69 ± 0,01 (8221) | 0,61 ± 0,01 (4758) | 0,54 ± 0,01 (1928) | 0,42 ± 0,02 (828) | 0,30 ± 0,03 (265) |

**4.3 Bity ortogonalne** (kolumny `x_tbp_L1_J12`: 428770 jest, 427929 jest; `x_tbp_L1_MBTS_1`: jest). softMB = `HLT_noalg_mb_L1MBTS_1`, `HLT_noalg_mb_L1MBTS_1_1`, `HLT_noalg_mb_L1MBTS_2`, `HLT_mb_mbts_L1MBTS_1`, `HLT_mb_mbts_L1MBTS_1_1`, `HLT_mb_mbts_L1MBTS_2`, `HLT_mb_sptrk_L1MBTS_1`, `HLT_mb_sptrk_L1MBTS_2`, `HLT_noalg_L1MBTS_2`, `HLT_mb_sptrk_L1RD0_FILLED`, `HLT_noalg_L1RD0_FILLED`, `HLT_noalg_mb_L1RD0_FILLED`, `HLT_mb_mbts_L1RD0_FILLED`, `HLT_mb_sptrk_pt2_L1MBTS_2`, `HLT_mb_sptrk_pt2_L1RD0_FILLED`. SF względem MC inkluzywnego (bez wag) i przeważonego w pT.

| klasa | N_sig | ε dane | metoda | SF (bez wag) | SF po pT | jakość |
|---|---|---|---|---|---|---|
| 8:J12only ∧ TBP_MBTS_1 | 2184 ± 73 | 0,304 ± 0,018 | fit | 0,46 ± 0,03 | 0,60 ± 0,03 | ok |
| 8:softMB (niezależnie od J12) | 39 ± 10 | 0,392 ± 0,118 | SBS | 0,59 ± 0,18 | 0,59 ± 0,19 | niemiarodajny |
| 8:softMB ∧ TBP_J12=0 | 28 ± 9 | 0,610 ± 0,160 | SBS | 0,92 ± 0,24 | 0,97 ± 0,28 | niemiarodajny |
| 8:softMB ∧ TBP_J12=1 | 11 ± 4 | -0,202 ± 0,099 | SBS | -0,30 ± 0,15 | -0,27 ± 0,13 | niemiarodajny |
| 8:sptrk_pt≥4 ∧ TBP_J12=0 | 278 ± 25 | 0,423 ± 0,051 | fit | 0,64 ± 0,08 | 0,69 ± 0,08 | ok |
| 8:sptrk_pt≥4 ∧ TBP_J12=1 | 51 ± 11 | 0,263 ± 0,095 | SBS | 0,40 ± 0,14 | 0,54 ± 0,20 | słaby |
| 8:all ∧ TBP_J12=0 | 562 ± 36 | 0,473 ± 0,038 | fit | 0,71 ± 0,06 | 0,73 ± 0,05 | ok |
| 8:all ∧ TBP_J12=1 | 2237 ± 74 | 0,304 ± 0,018 | fit | 0,46 ± 0,03 | 0,60 ± 0,03 | ok |
| 9:MB ∧ TBP_J12=0 | 928 ± 45 | 0,657 ± 0,028 | fit | 0,99 ± 0,04 | 0,99 ± 0,04 | ok |
| 9:MB ∧ TBP_J12=1 | 29 ± 9 | 0,220 ± 0,127 | SBS | 0,33 ± 0,19 | 0,54 ± 0,32 | niemiarodajny |
| 9:MB ∧ sptrk_pt≥4 ∧ TBP_J12=0 | 481 ± 31 | 0,529 ± 0,037 | fit | 0,80 ± 0,06 | 0,89 ± 0,05 | ok |

## 5. T3 — zastępnik energii neutralnej (spekulatywne)

**5.1 Średnie sygnału** (r_n = ξ_trk/ξ_cal; ξ_cal = A + C; ξ_trk z torów pT > 0,2 GeV)

| próbka | ułamek z ξ_cal > 0 | ⟨ξ_cal⟩ | ⟨ξ_trk⟩ | ⟨r_n⟩ |
|---|---|---|---|---|
| 9:MB | 1,000 | 0,3017 ± 0,0080 | 0,0196 ± 0,0008 | 0,0684 ± 0,0017 |
| 8:all | 1,000 | 0,3860 ± 0,0055 | 0,0379 ± 0,0015 | 0,1126 ± 0,0049 |
| 8:J12only | 1,000 | 0,3937 ± 0,0062 | 0,0411 ± 0,0018 | 0,1208 ± 0,0062 |
| 8:MBTSbased | 1,000 | 0,3780 ± 0,0152 | 0,0275 ± 0,0016 | 0,0802 ± 0,0050 |
| 8:other | 1,000 | 0,3289 ± 0,0199 | 0,0240 ± 0,0015 | 0,0860 ± 0,0041 |
| 8:notJ12 | 1,000 | 0,3571 ± 0,0126 | 0,0252 ± 0,0011 | 0,0803 ± 0,0035 |
| MC fd_nd | 0,999 | 0,0251 ± 0,0001 | 0,0188 ± 0,0001 | 1,4962 ± 0,0128 |

**5.2 SF po przeważeniu MC** (u = kwantyl r_n względem klasy odniesienia tego samego runu; MC: własny)

| schemat | 9:MB | 8:all | 8:J12only | 8:MBTSbased | 8:other | 8:notJ12 | cof. J12only |
|---|---|---|---|---|---|---|---|
| pT | 0,978 ± 0,036 | 0,632 ± 0,027 | 0,600 ± 0,031 | 0,695 ± 0,070 | 0,776 ± 0,078 | 0,730 ± 0,053 | 0,00 |
| pT4 | 0,977 ± 0,036 | 0,606 ± 0,025 | 0,570 ± 0,029 | 0,682 ± 0,068 | 0,773 ± 0,077 | 0,723 ± 0,053 | 0,00 |
| pT×r_n (wprost) | 0,964 ± 0,048 | 0,608 ± 0,030 | 0,573 ± 0,033 | 0,677 ± 0,071 | 0,768 ± 0,081 | 0,718 ± 0,057 | 0,32 |
| pT×u (kwantyl r_n) | 0,975 ± 0,036 | 0,573 ± 0,025 | 0,531 ± 0,028 | 0,676 ± 0,068 | 0,758 ± 0,076 | 0,716 ± 0,052 | 0,00 |
| pT×u×nSelBL | 0,992 ± 0,035 | 0,659 ± 0,028 | 0,623 ± 0,032 | 0,743 ± 0,073 | 0,802 ± 0,078 | 0,774 ± 0,055 | 0,05 |
| pT×u×perp_avg | 0,976 ± 0,036 | 0,590 ± 0,025 | 0,548 ± 0,028 | 0,688 ± 0,069 | 0,777 ± 0,078 | 0,732 ± 0,053 | 0,05 |

**5.3 Dane→dane w 428770: ε(J12only) / ε przewidziane z 8:notJ12**

| komórki | ε J12only | ε przewidziane | stosunek | cof. |
|---|---|---|---|---|
| pT | 0,293 | 0,404 | 0,725 ± 0,062 | 0,00 |
| pT×r_n | 0,293 | 0,397 | 0,737 ± 0,067 | 0,00 |
| pT×perp_avg | 0,293 | 0,404 | 0,725 ± 0,060 | 0,24 |

**5.4 ε_iso w binach u (kwantyl r_n; małe u = mało torów względem kalorymetru = duża frakcja neutralna)**

| próbka | 0–0.2 | 0.2–0.4 | 0.4–0.6 | 0.6–0.8 | 0.8–1 | χ² stałej |
|---|---|---|---|---|---|---|
| 9:MB | 0,82 ± 0,04 (191) | 0,69 ± 0,05 (200) | 0,64 ± 0,05 (190) | 0,59 ± 0,06 (184) | 0,54 ± 0,05 (194) | 19,8/4 |
| 8:J12only | 0,38 ± 0,05 (206) | 0,31 ± 0,04 (348) | 0,34 ± 0,04 (361) | 0,26 ± 0,03 (402) | 0,26 ± 0,02 (880) | 8,3/4 |
| 8:notJ12 | 0,63 ± 0,07 (112) | 0,45 ± 0,08 (110) | 0,40 ± 0,08 (109) | 0,26 ± 0,06 (122) | 0,50 ± 0,07 (112) | 17,3/4 |
| MC fd_nd | 0,60 ± 0,01 (8241) | 0,61 ± 0,01 (8241) | 0,64 ± 0,01 (8242) | 0,68 ± 0,01 (8241) | 0,78 ± 0,00 (8242) | |

## 6. Zastrzeżenia

- Stożki prostopadłe mierzą aktywność przy tej samej η, 90° od D⁰; w zdarzeniach dwudżetowych (D⁰ w dżecie, drugi dżet naprzeciw) to obszar poprzeczny — piedestał „czysty”. Stożek losowy może trafić w stronę przeciwną (dżet odrzutu) → zwykle ≥ perp.
- T2: dżety torowe ≠ dżety kalorymetryczne L1; „twardy dżet torowy naprzeciw” nie gwarantuje, że to on odpalił L1_J12 (frakcja neutralna, próg ~12 GeV EM). MC minbias ma mało takich topologii — patrz n_MC; SF przy n_MC < 100 nie wchodzi do werdyktu.
- T2 weto B [koło]: brak dżetu przy D⁰ zbudowanego m.in. z torów stożka bezpośrednio podnosi ε — w danych i w MC, ale niekoniecznie tak samo.
- T2 K2.3: trigTBP_L1_MBTS_* nie wybiera próbki nieobciążonej (zdarzenie i tak zapisał J12); nieobciążone są tylko zdarzenia zapisane łańcuchem minbias.
- T3: ξ_cal w MC ma inną skalę niż w danych (test ROOT 17.11: mediana ξ_cal_A 0,0045 w fd_nd vs 0,05 w 427929 i 0,14 w 428770), a w 428770 zawiera pile-up → r_n nie jest czystą frakcją naładowaną; schemat kwantylowy zakłada, że ranga r_n w MC odpowiada randze w klasie odniesienia; 8:notJ12 nie jest próbką minbias. Wynik T3 w obie strony jest wskazówką, nie dowodem.
- Przeważanie: komórki z n_MC < 25 cofają się do mniejszej liczby wymiarów („cof.”); schemat z cof. > 0,30 w J12only nie jest testem.
- SBS zakłada kształt tła Δm niezależny od zmiennej wewnątrz grupy (bin pT × pass/fail).
