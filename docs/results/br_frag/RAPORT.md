# BR_frag — zewnętrzne wejścia do przekroju fiducjalnego pojedynczych D: stosunki rozgałęzień i frakcje fragmentacji

Data: 2026-11-18. Źródła wyłącznie publiczne (PDG 2026, arXiv). Do sieci szły tylko ogólne zapytania i identyfikatory arXiv.
Lokalnie przeczytane (tylko odczyt): `GITLAB_2026-11-17/RAPORT.md` §2(d), §2(f); `RECON_program/RAPORT.md`; jedna liczba z `RECON_frag_mult/RAPORT.md`.

**Konwencja statusu**

- **[Ź]** — otwarte i przeczytane w pełnym tekście (PDF → `pdftotext`), liczba odczytana przeze mnie z tabeli/tekstu.
- **[S]** — tylko abstrakt / wynik wyszukiwarki / streszczenie strony.
- **[W]** — arytmetyka lub wniosek własny z liczb [Ź].
- **[niezweryfikowane]** — z pamięci, nieotwierane.

Wszystkie dostępy: 2026-11-18.

---

## 0. Najważniejsze

1. Przekrój fiducjalny na poziomie mezonu D potrzebuje z zewnątrz **tylko BR**: 0,76 % (D⁰), 1,50 % (D⁺), 1,06 % (łańcuch D*⁺). Frakcje fragmentacji nie wchodzą do wyniku, tylko do składu tła w MC, do porównań z teorią i do ewentualnej ekstrapolacji do σ(cc̄).
2. W R = σ_tag/σ_incl tego samego gatunku BR znosi się dokładnie. **Nie znosi się udział nieprompt (b→D)**, jeśli dyfrakcyjny ułamek b jest inny niż c (rachunki: 0,3–0,4 % wobec 1,6 %) — efekt na R rzędu −4…−12 % przy definicji „prompt + nieprompt” [W]. To wymaga decyzji o definicji sygnału.
3. Frakcje e⁺e⁻/ep (1509.01061) nie opisują pp przy niskim pT: ALICE widzi f(c→Λc) ≈ 17–19 % zamiast 6 %, a f(c→D⁰) ≈ 38–40 % zamiast 61 %. Stosunki mezon/mezon (D⁺/D⁰, D*⁺/D⁰) pozostają zgodne z e⁺e⁻.

---

## 1. Stosunki rozgałęzień

Źródło: PDG 2026 (F. Takahashi et al., Int. J. Mod. Phys. A 41, 2630011 (2026)), listingi z 1.06.2026:
`https://pdg.lbl.gov/2026/listings/rpp2026-list-D-zero.pdf`, `…-list-D-plus-minus.pdf`, `…-list-D-star-2010-plus-minus.pdf` [Ź];
te same wartości na `https://pdglive.lbl.gov/Particle.action?node=S032` (S031, M062) [S].

| Wielkość | Wartość | ± bezwzgl. | ± wzgl. | Uwagi | Status |
|---|---|---|---|---|---|
| B(D⁰→K⁻π⁺) | 3,936 % | 0,030 % | 0,76 % | „OUR FIT”, S = 1,3; fit 81 stosunków, 153 pomiary, χ² = 183,0/114. Średnia prosta: 3,909 ± 0,034 % | [Ź] |
| B(D⁺→K⁻π⁺π⁺) | 9,31 % | 0,14 % | 1,50 % | „OUR FIT”, S = 1,4; fit 41 stosunków, 57 pomiarów, χ² = 81,9/36. **Jedyny pomiar bezwzględny w ficie: CLEO-c** (BONVICINI 14: 9,224 ± 0,059 ± 0,157 %) | [Ź] |
| B(D*⁺→D⁰π⁺) | 67,7 % | 0,5 % | 0,74 % | fit 3 modów z więzem Σ = 1 (D⁺π⁰ 30,7 ± 0,5 %, D⁺γ 1,6 ± 0,4 %); korelacje: D⁰π⁺–D⁺π⁰ −62 %, D⁰π⁺–D⁺γ −43 %, D⁺π⁰–D⁺γ −44 %. Dominuje CLEO II 1998 | [Ź] |
| **B(D*⁺→D⁰π⁺) × B(D⁰→K⁻π⁺)** | **2,665 %** | **0,028 %** | **1,06 %** | iloczyn, składniki nieskorelowane (różne eksperymenty i fity) | [W] |
| B(D⁰→K⁺π⁻) (DCS + mieszanie) | 1,48 × 10⁻⁴ | 0,07 × 10⁻⁴ | 4,7 % | S = 3,1; = 0,376 % modu K⁻π⁺. R_D (samo DCS) = (3,45 ± 0,06) × 10⁻³ | [Ź] |
| B(D⁺→K⁺π⁺π⁻) | 4,88 × 10⁻⁴ | 0,08 × 10⁻⁴ | — | DCS; inny układ ładunków niż sygnał (piony przeciwnych znaków) | [Ź] |
| B(K⁻π⁺π⁺)/B(K⁻π⁺) | 2,365 | 0,040 | 1,68 % | przy braku korelacji; z korelacją ρ ≈ 0,3–0,5 → ≈ 1,2–1,4 % | [W] |

Pomiary wchodzące do B(D⁰→K⁻π⁺) [Ź]: BESIII 2018 (3,883 ± 0,006 ± 0,051 %), CLEO-c 2014 (3,934 ± 0,021 ± 0,061 %), BaBar 2008 (4,007 ± 0,037 ± 0,072 %), CLEO II (3,82 ± 0,07 ± 0,12 %), ALEPH, ARGUS.

### 1.1 Korelacje

- **D⁰ ↔ D*⁺:** ten sam B(D⁰→K⁻π⁺) w obu kanałach → jeden wspólny parametr, 100 % skorelowany. B(D*⁺→D⁰π⁺) osobny, nieskorelowany.
- **D⁰ ↔ D⁺:** PDG fituje D⁰ i D⁺ osobno i nie podaje korelacji między nimi. W CLEO-c (arXiv:1312.6775, tab. IV) współczynnik korelacji B(Kπ)–B(Kππ) wynosi **0,69** (wspólne tory, PID, N_DD̄) [Ź]. CLEO-c ma ~30 % wagi w średniej D⁰ i 100 % w D⁺, stąd efektywnie ρ ≈ 0,3–0,5 [W, zgrubnie]. Propozycja: nominalnie ρ = 0 (zachowawcze dla stosunku D⁺/D⁰), wariant ρ = 0,5 jako test.
- ALICE w stosunku D⁺/D⁰ przyjmuje ±1,8 % od BR (0,009/0,500), czyli brak korelacji [Ź, 2308.04877 tab. 5].

### 1.2 Promieniowanie końcowe (FSR)

- CLEO-c: „we include the effects of final-state radiation by allowing for additional unobserved photons in the final state”; FSR symulowane PHOTOS, systematyka FSR 0,2–0,8 % per mod [Ź, 1312.6775]. PDG nie komentuje FSR przy wartości z fitu; jedyna uwaga dotyczy starego AKERIB 93 („includes radiative corrections; without them 0,0391”) [Ź].
- Wniosek: wartości PDG to **B(D→Kπ(γ))**, inkluzywne w fotonach. Dla nas znaczy to tyle:
  - sygnał w MC musi obejmować rozpady z dodatkowymi fotonami z PHOTOS. Jeśli dopasowanie prawdy wymaga dokładnie dwóch (trzech) córek, część sygnału trafi do `DecayModeMismatch` i wydajność będzie zaniżona. **Do sprawdzenia w jego kodzie** (`src/DSelection.cxx`, etykiety `MatchedC`/`MatchedB`).
  - ogon radiacyjny w masie jest częścią szablonu sygnału; jego niepewność jest w PEAK_WIDTH / kształcie, nie w BR.
  - rozmiar: ułamek rozpadów z twardym fotonem to kilka procent [niezweryfikowane].

### 1.3 Rozpady podwójnie tłumione (D⁰→K⁺π⁻)

- PDG odjęło DCS od pomiaru BESIII (który mierzył sumę K∓π±), więc 3,936 % to **czyste K⁻π⁺** [Ź, przypis 1 w bloku danych].
- **D⁰ bez tagu D*** (nie znamy zapachu): pik zawiera D⁰→K⁻π⁺ i D̄⁰→K⁻π⁺ (DCS). Jeśli liczymy D⁰ + D̄⁰ razem, właściwy jest B(K⁻π⁺) + B(K⁺π⁻) = **3,951 %** (+0,38 %, połowa niepewności BR) [W]. Przy `keep_both` DCS przechodzi między szablonem sygnału a `Reflection` — musi być potraktowany tak samo w MC i w definicji. Sprawdzić, czy EvtGen w próbkach ma mod DCS.
- **D*⁺→D⁰π⁺ z tagiem ładunku pionu:** kombinacja „wrong sign” jest poza selekcją → czyste B(K⁻π⁺) = 3,936 %.
- **D⁺:** mod DCS ma inną konfigurację ładunków, nie wchodzi do sygnału.

---

## 2. Frakcje fragmentacji powabu

### 2(a) Średnie e⁺e⁻ / ep / (stare) pp — Lisovyi, Verbytskyi, Zenaiev

**arXiv:1509.01061**, Eur. Phys. J. C 76 (2016) 397 [Ź, PDF]. Tabela 12 (kombinacja globalna, więz S = 1) i tabela 13 (korelacje):

| f(c→H) | Wartość domyślna | Wariant (stałe σ(e⁺e⁻→c), Γ_cc̄/Γ_had) |
|---|---|---|
| D⁰ | 0,6086 ± 0,0076 | 0,6141 ± 0,0073 |
| D⁺ | 0,2404 ± 0,0067 | 0,2439 ± 0,0067 |
| D_s⁺ | 0,0802 ± 0,0040 | 0,0797 ± 0,0040 |
| Λ_c⁺ | 0,0623 ± 0,0041 | 0,0549 ± 0,0026 |
| D*⁺ | 0,2429 ± 0,0049 | 0,2386 ± 0,0046 |
| D*⁰ | 0,2306 ± 0,0315 | 0,2250 ± 0,0299 |
| D_s*⁺ | 0,0548 ± 0,0074 | 0,0537 ± 0,0072 |

Korelacje (tab. 13): D⁰–D⁺ −0,66; D⁰–D_s −0,32; D⁰–Λc −0,41; D⁺–D_s −0,19; D⁺–Λc −0,19; D_s–Λc −0,07; D*⁺–D⁰ +0,19; D*⁺–D⁺ −0,08; D*⁺–D_s −0,07; D*⁺–Λc −0,12.

- Więz: S = f(D⁰) + f(D⁺) + f(D_s) + (1 + λ) f(Λc) = 1, λ = 0,136 ± 0,006 (Ξc, Ωc ze skalowania sektora dziwnego LEP). D*⁺ **nie** wchodzi do sumy (rozpada się na D⁰/D⁺).
- Cztery gatunki słabo rozpadające się + jeden więz = 3 stopnie swobody — stąd „trzy wariacje wektorów własnych” w procedurze ATLAS (§3.3).
- BR użyte w tej pracy to PDG 2014/15: D⁰→Kπ 3,93 ± 0,04 %, D⁺→Kππ 9,46 ± 0,24 %, D*⁺ 67,7 ± 0,5 %, Λc→pKπ 6,84 % [Ź, tab. 1; wartość D*⁺ tylko ze streszczenia]. D⁺ przesunęło się od tego czasu o −1,6 %; frakcje nie są przeliczone na PDG 2026.
- Stosunki [W]: D⁺/D⁰ = 0,395; D*⁺/D⁰ = 0,399.

### 2(b) pp przy LHC — ALICE, |y| < 0,5, pT > 0

- 5,02 TeV: **arXiv:2105.06335**, Phys. Rev. D 105, L011103 [Ź, PDF; czasopismo S].
- 13 TeV: **arXiv:2308.04877**, JHEP 12 (2023) 086 [Ź, PDF]. Tabela 7 tej pracy podaje też zaktualizowane wartości 5,02 TeV — te cytuję.

| f(c→H) [%] | pp 5,02 TeV (zaktualizowane) | pp 13 TeV | e⁺e⁻/ep (§2a) | pp 13 TeV / e⁺e⁻ [W] |
|---|---|---|---|---|
| D⁰ | 39,6 ± 1,7 +2,6/−3,8 | **38,2 ± 1,3 +2,3/−4,3** | 60,9 ± 0,8 | 0,63 |
| D⁺ | 17,5 ± 1,8 +1,7/−2,1 | **19,1 ± 1,4 +1,5/−2,3** | 24,0 ± 0,7 | 0,79 |
| D_s⁺ | 7,4 ± 1,0 +1,9/−1,1 | **6,1 ± 0,5 +1,2/−0,9** | 8,0 ± 0,4 | 0,76 |
| Λ_c⁺ | 18,9 ± 1,3 +1,5/−2,0 | **16,8 ± 0,8 +1,5/−2,1** | 6,2 ± 0,4 | 2,7 |
| Ξ_c⁰ | 8,1 ± 1,2 +2,5/−2,5 | 9,9 ± 1,3 +2,3/−2,4 | (w λ) | — |
| Ξ_c⁺ | przyjęte = Ξ_c⁰ | 9,6 ± 1,2 +3,9/−4,8 | (w λ) | — |
| J/ψ | 0,44 ± 0,03 +0,04/−0,06 | 0,37 ± 0,02 +0,04/−0,05 | — | — |
| D*⁺ (poza sumą) | 15,7 ± 1,2 +4,1/−1,9 | **15,6 ± 0,7 +2,5/−2,2** | 24,3 ± 0,5 | 0,64 |
| Σ_c⁰,⁺,⁺⁺ (poza sumą) | — | 7,2 ± 1,2 +1,6/−1,9 | — | — |

(Pierwsza niepewność statystyczna, druga systematyczna z ekstrapolacją. Pierwotne wartości 5,02 TeV z 2105.06335: D⁰ 39,1; D⁺ 17,3; D_s 7,3; Λc 20,4; Ξc⁰ 8,0; D*⁺ 15,5 [Ź].)

- Stosunki mezonowe, pT > 0 [Ź, tab. 5]: D⁺/D⁰ = 0,442 ± 0,055 ± 0,033 ± 0,008(BR) przy 5,02 TeV; **0,500 ± 0,047 ± 0,033 ± 0,009(BR)** przy 13 TeV. D*⁺/D⁰ = 306/749 = 0,41 z tab. 4 [W]. Zgodne z e⁺e⁻ (0,395 i 0,399) w granicach ~1–2σ.
- ALICE: wzrost f(Λc) ~×3,3 względem e⁺e⁻/ep (istotność ~5σ), spadek mezonów D ×1,2–1,4; brak istotnej zależności od energii między 5,02 a 13 TeV [Ź].
- Korelacje: ALICE **nie podaje macierzy**. Opis [Ź]: tory i poprawka na frakcję prompt w pełni skorelowane między gatunkami, ekstrakcja sygnału i statystyka nieskorelowane, ekstrapolacja częściowo; Ωc jako asymetryczna systematyka mianownika. Suma = 1 z konstrukcji (normalizacja do sumy zmierzonych gatunków), więc frakcje są antyskorelowane, ale liczb nie ma.
- σ(cc̄), dσ/dy w |y| < 0,5: 5,02 TeV 1165 ± 44 +63/−67 +98/−38 ± 43(BR) ± 42(RS) ± 24(lumi) µb [Ź, 2105.06335]; 13 TeV 2031 ± 61 +135/−141 +196/−63 ± 97(BR) ± 33(lumi) ± 73(y) µb [centralna wartość i pierwsze składniki Ź, 2308.04877; składniki BR/lumi/y ze streszczenia].

### 2(c) Która wersja jest właściwa

- **Do pp przy niskim pT i centralnej rapidity: ALICE 13 TeV** (najbliżej 13,6 TeV). Wartości e⁺e⁻ dają f(c→D⁰) zawyżone o ~60 % względne, więc σ(cc̄) z samego D⁰ byłoby zaniżone ×0,63.
- Zastrzeżenia: (i) wartości ALICE są całkowane po pT od zera i dla |y| < 0,5; nasz obszar to pT 2–30 GeV, |η| < 2,2. Λc/D⁰ maleje z pT, więc przy wyższych pT frakcje zbliżają się do e⁺e⁻ [niezweryfikowane ilościowo]. Frakcje „właściwe dla binu” nie istnieją jako jedna tabela. (ii) Niepewności ALICE są 5–10× większe niż e⁺e⁻.
- Dla stosunków mezon/mezon wybór wersji ma małe znaczenie (zgodne w obu zestawach).
- Pythia 8 A3 w naszych próbkach: f(c→D*⁺) = 0,213 wobec 0,243 w e⁺e⁻ (`RECON_frag_mult/RAPORT.md` §rekomendacja — liczba lokalna, nie sprawdzana ponownie).

---

## 3. Jak to wchodzi do pomiaru

### 3.1 Kontekst: jego fit (z `GITLAB_2026-11-17/RAPORT.md` §2(d), (f))

- 12 regionów (tag × {D⁺, D*⁺, D⁰} × 4 biny pT), jeden POI na gatunek × bin pT prawdy; swobodna tylko normalizacja `Comb`; **wszystkie pozostałe tła z MC** (σ·L/N_gen): `Reflection`, `CharmMismatch`, `DecayModeMismatch`, `PartComb`, `SigNonFid`, `SigOutFid`.
- Systematyki BR i frakcji produkcji: „do dodania”, wzorzec w repozytorium `hf-production-fraction-uncertainties`. Lokalny raport nie mówi o nim nic więcej (bez zmian od ostatniego przeglądu). Przypuszczam, że to implementacja procedury z noty ATL-PHYS-PUB-2022-035 [W; nota nieotwierana].

### 3.2 Gdzie co wchodzi

| Wejście | σ_fid(D), poziom mezonu | σ(cc̄) / porównanie na poziomie kwarków | Przez MC (akceptancja, tła) |
|---|---|---|---|
| BR sygnału | **tak, 1:1** (σ = N/(ε·L·B)) | tak | nie (o ile MC ma ten sam mod wymuszony lub wagi BR) |
| BR innych modów powabu | nie | nie | tak: `DecayModeMismatch`, `CharmMismatch`, `PartComb` (normalizacja z MC) |
| f(c→H) | **nie** | **tak, 1:1** (σ_cc̄ = σ_D / (2 f)) | tak: skład `CharmMismatch` (D_s→KKπ, Λc→pKπ pod D⁺), udział D⁰ z D*⁺, pośrednio cechy klasyfikatora |
| Udział nieprompt | zależy od definicji (prompt vs prompt + b→D) | tak | tak: ε różne dla prompt i nieprompt (L_xy, d₀, izolacja) |

- Jeśli POI jest zdefiniowany jako μ = σ/σ_MC, to BR wchodzi przy przeliczeniu μ → σ, a BR w generatorze (EvtGen) musi być przeskalowany do PDG. Niepewność BR sygnału: `OverallSys` na `Sig*` danego gatunku albo dodana po ficie — równoważne, bo jest czysto multiplikatywna.

### 3.3 Standardowa procedura przeważania frakcji w MC

Wzorzec publiczny: ATLAS W+D, **arXiv:2302.00336** [Ź, PDF], odsyłający do ATL-PHYS-PUB-2022-035 (CDS 2816367; **nieotwierana**) i do 1509.01061.

1. Każde zdarzenie MC dostaje wagę = iloczyn po słabo rozpadających się hadronach powabowych f_ref(H)/f_MC(H), z f_ref = średnie światowe; f_MC liczone z prawdy generatora.
2. Niepewność: trzy wariacje wektorów własnych macierzy kowariancji frakcji (skorelowane, z więzem Σ = 1), każda jako osobny parametr uciążliwy z wagami zdarzeń.
3. BR innych modów: próbki z BR przesuniętymi jednocześnie w ramach niepewności PDG; zmiana normalizacji i kształtu tła propagowana do fitu.

Rozmiary w W+D [Ź]: zmiana pojedynczych frakcji mezonów do **20 %** zależnie od generatora; normalizacja szablonów tła do **3 %**; niepewność frakcji na tle ≤ 3 % (plus kształt masy); BR innych modów ≤ 5 % na tle; wpływ obu na obserwable „pomijalny”; BR sygnału 1,7 % (D⁺) i 1,1 % (D*⁺) na wydajku — zgodne z §1 (1,50 % i 1,06 % w PDG 2026).

U nas [W]:
- Sygnał: waga frakcji zmienia normalizację przewidywania MC, nie wydajność ε per gatunek → przy pomiarze σ_fid(D) znika. Zostaje efekt drugiego rzędu przez udział D⁰ z D*⁺ (jeśli ε zależy od pochodzenia).
- Tła pikujące i półpikujące z MC o stałej normalizacji: tu frakcje i BR wchodzą bezpośrednio. Przy przejściu e⁺e⁻ → pp wkład Λc rośnie ×2,7, D_s spada ×0,76. Efekt na σ = (zmiana frakcji) × (udział danego tła pod pikiem); przy udziale kilku procent to ≲ 1–3 %.
- Propozycja: nominał = wagi do ALICE 13 TeV; wariacja „dwupunktowa” do e⁺e⁻ (1509.01061) jako jeden parametr; plus wewnętrzne niepewności ALICE (uproszczone: D⁰, D⁺, D_s, Λc niezależnie i renormalizacja do Σ = 1, bo macierzy nie ma).

### 3.4 Co się znosi

| Obserwabla | BR sygnału | f(c→H) | Nieprompt | Uwagi |
|---|---|---|---|---|
| σ_fid(D) | nie (0,76 / 1,50 / 1,06 %) | tylko przez tła MC | zależy od definicji | — |
| **R = σ_tag(D)/σ_incl(D)**, ten sam gatunek | **znosi się dokładnie** | znosi się w pierwszym rzędzie; tła `CharmMismatch` w liczniku i mianowniku skalują się razem, o ile skład jest ten sam | **nie znosi się** | DCS, FSR też się znoszą |
| D⁺/D⁰ w jednej próbce | zostaje B(Kππ)/B(Kπ): 1,2–1,7 % | to jest mierzona wielkość | częściowo (różne f_np per gatunek, kilka %) | — |
| D*⁺/D⁰ w jednej próbce | B(D⁰→Kπ) znosi się; zostaje B(D*⁺→D⁰π⁺): 0,74 % | jw. | częściowo | — |
| Podwójny stosunek (D⁺/D⁰)_tag / (D⁺/D⁰)_incl | **znosi się dokładnie** | to jest test uniwersalności | resztkowo | — |

Nieprompt w R [W]: R_obs = R_c · [1 − f_np · (1 − R_b/R_c)]. Z R_b/R_c ≈ 0,2 (dipol: bb̄ 0,3–0,4 %, cc̄ 1,6 % — `RECON_program` B1) i f_np = 5–15 % daje to −4…−12 %. Liczby modelowe, nie pomiar; pokazują jednak, że definicja sygnału (tylko prompt albo prompt + nieprompt) musi być ustalona i taka sama w liczniku i mianowniku.

### 3.5 Udział nieprompt jako osobne wejście

- ALICE 13 TeV, |y| < 0,5, **arXiv:2402.16417** [Ź, PDF]: naturalny udział nieprompt „od 5 % przy niskim pT do 15 % przy wysokim pT”; prompt/nieprompt ≈ 20 przy niskim pT, stosunek nieprompt/prompt rośnie do pT ≈ 12 GeV. Wartości per bin tylko na rysunkach.
- U niego etykiety `MatchedC` / `MatchedB` istnieją w kodzie (raport §2(d)); z lokalnego raportu nie wynika, czy `Sig` zawiera oba. **Do ustalenia.**
- Opcje: (A) sygnał = prompt + nieprompt (definicja bez etykiety, spójna z fiducjalną definicją protonu), f_np z MC wchodzi tylko przez różnicę wydajności; (B) sygnał = prompt, f_np odejmowany z MC/FONLL — wtedy f_np jest pełnym parametrem uciążliwym. Przy cięciach na L_xy i d₀ wydajność nieprompt jest inna, więc w obu opcjach potrzebna wariacja f_np. Rozmiar wariacji (±30 % względne poniżej) to moja propozycja, nie liczba ze źródła.

---

## 4. Tabela parametrów uciążliwych

| Parametr uciążliwy | Wartość centralna | ±1σ | Źródło | Status | Znosi się w R? |
|---|---|---|---|---|---|
| `BR_D0_Kpi` (D⁰ i D*⁺, wspólny) | 0,03936 | ±0,00030 (0,76 %) | PDG 2026, listing D⁰, Γ47 | [Ź] | tak |
| `BR_Dp_Kpipi` | 0,0931 | ±0,0014 (1,50 %) | PDG 2026, listing D±, Γ68 | [Ź] | tak |
| `BR_Dst_D0pi` | 0,677 | ±0,005 (0,74 %) | PDG 2026, listing D*(2010)± | [Ź] | tak |
| (pochodny) łańcuch D*⁺→(K⁻π⁺)π⁺ | 0,02665 | ±0,00028 (1,06 %) | iloczyn dwóch powyższych | [W] | tak |
| (pochodny) D⁰ bez tagu, K⁻π⁺ + DCS | 0,03951 | ±0,00030 | PDG 2026, Γ47 + Γ344 | [W] | tak |
| korelacja `BR_D0_Kpi`–`BR_Dp_Kpipi` | 0 (wariant 0,5) | — | CLEO-c 1312.6775 tab. IV: 0,69 | [Ź] / [W] | — |
| `FF_c_D0` | 0,382 (pp) | +0,026/−0,045 | ALICE 2308.04877 tab. 7 | [Ź] | tak (1. rząd) |
| `FF_c_Dp` | 0,191 (pp) | +0,021/−0,027 | jw. | [Ź] | tak (1. rząd) |
| `FF_c_Dst` (poza sumą) | 0,156 (pp) | +0,026/−0,023 | jw. | [Ź] | tak (1. rząd) |
| `FF_c_Ds` | 0,061 (pp) | +0,013/−0,010 | jw. | [Ź] | tak (1. rząd) |
| `FF_c_Lc` | 0,168 (pp) | +0,017/−0,022 | jw. | [Ź] | tak (1. rząd) |
| `FF_model` (pp ↔ e⁺e⁻, dwupunktowy) | pp | e⁺e⁻: D⁰ 0,6086 ± 0,0076; D⁺ 0,2404 ± 0,0067; D*⁺ 0,2429 ± 0,0049; D_s 0,0802 ± 0,0040; Λc 0,0623 ± 0,0041 | 1509.01061 tab. 12, korelacje tab. 13 | [Ź] | tak (1. rząd) |
| `BR_other_charm` (tła MC) | BR z EvtGen → PDG | per mod z PDG; wzorzec W+D: ≤ 5 % na tle | 2302.00336 §9 | [Ź] (procedura) | w dużej części |
| `f_nonprompt` | z MC, per bin pT (rząd 5–15 %) | ±30 % względne (propozycja) | ALICE 2402.16417 (zakres); wariacja własna | [Ź] / [W] | **nie** |

Niepewności FF w tabeli: stat ⊕ syst w kwadraturze [W]. FF wchodzą tylko jako wagi MC (tła) oraz przy ekstrapolacji do σ(cc̄); do σ_fid(D) nie wchodzą multiplikatywnie.

---

## 5. Czego nie udało się zweryfikować

- **ATL-PHYS-PUB-2022-035** (procedura przeważania, wartości referencyjne, wektory własne) — nieotwierana (CDS poza zakresem źródeł); znam ją tylko z cytowania w 2302.00336. Zawartości `hf-production-fraction-uncertainties` nie znam.
- Rozmiar efektu FSR (ułamek rozpadów z fotonem powyżej progu) — [niezweryfikowane].
- Zależność frakcji od pT i rapidity w naszym obszarze (2–30 GeV, |η| < 2,2) — brak tabeli; trend Λc/D⁰(pT) z pamięci.
- Macierz korelacji frakcji ALICE — nie istnieje publicznie w tych pracach.
- Korelacja B(D⁰→Kπ)–B(D⁺→Kππ) na poziomie fitu PDG — nie jest podana; 0,3–0,5 to oszacowanie.
- Rozbicie niepewności dσ_cc̄/dy przy 13 TeV (BR, lumi, y) oraz czasopismo pracy 2105.06335 — odczytane tylko ze streszczenia [S].
- Strony pdglive ładują bloki pomiarów dynamicznie; wartości sprawdziłem w PDF-ach listingów.
- Nie sprawdzałem BR i frakcji faktycznie użytych w naszych próbkach MC (EvtGen, A3) ani obecności modu DCS i ustawień PHOTOS.

## 6. Odnośniki

- PDG 2026: https://pdg.lbl.gov/2026/listings/rpp2026-list-D-zero.pdf ; https://pdg.lbl.gov/2026/listings/rpp2026-list-D-plus-minus.pdf ; https://pdg.lbl.gov/2026/listings/rpp2026-list-D-star-2010-plus-minus.pdf ; https://pdglive.lbl.gov/Particle.action?node=S032 (S031, M062) — dostęp 2026-11-18.
- arXiv:1509.01061 (Lisovyi, Verbytskyi, Zenaiev, EPJC 76 (2016) 397) [Ź].
- arXiv:2105.06335 (ALICE, 5,02 TeV, PRD 105, L011103) [Ź].
- arXiv:2308.04877 (ALICE, 13 TeV, JHEP 12 (2023) 086) [Ź].
- arXiv:2402.16417 (ALICE, nieprompt D, 13 TeV) [Ź].
- arXiv:1312.6775 (CLEO-c, PRD 89, 072002) [Ź].
- arXiv:2302.00336 (ATLAS W+D, 13 TeV) [Ź].
