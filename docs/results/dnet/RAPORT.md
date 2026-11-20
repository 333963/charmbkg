# DNET — „odpowiednia sieć”: klasyfikator sieciowy dla D⁰→Kπ i D⁺→Kππ wobec drzewa XGB H

Projekt B (charmbkg), tylko pojedyncze D (bez par D–D). Projekt, uzasadnienie technologii i kryteria zapisane **17.11.2026 przed
wynikami pełnej statystyki** (dotąd tylko test mechaniki na 1 pliku h5 na próbkę, §8). Wyniki dopisuje `DNET/analyze.py report` w §10.

## 1. Pytanie i odpowiedź projektowa (krótko)

Pytanie: czy i która sieć jest lepsza od obecnego drzewa XGB „H” (18 cech, bez izolacji) dla tej analizy.

Rekomendowana technologia (do potwierdzenia liczbami z §10):

- **Model: mały MLP tabelaryczny** (rezydualny, 128 jednostek, 3 bloki, ~56 tys. wag) z przetwarzaniem wstępnym w grafie, a nie duża sieć na torach.
  Powód: ~10⁵ sygnału MC na bin pT i kilkanaście–kilkadziesiąt cech tabelarycznych to obszar, w którym drzewa i dobrze ustawione MLP są
  równoważne, a transformery tabelaryczne nie dają zysku (§3); sufit informacyjny tła dla cech H jest już osiągnięty przez XGB (nota §5.1).
- **Źródło zysku: cechy, nie model.** Sieć ma sens tylko dlatego, że **strata z dekorelacją od masy pozwala użyć kinematyki rozpadu**
  (cosθ* ze znakiem, pT córek, asymetria, ΔR, α), której drzewo nie może użyć bez rzeźbienia. Dlatego porównanie jest dwuosiowe: model
  (XGB ↔ MLP na tych samych cechach) i cechy (H ↔ H + kinematyka).
- **Dekorelacja: DisCo** (dCor² w stracie) jako wybór podstawowy, MoDe[1] jako wariant kontrolny; adwersarz i planing odrzucone (§2.1).
- **Bez izolacji w wariancie nominalnym**; izolacja bezwzględna i względna tylko jako warianty z osobno zmierzonym kosztem dane/MC.
- **Tło treningowe z MC (truth ≥ 8) w całym zakresie mas**; tło z pasm bocznych danych tylko w wariancie hybrydowym „sbmc”
  (BCE na SB danych, dekorelacja na tle MC) — powód w §2.3.
- **Zbiór torów (DeepSets / transformer wzorca GN2) — etap 2, warunkowy.** Otoczenie = wyuczona izolacja, czyli dokładnie ta informacja,
  która w 428770 nie przenosi się z MC na dane; „tylko córki” nie niesie informacji ponad cechy tabelaryczne poza korelacjami parametrów torów.

## 2. Wymagania → wybory

### 2.1 Brak rzeźbienia masy tła → dekorelacja w stracie (DisCo)

| metoda | za | przeciw | decyzja |
|---|---|---|---|
| **DisCo** (dCor² score–masa na tle, 2001.05310) | jeden hiperparametr λ, stabilny trening, jedna sieć; ta sama postać co w CharmNet (`learn.py`, dCor na tle kombinatorycznym) → wdrożenie bez zmian w charmpp; zero dCor ⇔ niezależność (także nieliniowa) | wymusza pełną niezależność, czyli więcej niż wymaga R_sculpt (który dopuszcza liniową ε_B(m)); koszt O(n²) na wsad → podpróbka 1024 wierszy tła | **wybrana** |
| **MoDe[1]** (2010.09745) | dopuszcza liniową zależność ε_B(m) — dokładnie to, co dopuszcza nasze kryterium R_sculpt i model tła fitu → potencjalnie mniejsza strata mocy | biny masy, progi i temperatura indykatora = dodatkowe ustawienia; nie ma jej w charmpy; nasza implementacja ma gładki indykator (nie oryginalny) | **wariant kontrolny** (1 konfiguracja) |
| adwersarz (1611.01046, 1703.03507) | najlepszy w części porównań | dwie sieci, gra min–max, niestabilny, wiele ustawień; przy naszej statystyce i wymaganiu ≥ 3 seedów zbyt kosztowny w strojeniu | odrzucony |
| planing / wagi spłaszczające masę | trywialny | usuwa tylko brzegową zależność od masy; nie działa, gdy masa jest funkcją kilku wejść (pT córek + ΔR + cosθ* odtwarzają m(Kπ)) (1908.08959) | odrzucony |
| uBoost / strata płaskości dla drzew (1305.7248, 1410.4140) | dekorelacja dla GBDT | nie ma jej w XGBoost; osobna biblioteka, brak eksportu ONNX w środowisku | odnotowana, nie testowana |

Uwaga do CharmNet: w gałęzi minimum-bias charmpy λ = 0,001 przy stracie dCor (nie dCor²) z podpróbką ≤ 1000 — człon jest rzędu 10⁻⁴ straty,
czyli praktycznie wyłączony (w innej wersji CharmNet wyłączony jawnie). U nas λ jest skanowane (1 / 5 / 25 / 100 przy dCor²), a nominał (5 i 25) ma 3 seedy.
Skala: dCor² niewytrenowanej sieci ≈ 0,006–0,01, różnice BCE między modelami ≈ 0,01–0,05 → λ rzędu 1–25; test mechaniki (§8) pokazał
zapaść do stałej przy λ = 400 i silną stratę przy 100, stąd zakres skanu.

### 2.2 Kalibrowalność z danych → wejścia istniejące dla D⁰←D*

Wszystkie wejścia etapu 1 są funkcjami kandydata D⁰ i jego dwóch córek, więc istnieją identycznie dla D⁰ z rozpadu D*; ε_S mierzona fitem
Δm pass/fail per run i klasa trygera (428770: all, J12only; 427929: MB427929) tym samym narzędziem co B_trig/CN_clone
(`btrig_common.fit_eps` → `charmtag.dstar_refit.fit_passfail`). Otoczenie (etap 2) też istnieje dla D⁰←D*, ale wolny pion jest w stożku
(ΔR małe) → uczyłby sieć „jestem z D*”; w ekstrakcji NIE jest usuwany (tak jak w ptcone gałęzi), co jest jawnym źródłem różnicy
D*-tag ↔ inkluzywne D⁰ (nota §5.8: transfer zawodzi dla cięcia bezwzględnego przy 8–30 GeV). Dla D⁺ nie ma tagu → SF nie do zmierzenia;
dla D⁺ rozstrzygają Zproxy, transfer na SB danych i rzeźbienie.

### 2.3 Odporność dane/MC → warianty izolacji i źródła tła

- Izolacja: `HK` (bez), `HKia` (ptcone40, ptcone100 [GeV], jak CharmNet), `HKir` (ptcone/pT). Znany koszt: twarde cięcie SF 0,51; sieć z izolacją
  R_SF 0,91–0,94; mechanizm: 56 % D⁰ ma dżet torowy w ΔR < 0,4 w J12only wobec 24 % w MC (nota §5.8).
- Tło treningowe: `mc` (mb truth ≥ 8, cały zakres mas), `sb` (dalekie SB danych 428770, jak nominalny XGB H), `sbmc` (BCE na SB danych,
  DisCo na tle MC). **Pułapka `sb` dla sieci z kinematyką:** w SB nie ma zdarzeń pod pikiem, więc sieć, która z pT córek, ΔR i cosθ* odtwarza
  masę, może nauczyć się „m poza oknem → tło” przy zerowym dCor na próbce SB. Dlatego `sb` dla zestawu HK jest tylko kontrolą (oczekiwane
  rzeźbienie na MC truth ≥ 8 w oknie), a wariantem roboczym z danymi jest `sbmc`.
- Ocena wyłącznie out-of-fold: 5 foldów po zdarzeniach (hash `event_uid`), ten sam podział dla wszystkich modeli i dla danych.
- Cech dE/dx ponad H (zestaw `HKP`) nie ma w nominale: FEAT_charmnet pokazał C2ST 0,915 dla dE/dx córek D⁺ i artefakt „MC vs dane”.

### 2.4 Statystyka i typ danych → mała sieć, kontrola pojemności

Dwie architektury na tych samych cechach: `m64` (3×64 ReLU = przepis CharmNet, 9,6 tys. wag dla 18 wejść) i `r128` (rezydualna, LayerNorm,
SiLU, 56 tys. wag). Nie zakładamy, że większa wygra: obie są w paczce, a zysk z pojemności odczytujemy z pary `m64.HK` ↔ `r128.HK`.
Ustawienia „dobre z góry” zamiast strojenia (RealMLP 2407.04491): odporne skalowanie (mediana/IQR) + asinh w grafie, AdamW, harmonogram
one-cycle, 30 epok, wybór epoki na 10 % zdarzeń walidacyjnych, wsad 4096. Limit próbki 4·10⁵ sygnału / 8·10⁵ tła (ten sam dla drzew).

### 2.5 Wdrożenie → ONNX

Eksport ONNX (opset 17) modelu trenowanego na całości, wejście `features` [N, n_cech] w kolejności z `meta/<tag>.json`, wyjście `score` ∈ (0, 1),
skalowanie i uzupełnianie braków wewnątrz grafu (żadnego zewnętrznego skalera), kontrola onnxruntime ↔ torch zapisana w meta. Różnica względem
CharmNet w charmpp: tam jeden graf zawiera `n_folds` sieci, a fold wybiera `event_number % n_folds` (drugie wejście). Nasze foldy to hash
`event_uid` (MC: (run, event) nie jest kluczem) — do produkcji model trzeba wytrenować ponownie z podziałem po `event_number` i wyeksportować
w układzie CharmNet (zmiana ~20 linii w `dn_models.export_onnx`). XGB nie ma w środowisku eksportera ONNX (brak onnxmltools) — to argument
praktyczny za siecią przy równej jakości.

## 3. Przegląd literatury (co z tego pasuje do nas)

| # | pozycja | treść | wniosek dla nas |
|---|---|---|---|
| 1 | ATLAS, „Transforming jet flavour tagging at ATLAS” (GN2), arXiv:2505.19689 | transformer na torach dżetu z zadaniami pomocniczymi (pochodzenie toru, wierzchołki); c-jet rejection ×3,5 względem DL1d | zysk pochodzi z dziesiątek torów na dżet, milionów dżetów treningowych i zadań pomocniczych; nasz obiekt to 2–3 tory z gotowym wierzchołkiem i 10⁵ sygnału → wzorzec tylko dla wariantu „otoczenie” (etap 2), bez zadań pomocniczych |
| 2 | Qu, Gouskos, ParticleNet, arXiv:1902.08570; Qu, Li, Qian, Particle Transformer, arXiv:2202.03772 | grafy / transformer na chmurze cząstek; ParT wymaga ~10⁸ dżetów do pełnej przewagi | za duże na naszą statystykę; uwaga parami (ΔR, m_ij) dla 2 torów = dokładnie cechy kinematyczne z zestawu HK |
| 3 | Zaheer i in., Deep Sets, arXiv:1703.06114; Komiske, Metodiev, Thaler, Energy Flow Networks, arXiv:1810.05165; ATLAS DIPS (ATL-PHYS-PUB-2020-014) | suma po elementach zbioru + sieć końcowa; tańsze i stabilniejsze od uwagi przy małej próbce | kandydat `ds` w etapie 2 obok małego transformera `tr` |
| 4 | Kasieczka, Shih, DisCo, arXiv:2001.05310; Kasieczka i in., ABCDisCo, arXiv:2007.14400 | dCor² w stracie, jakość jak adwersarz przy jednym hiperparametrze; ABCD z dwóch zdekorelowanych zmiennych | metoda podstawowa; ABCDisCo (dwie sieci zdekorelowane od siebie) nie jest potrzebne — tło bierzemy z fitu masy, nie z ABCD |
| 5 | Kitouni i in., MoDe, arXiv:2010.09745 | kontrola rzędu zależności ε_B(m): płaska, liniowa, kwadratowa | wariant kontrolny MoDe[1], bo nasz fit i R_sculpt dopuszczają liniowość |
| 6 | Bradshaw i in., „Mass agnostic jet taggers”, arXiv:1908.08959; Louppe i in., arXiv:1611.01046; Shimmin i in., arXiv:1703.03507 | porównanie planingu, PCA, adwersarza, uBoost: metody w stracie > metody na danych, gdy masa jest funkcją wejść | planing odrzucony; adwersarz odrzucony ze względu na koszt strojenia |
| 7 | Stevens, Williams, uBoost, arXiv:1305.7248; Rogozhnikov i in., arXiv:1410.4140 | jednorodność wydajności w masie dla drzew (praktyka LHCb) | pokazuje, że dekorelacja dla drzew istnieje, ale poza XGBoost i bez ONNX → nie testujemy |
| 8 | Grinsztajn, Oyallon, Varoquaux, arXiv:2207.08815; Shwartz-Ziv, Armon, arXiv:2106.03253; McElfresh i in., arXiv:2305.02997 | na średnich danych tabelarycznych (10⁴–10⁶ wierszy) GBDT ≥ sieci; różnice małe, sieci przegrywają na cechach nieinformatywnych i rozkładach z ogonami | nie oczekiwać zysku z samej zamiany XGB → MLP na cechach H (hipoteza B: ~0); skalowanie odporne + asinh adresuje ogony |
| 9 | Gorishniy i in., FT-Transformer, arXiv:2106.11959; Holzmüller i in., RealMLP, arXiv:2407.04491 | FT-Transformer ≈ ResNet-MLP ≈ GBDT po strojeniu; dobrze ustawiony MLP dorównuje GBDT bez strojenia | wybór: MLP rezydualny z ustawieniami domyślnymi; FT-Transformer nie testowany (koszt bez oczekiwanego zysku) |
| 10 | Hollmann i in., TabPFN, arXiv:2207.01848 (v2: Nature 2025) | sieć wstępnie wytrenowana, wnioskowanie w kontekście; limit ~10⁴ wierszy (v2: ~10⁵), bez własnej straty | odpada: brak członu dekorelacji, brak eksportu ONNX małego modelu, próbka za duża |
| 11 | García Pardiñas i in. (LHCb), DFEI, arXiv:2304.08610 | GNN na wszystkich torach zdarzenia do rekonstrukcji rozpadów ciężkich hadronów | kierunek dla selekcji D na poziomie zdarzenia; wymaga d0/z0 wszystkich torów — nasze drzewo `Tracks` ich nie ma |

Praktyka selekcji D/B w eksperymentach (ALICE, LHCb, CMS BPH): GBDT lub mały MLP na cechach wierzchołka; sieci na torach stosuje się do
tagowania dżetów, nie do ekskluzywnych kandydatów. Wniosek: dla nas właściwa jest sieć tabelaryczna z dekorelacją; architektury
zbiorowe są uzasadnione tylko, jeśli chcemy uczyć izolację z torów — a to jest informacja o znanym koszcie dane/MC.

## 4. Kryteria „sieć przyjęta zamiast drzewa” (zapisane 17.11.2026 przed wynikami pełnej statystyki)

Odniesienie: `D0.xgb.H.sb.n.s1` (XGB H, tło z dalekich SB danych — jak nominał) i odpowiednio `Dplus.xgb.H.sb.n.s1`. Punkt pracy WP70
(pomocniczo WP80, WP50), biny pT 3–5 / 5–8 / 8–30 GeV. G = Z/Z_H, Z = ε_S/√ε_B, błąd z bootstrapu sparowanego po zdarzeniach.

- **K1 (zysk).** G ≥ 1,10 w co najmniej 2 binach pT **jednocześnie** na MC (tło truth ≥ 8 w oknie) i po przeniesieniu na dane (ε_B z dalekich SB
  428770), z istotnością zysku na MC ≥ 2σ.
- **K2 (rzeźbienie).** |R_sculpt − 1| (MC truth ≥ 8 w oknie, WP70) nie większe niż dla H w granicach 2σ; kontrola na paśmie danych
  [okno, dalekie SB) i nachylenie ε_B(m) w granicach ±0,10/100 MeV od H (bramka G3 z DLZ).
- **K3 (kalibracja, tylko D⁰).** |SF − 1| w binach 5–8 i 8–30 nie większe niż dla H w granicach błędu, w każdej z klas: 428770 all, 428770
  J12only, 427929 MB427929. Brak pomiaru (fit nieokreślony) = nierozstrzygnięte, nie „spełnione”.
- **K4 (stabilność).** Odchylenie standardowe G między ≥ 3 seedami < 1/3 zysku (G − 1) w binach 3–5 i 5–8.
- **Werdykt.** Wszystkie cztery → sieć przyjęta. K1 + K2 spełnione, K3 lub K4 nierozstrzygnięte → warunkowo (potrzebna statystyka / seedy).
  K1 niespełnione → zostaje drzewo; wtedy jedynym argumentem za siecią jest wdrożenie (ONNX, zgodność z CharmNet) przy G ≥ 0,98.
- **Dodatkowo (bez wpływu na werdykt, do decyzji o etapie 2).** Etap 2 jest wart wysłania, jeśli w etapie 1: (a) `r128.HK` spełnia K1 i K2
  (jest co poprawiać architekturą) **albo** (b) wariant z izolacją (`HKia`/`HKir`) daje G_tr ≥ 1,10 ponad `HK` przy R_SF ≥ 0,95 (izolacja opłaca
  się mimo kosztu) — wtedy otoczenie z torów może dać więcej niż dwa stożki. W przeciwnym razie etap 2 nie ma uzasadnienia.
- Wariant z izolacją jest oceniany osobno: do K3 dochodzi R_SF = SF_iso/SF_noiso (kryterium K1 z CN_clone: koszt, jeśli < 0,90 przy ≥ 2σ).

## 5. Kandydaci i konfiguracje

Tag = `<gatunek>.<model>.<cechy>.<tło>.<dekorelacja><λ>.s<seed>`; definicje w `Q_lxplus/dnet1/dn_feat.py` (`SETS`, `BLOCKS`).

| zestaw cech (D⁰) | liczba | zawartość |
|---|---|---|
| H | 18 | odniesienie (`charmbkg.features`, bez izolacji) |
| HX | 29 | H + parametr zderzenia / czas: d0·d0, |d0|/σ (param.), Σ(d0/σ)², |Δz0|, |impact|, |z0| córek, cτ (= DLZ_d0 „NOKIN”) |
| HK | 37 | HX + kinematyka rozpadu: cosθ* ze znakiem, asymetria pT, min/max pT córek i ich stosunki do pT(D), ΔR(K, π), α Armenterosa |
| HKP | 44 | HK + dE/dx córek względem Bethe–Blocha (4 odległości, log p_π, npts_π, dE/dx_π) |
| HKia / HKir | 39 | HK + ptcone40, ptcone100 bezwzględne [GeV] / względne (ptcone/pT) |

D⁺: H (18), HDLZ (26 = H + DLZ8 z prerejestracji DLZ_dplus), HDLZia/HDLZir (zdefiniowane, nieużyte w blokach etapu 1).
„XGB H + nowe cechy”: wynik `dlz_data`/DLZ_d0 nie był gotowy przy składaniu paczki → `xgb.HX` (bez kinematyki) i `xgb.HK` (z kinematyką — kontrola
hipotezy „drzewo z kinematyką rzeźbi”); po wyniku DLZ_d0 wystarczy dopisać zestaw do `SETS["D0"]` i blok do `BLOCKS` (DN_SKIP_EXT=1).

| blok (job) | konfiguracje | pytanie |
|---|---|---|
| xgbA | xgb.H.sb ×3 seedy, xgb.H.mc ×3 | A: odniesienie i jego szum |
| xgbB | xgb.HX.{mc,sb}, xgb.HK.{mc,sb}, xgb.HKP.mc | A: nowe cechy w drzewie; czy kinematyka rzeźbi bez dekorelacji |
| mlpB1, mlpB2 | m64.H.mc (λ 0), m64.H.mc.d5 ×3 seedy, m64.H.sb.d5, r128.H.mc.d5 | B: sama zmiana modelu |
| mlpC1 | r128.HK.mc: λ = 0, 1, 5, 25, 100 | C: skan dekorelacji (krzywa moc ↔ rzeźbienie) |
| mlpC2 | r128.HK.mc.d5, .d25: seedy 2, 3 | C: stabilność |
| mlpC3 | r128.HK.mc.m5 (MoDe[1]), m64.HK.mc.d5, r128.HX.mc.d5, r128.HKP.mc.d5 | metoda dekorelacji, pojemność, rozdzielenie „cechy ↔ kinematyka”, dE/dx |
| mlpC4 | r128.HK.sbmc.d5, r128.HK.sb.d5 (kontrola pułapki), r128.HKia.mc.d5, r128.HKir.mc.d5 | tło z danych, izolacja |
| dpA, dpC | D⁺: xgb.{H,HDLZ}.{mc,sb}, m64.H.mc.d5, r128.HDLZ.mc (λ 0; d5 ×3 seedy), r128.HDLZ.sbmc.d5 | D⁺ |
| (analiza) | E: zespół = średnia logitów drzewa (H.sb; HX.mc) i sieci nominalnej | czy modele się uzupełniają |

Etap 2 (`Q_lxplus/dnet2/dn2_train.py`, tylko D⁰): tokeny = córki (log pT, Δη, Δφ, ładunek względem K, dE/dx, npts, d0, z0sinθ, rola) + do 20
torów otoczenia o największym pT w ΔR < 1,0 bez córek (te same pola bez d0/z0 — drzewo `Tracks` ich nie ma); cechy globalne = 18 cech H.
`ds` = DeepSets (średnia + maksimum), `tr` = transformer 2 warstwy, d = 64, 4 głowy, odczyt z tokenu globalnego. Konfiguracje:
ds/tr × „dau” (same córki), ds/tr × e04, tr × e10, tr.e04 bez dekorelacji, tr.e04 z tłem sbmc; 3 seedy dla tr.dau i tr.e04.

## 6. Protokół i metryki (`DNET/analyze.py`)

- Te same wiersze dla wszystkich modeli: kandydaci pT 3–30 GeV, zakres mas gatunku, L_xy ≥ 0,25 mm; sygnał = fd_nd truth 0/1 (D⁰: hypo_swap ≠ 1;
  D⁺: weto D*); tło MC = mb_nd truth ≥ 8; SB danych = |m − m₀| > 100 MeV (D⁰, + weto refleksji) / 110 MeV (D⁺, + weto D*).
  Wagi treningu: tło przeważone do widma (pT, |η|) sygnału, klasy zrównoważone — identyczne dla drzew i sieci.
- Progi WP z sygnału MC OOF w oknie sygnału per bin pT. Z_mc (tło MC w oknie), Z_tr (ε_B z SB danych 428770; osobno 427929), G = Z/Z_H z
  bootstrapem Poissona po zdarzeniach (200 replik, wspólne dla modeli).
- R_sculpt, nachylenie ε_B(m) (definicja FEAT `sculpt.py`), corr i dCor²(score, m | tło MC).
- SF z tagu D* (fit Δm pass/fail, profil) przy WP80/70/50; C2ST_1D sygnału = 0,5 + 0,5·max_WP |ε_dane − ε_MC| (dokładność najlepszego
  klasyfikatora progowego na samym score; pełne C2ST na wejściach jest w FEAT/DLZ).
- Rozrzut seedów, liczba parametrów, czas treningu z `meta/<tag>.json`. Modele etapu 2 oceniane na wierszach z wyciągniętym otoczeniem,
  odniesienie przeliczane na tych samych wierszach (kolumna „podzbiór”).

## 7. Paczki i obciążenie

| etap | paczka | DAG | węzły | zasoby | szacowany czas |
|---|---|---|---|---|---|
| 1 | `Q_lxplus/dnet1` (`submit_dn.sh`) | 24 ext → merge → 10 train → ana | **36** | ext 1 cpu / 6–10 GB (longlunch); train 4 cpu / 16 GB (workday); ana 8 cpu | ext ~10 min, merge ~10 min, train 0,3–2,5 h (DisCo ≈ 25 min na konfigurację: 5 foldów + model pełny), ana ~0,5–1 h → **~3–4 h + kolejka** |
| 2 | `Q_lxplus/dnet2` (`submit_dn2.sh`) | 24 ext → 9 train → ana | **34** | ext 1 cpu / 12 GB; train `request_gpus = 1` + 4 cpu / 24 GB (albo CPU: 8 cpu, tomorrow) | GPU: ~0,3 h/blok; CPU: ~1–3 h/blok (na 4 rdzeniach przy braku CUDA w slocie: do ~5 h) |

GPU: pula jest (`condor_status`: kilkanaście maszyn H100/H200/V100, część slotów wolna 17.11 rano), widok `LCG_110a_cuda` istnieje. Etap 1 —
**bez GPU**: MLP 56 tys. wag na 1,2·10⁶ wierszy liczy się ~3–6 min na model na 4 rdzeniach, czekanie na slot GPU byłoby dłuższe niż zysk. Etap 2 —
GPU opłaca się umiarkowanie (transformer, 14–22 tokeny: ~15 min/model na 8 rdzeniach, 5 foldów × 11 konfiguracji ≈ 14 h rdzenio-ósemek wobec
~1 h GPU); `DN2_GPU=1` jest domyślne, job sam spada na CPU, gdy w slocie nie ma CUDA; `DN2_GPU=0` → 8 rdzeni, flavour tomorrow.
Środowisko GPU (widok LCG CUDA) nie było testowane na węźle z GPU — tylko ścieżka CPU.
Wyjścia: EOS `charm_q/dnet` (parts ~1,5 GB, data ~4 GB, sc ~0,5 GB) i `charm_q/dnet2`; lokalnie pobierane tylko JSON + ONNX (< 20 MB).

## 8. Test mechaniki (17.11, lxplus `/tmp/klyczak_dnet`, 1 plik h5 na próbkę; bez condora, bez EOS/AFS)

Zakres: sprawdzenie, że kod działa od h5/ROOT do raportu. **Liczby z testu nie są wynikiem** (1 plik na próbkę: ~10⁴ sygnału, 2 foldy, 3–10 epok,
XGB w trybie szybkim 40 drzew; fity Δm w większości nieokreślone). Pliki: `DNET/test_2026-11-17/`.

| krok | wynik |
|---|---|
| `dn_job.sh ext` (fd, mb, 428770, 427929; D⁰ + D⁺ + pary D⁰←D*) | działa; 15–52 s na plik; na plik: fd D⁰ 16 168 / D⁺ 7 586 sygnału, mb D⁰ 39 300 / D⁺ 138 575, 428770 (15 % zdarzeń) D⁰ 27 437 + 15 471 par |
| `dn_job.sh merge` | działa |
| `dn_train.py train` — XGB (H sb/mc ×3 seedy; D⁺ H/HDLZ), MLP m64 i r128, DisCo λ 0/5/25/100/400, MoDe, tła mc/sb/sbmc, D⁰ i D⁺ | działa (25 konfiguracji); r128 = 55 553 wag (37 wejść), m64 = 9 601 (18 wejść) |
| eksport ONNX (MLP; eksport na 32 wierszach, kontrola onnxruntime na 256) | zgodność z torch 6·10⁻⁸; pierwsza wersja z `asinh` nie eksportowała się → zapis przez log/sqrt |
| `analyze.py compute` (bootstrap, R_sculpt, zespół, 216 fitów Δm pass/fail) + `report` | działa; ~4 s CPU na fit bez profilu → ~650 fitów z profilem na 8 rdzeniach ≈ 20–40 min |
| `make_dag_dn.sh` / `make_dag_dn2.sh` (bez wysyłania) | 36 i 34 węzły; listy 97 + 458 + 78 + 24 plików |
| `dn2_job.sh ext` (fd, mb; drzewo `Tracks`) | działa, 23–35 s na plik; własne Σ pT w ΔR < 0,4 koreluje z gałęzią ptcone40 na 0,996–0,9995 (mediana różnicy −0,16…−0,24 GeV: inna definicja progu/odjęcia córek) |
| `dn2_train.py` — ds/tr × dau/e04 | działa; ds 39 169 wag, tr 89 857; ONNX transformera (3 wejścia) zgodny na 6·10⁻⁸ przy innym rozmiarze wsadu |

Co test powiedział o ustawieniach (zmienione przed oddaniem): DisCo przy λ = 400 sprowadza sieć do stałej (G ≈ 0,3), przy λ = 100 traci ~35 %,
a sieć m64 zapadała się już przy λ = 25 w krótkim treningu → skan przesunięty na 1 / 5 / 25 / 100, nominał 5 i 25; podpróbka dekorelacji
2048 → 1024 (koszt DisCo ~0,1 s na wsad na 2 rdzeniach dominował czas). Nie testowano: uruchomienia pod condorem, zapisu na EOS (xrdcp),
środowiska GPU, pełnej pamięci (merge ~4 GB na próbkę 428770 — węzeł ma 16 GB).

## 9. Zastrzeżenia

- λ DisCo dobrane bez pełnej statystyki; jeśli skan pokaże, że optimum leży poza 5–25, seedy trzeba dołożyć dla innej wartości (`DN_LAM`, `DN_SKIP_EXT=1`).
  MoDe ma inną skalę kary (~10⁻³) — jedna wartość λ nie rozstrzyga o metodzie, tylko czy warto ją skanować.
- Sygnał D⁰: wiersze z poprawną hipotezą K/π (hypo_swap ≠ 1); refleksja nie jest ani sygnałem, ani tłem w treningu (jak FEAT).
- ε_S(MC) bez wag emulacji trygera; SF w J12only zawiera więc także niedopasowanie składu trygerowego (jak w CN_clone).
- Podpróbki zdarzeń: dane inkluzywne D⁰ 10 %, D⁺ 3 %, mb dla D⁺ 30 %; pary D⁰←D* 100 %.
- Etap 2: drzewo `Tracks` nie ma d0/z0 → otoczenie niesie tylko kinematykę, ładunek i dE/dx; wolny pion z D* nie jest usuwany z otoczenia.
- ONNX: eksportowany jest model trenowany na całości (jeden), nie układ foldów CharmNet (§2.5).
- XGB „H + nowe cechy” użyje listy z DLZ_d0 dopiero po zakończeniu etapu `dlz_data`.

## 10. Wyniki (po jobach; generowane)

<!-- WYNIKI:BEGIN -->
Wygenerowano 2026-11-18T09:09:19 przez `DNET/analyze.py report` z `charm_work/DNET/runs/lxplus/dnet_results.json,charm_work/DNET/runs/lxplus2/dnet2_results.json`.

### D0 — odniesienie `D0.xgb.H.sb.n.s1`; wierszy: sig 195639, bmc 1355464, bmc_all 3058753, sb428770 281918, dall428770 1442784, sb427929 97611, dall427929 501903

#### Zproxy względem XGB H (G = Z/Z_H, bootstrap sparowany), WP70; R_sculpt i corr na tle MC

| konfiguracja | cech | G_mc 3–5 | G_mc 5–8 | G_mc 8–30 | G_tr(SB 428770) 3–5 | 5–8 | 8–30 | R_sculpt MC (WP70) | R dane (pasmo) | nachylenie /100 MeV | corr(s, m) | par. | czas [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `D0.xgb.H.mc.n.s1` | 18 | 1.037 ± 0.002 | 1.059 ± 0.005 | 1.077 ± 0.018 | 0.971 ± 0.004 | 0.984 ± 0.009 | 0.929 ± 0.015 | 0.945 ± 0.006 | 0.971 ± 0.016 | -0.28 | -0.071 | 995 | 493.4 |
| `D0.xgb.H.mc.n.s2` | 18 | 1.033 ± 0.002 | 1.056 ± 0.005 | 1.063 ± 0.018 | 0.969 ± 0.004 | 0.994 ± 0.010 | 0.923 ± 0.015 | 0.952 ± 0.006 | 0.970 ± 0.016 | -0.27 | -0.071 | 952 | 470.9 |
| `D0.xgb.H.mc.n.s3` | 18 | 1.034 ± 0.002 | 1.055 ± 0.005 | 1.078 ± 0.018 | 0.967 ± 0.004 | 0.984 ± 0.009 | 0.918 ± 0.015 | 0.948 ± 0.006 | 0.966 ± 0.016 | -0.28 | -0.071 | 948 | 473.2 |
| `D0.xgb.H.sb.n.s1` | 18 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 0.942 ± 0.006 | 0.992 ± 0.017 | -0.30 | -0.085 | 1346 | 374.2 |
| `D0.xgb.H.sb.n.s2` | 18 | 0.996 ± 0.001 | 0.998 ± 0.004 | 0.976 ± 0.009 | 0.998 ± 0.003 | 0.996 ± 0.007 | 0.998 ± 0.010 | 0.946 ± 0.006 | 0.990 ± 0.017 | -0.30 | -0.085 | 1500 | 368.8 |
| `D0.xgb.H.sb.n.s3` | 18 | 0.998 ± 0.001 | 1.001 ± 0.003 | 0.987 ± 0.009 | 0.999 ± 0.003 | 1.005 ± 0.007 | 1.011 ± 0.011 | 0.945 ± 0.006 | 0.993 ± 0.017 | -0.30 | -0.085 | 1500 | 374.9 |
| `D0.xgb.HK.mc.n.s1` | 37 | 1.151 ± 0.003 | 1.156 ± 0.008 | 1.105 ± 0.019 | 1.108 ± 0.007 | 1.082 ± 0.015 | 0.957 ± 0.017 | 1.071 ± 0.008 | 1.017 ± 0.019 | -0.09 | -0.039 | 1448 | 565.0 |
| `D0.xgb.HK.sb.n.s1` | 37 | 1.038 ± 0.002 | 1.008 ± 0.007 | 0.971 ± 0.014 | 1.533 ± 0.014 | 1.246 ± 0.018 | 1.085 ± 0.016 | 1.608 ± 0.013 | 1.214 ± 0.028 | -0.42 | -0.153 | 1500 | 442.2 |
| `D0.xgb.HKP.mc.n.s1` | 44 | 1.154 ± 0.003 | 1.156 ± 0.008 | 1.138 ± 0.020 | 1.109 ± 0.007 | 1.081 ± 0.015 | 0.959 ± 0.017 | 1.071 ± 0.008 | 1.015 ± 0.019 | -0.09 | -0.038 | 924 | 708.5 |
| `D0.xgb.HX.mc.n.s1` | 29 | 1.136 ± 0.003 | 1.139 ± 0.008 | 1.137 ± 0.021 | 1.044 ± 0.005 | 1.046 ± 0.014 | 0.990 ± 0.017 | 0.969 ± 0.007 | 0.962 ± 0.017 | -0.17 | -0.048 | 1161 | 425.6 |
| `D0.xgb.HX.sb.n.s1` | 29 | 1.101 ± 0.002 | 1.074 ± 0.006 | 1.041 ± 0.013 | 1.081 ± 0.006 | 1.058 ± 0.012 | 1.055 ± 0.015 | 0.959 ± 0.007 | 0.963 ± 0.017 | -0.22 | -0.066 | 1270 | 295.9 |
| `D0.ds.H+dau.mc.d5.s1` (podzbiór) | 18 + 2×11 | 1.002 ± 0.003 | 1.015 ± 0.011 | 1.036 ± 0.026 | 0.881 ± 0.007 | 0.950 ± 0.015 | 0.997 ± 0.022 | 0.970 ± 0.010 | 0.966 ± 0.021 | -0.02 | 0.026 | 39169 | 245.5 |
| `D0.ds.H+e04.mc.d5.s1` (podzbiór) | 18 + 14×11 | 1.197 ± 0.005 | 1.259 ± 0.019 | 1.263 ± 0.047 | 1.285 ± 0.016 | 1.470 ± 0.042 | 1.412 ± 0.050 | 0.974 ± 0.012 | 0.984 ± 0.031 | -0.03 | 0.017 | 39169 | 258.9 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | – | 1.139 ± 0.002 | 1.133 ± 0.006 | 1.118 ± 0.014 | 1.081 ± 0.005 | 1.091 ± 0.012 | 1.058 ± 0.014 | 0.976 ± 0.007 | 0.984 ± 0.018 | -0.15 | -0.033 | – | 0 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | – | 1.172 ± 0.003 | 1.179 ± 0.008 | 1.184 ± 0.021 | 1.061 ± 0.006 | 1.071 ± 0.014 | 1.035 ± 0.017 | 0.984 ± 0.007 | 0.973 ± 0.017 | -0.07 | -0.013 | – | 0 |
| `D0.m64.H.mc.d5.s1` | 18 | 0.975 ± 0.002 | 1.010 ± 0.006 | 1.018 ± 0.015 | 0.894 ± 0.005 | 0.943 ± 0.011 | 0.957 ± 0.015 | 0.950 ± 0.006 | 0.981 ± 0.015 | -0.20 | -0.016 | 9601 | 1588.8 |
| `D0.m64.H.mc.d5.s2` | 18 | 0.980 ± 0.002 | 1.011 ± 0.006 | 1.047 ± 0.016 | 0.896 ± 0.005 | 0.945 ± 0.011 | 0.931 ± 0.015 | 0.948 ± 0.006 | 0.985 ± 0.015 | -0.20 | -0.015 | 9601 | 1526.3 |
| `D0.m64.H.mc.d5.s3` | 18 | 0.984 ± 0.002 | 1.016 ± 0.006 | 1.052 ± 0.016 | 0.904 ± 0.005 | 0.957 ± 0.011 | 0.957 ± 0.015 | 0.945 ± 0.006 | 0.979 ± 0.015 | -0.20 | -0.017 | 9601 | 1540.0 |
| `D0.m64.H.mc.n.s1` | 18 | 1.064 ± 0.002 | 1.083 ± 0.007 | 1.095 ± 0.019 | 0.997 ± 0.005 | 1.017 ± 0.013 | 0.985 ± 0.016 | 0.931 ± 0.006 | 0.973 ± 0.016 | -0.32 | -0.075 | 9601 | 189.8 |
| `D0.m64.H.sb.d5.s1` | 18 | 0.971 ± 0.002 | 0.966 ± 0.005 | 0.993 ± 0.014 | 0.935 ± 0.004 | 0.953 ± 0.010 | 0.954 ± 0.013 | 0.961 ± 0.006 | 0.985 ± 0.016 | -0.21 | -0.020 | 9601 | 684.4 |
| `D0.m64.HK.mc.d5.s1` | 37 | 1.133 ± 0.003 | 1.132 ± 0.008 | 1.127 ± 0.019 | 1.000 ± 0.006 | 1.009 ± 0.013 | 1.004 ± 0.018 | 0.988 ± 0.007 | 0.987 ± 0.017 | 0.06 | 0.022 | 10817 | 2213.8 |
| `D0.r128.H.mc.d5.s1` | 18 | 1.037 ± 0.002 | 1.078 ± 0.006 | 1.069 ± 0.017 | 0.960 ± 0.005 | 0.994 ± 0.011 | 0.988 ± 0.017 | 0.953 ± 0.006 | 0.973 ± 0.016 | -0.22 | -0.008 | 53121 | 2077.8 |
| `D0.r128.HK.mc.d1.s1` | 37 | 1.192 ± 0.003 | 1.182 ± 0.009 | 1.135 ± 0.023 | 1.880 ± 0.024 | 1.390 ± 0.026 | 1.107 ± 0.023 | 2.604 ± 0.027 | 1.466 ± 0.040 | -0.02 | 0.015 | 55553 | 2734.1 |
| `D0.r128.HK.mc.d100.s1` | 37 | 0.562 ± 0.002 | 0.302 ± 0.002 | 0.303 ± 0.007 | 0.510 ± 0.004 | 0.238 ± 0.004 | 0.256 ± 0.006 | 0.996 ± 0.003 | 0.990 ± 0.006 | 0.04 | 0.006 | 55553 | 2808.1 |
| `D0.r128.HK.mc.d25.s1` | 37 | 0.984 ± 0.002 | 1.014 ± 0.007 | 0.998 ± 0.016 | 0.875 ± 0.006 | 0.924 ± 0.011 | 0.927 ± 0.016 | 0.985 ± 0.006 | 0.985 ± 0.015 | 0.00 | 0.018 | 55553 | 2680.5 |
| `D0.r128.HK.mc.d25.s2` | 37 | 0.994 ± 0.002 | 1.017 ± 0.007 | 0.991 ± 0.016 | 0.875 ± 0.006 | 0.933 ± 0.012 | 0.932 ± 0.015 | 0.982 ± 0.006 | 0.976 ± 0.015 | 0.02 | 0.020 | 55553 | 2072.3 |
| `D0.r128.HK.mc.d25.s3` | 37 | 1.014 ± 0.003 | 1.038 ± 0.007 | 1.007 ± 0.016 | 0.894 ± 0.006 | 0.957 ± 0.012 | 0.915 ± 0.015 | 0.983 ± 0.006 | 0.971 ± 0.015 | 0.03 | 0.021 | 55553 | 1826.3 |
| `D0.r128.HK.mc.d5.s1` | 37 | 1.161 ± 0.003 | 1.163 ± 0.009 | 1.152 ± 0.019 | 1.029 ± 0.007 | 1.047 ± 0.014 | 1.037 ± 0.018 | 0.992 ± 0.007 | 0.979 ± 0.017 | 0.03 | 0.021 | 55553 | 2850.3 |
| `D0.r128.HK.mc.d5.s2` | 37 | 1.158 ± 0.003 | 1.167 ± 0.009 | 1.137 ± 0.019 | 1.022 ± 0.006 | 1.031 ± 0.014 | 1.009 ± 0.017 | 0.994 ± 0.007 | 0.976 ± 0.017 | 0.04 | 0.021 | 55553 | 2041.0 |
| `D0.r128.HK.mc.d5.s3` | 37 | 1.158 ± 0.003 | 1.170 ± 0.009 | 1.138 ± 0.019 | 1.026 ± 0.006 | 1.052 ± 0.014 | 1.041 ± 0.017 | 0.995 ± 0.007 | 0.980 ± 0.017 | 0.04 | 0.022 | 55553 | 2131.5 |
| `D0.r128.HK.mc.m5.s1` | 37 | 1.210 ± 0.004 | 1.188 ± 0.010 | 1.124 ± 0.022 | 2.966 ± 0.056 | 1.884 ± 0.053 | 1.207 ± 0.028 | 5.980 ± 0.089 | 1.698 ± 0.064 | -0.06 | -0.057 | 55553 | 1357.3 |
| `D0.r128.HK.mc.n.s1` | 37 | 1.254 ± 0.004 | 1.216 ± 0.010 | 1.123 ± 0.021 | 4.116 ± 0.118 | 2.443 ± 0.098 | 1.420 ± 0.040 | 10.896 ± 0.222 | 1.864 ± 0.090 | -0.09 | -0.037 | 55553 | 1028.2 |
| `D0.r128.HK.sb.d5.s1` | 37 | 1.114 ± 0.003 | 1.085 ± 0.007 | 1.078 ± 0.016 | 1.258 ± 0.008 | 1.149 ± 0.016 | 1.093 ± 0.019 | 1.253 ± 0.010 | 1.089 ± 0.022 | -0.17 | -0.003 | 55553 | 966.3 |
| `D0.r128.HK.sbmc.d5.s1` | 37 | 1.040 ± 0.003 | 1.062 ± 0.007 | 1.042 ± 0.015 | 1.766 ± 0.019 | 1.386 ± 0.025 | 1.118 ± 0.019 | 2.192 ± 0.019 | 1.458 ± 0.038 | -0.30 | 0.003 | 55553 | 2216.8 |
| `D0.r128.HKP.mc.d5.s1` | 44 | 1.161 ± 0.003 | 1.165 ± 0.009 | 1.153 ± 0.019 | 1.017 ± 0.006 | 1.046 ± 0.013 | 1.024 ± 0.018 | 0.988 ± 0.007 | 0.965 ± 0.017 | 0.03 | 0.021 | 56449 | 1910.2 |
| `D0.r128.HKia.mc.d5.s1` | 39 | 1.342 ± 0.004 | 1.345 ± 0.012 | 1.340 ± 0.028 | 1.609 ± 0.016 | 1.593 ± 0.035 | 1.398 ± 0.036 | 0.995 ± 0.009 | 0.983 ± 0.027 | 0.03 | 0.018 | 55809 | 1989.1 |
| `D0.r128.HKir.mc.d5.s1` | 39 | 1.334 ± 0.004 | 1.344 ± 0.012 | 1.342 ± 0.030 | 1.589 ± 0.016 | 1.612 ± 0.035 | 1.441 ± 0.038 | 0.999 ± 0.009 | 0.986 ± 0.027 | 0.03 | 0.018 | 55809 | 1993.0 |
| `D0.r128.HX.mc.d5.s1` | 29 | 1.148 ± 0.003 | 1.165 ± 0.008 | 1.161 ± 0.021 | 1.036 ± 0.007 | 1.058 ± 0.015 | 1.020 ± 0.018 | 0.980 ± 0.007 | 0.977 ± 0.017 | -0.09 | 0.001 | 54529 | 3671.1 |
| `D0.tr.H+dau.mc.d5.s1` (podzbiór) | 18 + 2×11 | 1.048 ± 0.004 | 1.069 ± 0.011 | 1.036 ± 0.025 | 0.939 ± 0.008 | 0.978 ± 0.017 | 0.969 ± 0.021 | 0.969 ± 0.010 | 0.958 ± 0.022 | -0.09 | 0.016 | 89857 | 355.9 |
| `D0.tr.H+dau.mc.d5.s2` (podzbiór) | 18 + 2×11 | 1.034 ± 0.003 | 1.060 ± 0.012 | 1.050 ± 0.029 | 0.937 ± 0.007 | 0.977 ± 0.017 | 0.982 ± 0.022 | 0.972 ± 0.010 | 0.964 ± 0.022 | -0.10 | 0.017 | 89857 | 303.4 |
| `D0.tr.H+dau.mc.d5.s3` (podzbiór) | 18 + 2×11 | 1.046 ± 0.004 | 1.053 ± 0.011 | 1.065 ± 0.029 | 0.940 ± 0.008 | 0.979 ± 0.016 | 0.977 ± 0.023 | 0.965 ± 0.010 | 0.972 ± 0.022 | -0.10 | 0.016 | 89857 | 291.0 |
| `D0.tr.H+e04.mc.d5.s1` (podzbiór) | 18 + 14×11 | 1.473 ± 0.008 | 1.652 ± 0.035 | 1.551 ± 0.068 | 1.580 ± 0.024 | 1.789 ± 0.062 | 1.600 ± 0.066 | 0.975 ± 0.015 | 1.014 ± 0.039 | -0.10 | 0.005 | 89857 | 19127.3 |
| `D0.tr.H+e04.mc.d5.s2` (podzbiór) | 18 + 14×11 | 1.452 ± 0.008 | 1.601 ± 0.029 | 1.538 ± 0.064 | 1.547 ± 0.023 | 1.774 ± 0.062 | 1.628 ± 0.072 | 0.979 ± 0.015 | 0.998 ± 0.038 | -0.10 | 0.003 | 89857 | 428.3 |
| `D0.tr.H+e04.mc.d5.s3` (podzbiór) | 18 + 14×11 | 1.425 ± 0.007 | 1.592 ± 0.031 | 1.551 ± 0.070 | 1.540 ± 0.022 | 1.721 ± 0.057 | 1.607 ± 0.073 | 0.973 ± 0.014 | 1.007 ± 0.038 | -0.13 | 0.003 | 89857 | 15158.5 |
| `D0.tr.H+e04.mc.n.s1` (podzbiór) | 18 + 14×11 | 1.579 ± 0.010 | 1.711 ± 0.035 | 1.585 ± 0.073 | 1.692 ± 0.029 | 1.789 ± 0.060 | 1.554 ± 0.068 | 0.969 ± 0.016 | 0.973 ± 0.039 | -0.25 | -0.053 | 89857 | 16186.6 |
| `D0.tr.H+e04.sbmc.d5.s1` (podzbiór) | 18 + 14×11 | 1.374 ± 0.007 | 1.450 ± 0.025 | 1.325 ± 0.050 | 1.653 ± 0.027 | 1.848 ± 0.062 | 1.768 ± 0.089 | 0.971 ± 0.014 | 1.050 ± 0.043 | -0.18 | -0.016 | 89857 | 12700.4 |
| `D0.tr.H+e10.mc.d5.s1` (podzbiór) | 18 + 22×11 | 1.620 ± 0.010 | 1.712 ± 0.034 | 1.662 ± 0.078 | 1.940 ± 0.033 | 1.889 ± 0.073 | 1.750 ± 0.087 | 0.971 ± 0.016 | 1.001 ± 0.046 | -0.08 | 0.001 | 89857 | 25007.3 |

#### WP80 / WP50 (G_mc; G_tr)

| konfiguracja | WP80 3–5 | 5–8 | 8–30 | WP50 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|---|---|
| `D0.ds.H+dau.mc.d5.s1` | 1.00; 0.89 | 1.01; 0.94 | 1.04; 0.98 | 1.01; 0.90 | 1.04; 0.92 | 1.02; 0.94 |
| `D0.ds.H+e04.mc.d5.s1` | 1.18; 1.29 | 1.24; 1.49 | 1.27; 1.43 | 1.20; 1.24 | 1.24; 1.25 | 1.07; 1.03 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 1.11; 1.07 | 1.10; 1.06 | 1.11; 1.06 | 1.17; 1.10 | 1.14; 1.07 | 1.12; 1.11 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 1.14; 1.05 | 1.14; 1.05 | 1.16; 1.02 | 1.21; 1.07 | 1.17; 1.03 | 1.21; 0.97 |
| `D0.m64.H.mc.d5.s1` | 0.96; 0.88 | 0.98; 0.91 | 1.02; 0.94 | 0.99; 0.93 | 1.04; 1.00 | 1.07; 0.96 |
| `D0.m64.H.mc.d5.s2` | 0.97; 0.88 | 0.99; 0.90 | 1.02; 0.94 | 1.00; 0.93 | 1.04; 0.96 | 1.05; 0.94 |
| `D0.m64.H.mc.d5.s3` | 0.97; 0.89 | 0.99; 0.91 | 1.03; 0.95 | 1.01; 0.95 | 1.05; 0.97 | 1.10; 0.95 |
| `D0.m64.H.mc.n.s1` | 1.04; 0.99 | 1.07; 1.00 | 1.07; 0.99 | 1.10; 1.01 | 1.11; 1.00 | 1.11; 0.95 |
| `D0.m64.H.sb.d5.s1` | 0.96; 0.92 | 0.98; 0.94 | 0.99; 0.96 | 0.96; 0.94 | 0.99; 0.97 | 0.97; 0.98 |
| `D0.m64.HK.mc.d5.s1` | 1.10; 0.98 | 1.10; 1.00 | 1.11; 0.98 | 1.15; 0.98 | 1.14; 1.00 | 1.10; 0.94 |
| `D0.r128.H.mc.d5.s1` | 1.01; 0.92 | 1.03; 0.94 | 1.07; 0.97 | 1.09; 1.01 | 1.10; 1.01 | 1.16; 0.97 |
| `D0.r128.HK.mc.d1.s1` | 1.15; 1.72 | 1.14; 1.39 | 1.14; 1.05 | 1.22; 1.79 | 1.16; 1.23 | 1.15; 1.02 |
| `D0.r128.HK.mc.d100.s1` | 0.53; 0.46 | 0.38; 0.30 | 0.36; 0.29 | 0.52; 0.48 | 0.73; 0.72 | 0.21; 0.18 |
| `D0.r128.HK.mc.d25.s1` | 0.96; 0.87 | 0.98; 0.91 | 0.98; 0.92 | 0.99; 0.88 | 1.01; 0.91 | 0.99; 0.84 |
| `D0.r128.HK.mc.d25.s2` | 0.94; 0.85 | 0.96; 0.91 | 0.98; 0.93 | 1.04; 0.92 | 1.05; 0.95 | 1.03; 0.88 |
| `D0.r128.HK.mc.d25.s3` | 0.97; 0.88 | 0.99; 0.91 | 0.99; 0.92 | 1.05; 0.92 | 1.04; 0.93 | 1.01; 0.86 |
| `D0.r128.HK.mc.d5.s1` | 1.12; 1.02 | 1.13; 1.03 | 1.12; 1.01 | 1.19; 1.03 | 1.16; 1.02 | 1.15; 0.99 |
| `D0.r128.HK.mc.d5.s2` | 1.12; 1.01 | 1.13; 1.02 | 1.13; 1.02 | 1.19; 1.03 | 1.17; 1.03 | 1.16; 0.95 |
| `D0.r128.HK.mc.d5.s3` | 1.12; 1.01 | 1.12; 1.02 | 1.12; 1.01 | 1.19; 1.04 | 1.17; 1.03 | 1.19; 0.97 |
| `D0.r128.HK.mc.m5.s1` | 1.19; 2.72 | 1.17; 1.96 | 1.14; 1.19 | 1.23; 2.77 | 1.14; 1.57 | 1.17; 1.08 |
| `D0.r128.HK.mc.n.s1` | 1.21; 4.20 | 1.21; 2.79 | 1.15; 1.48 | 1.28; 3.61 | 1.17; 1.97 | 1.15; 1.25 |
| `D0.r128.HK.sb.d5.s1` | 1.10; 1.19 | 1.09; 1.13 | 1.09; 1.06 | 1.11; 1.26 | 1.07; 1.13 | 1.03; 1.04 |
| `D0.r128.HK.sbmc.d5.s1` | 1.03; 1.57 | 1.06; 1.31 | 1.06; 1.07 | 1.03; 1.74 | 1.03; 1.32 | 0.97; 1.05 |
| `D0.r128.HKP.mc.d5.s1` | 1.12; 1.01 | 1.13; 1.02 | 1.12; 1.00 | 1.20; 1.04 | 1.17; 1.02 | 1.16; 0.97 |
| `D0.r128.HKia.mc.d5.s1` | 1.28; 1.53 | 1.31; 1.65 | 1.32; 1.43 | 1.37; 1.64 | 1.32; 1.44 | 1.48; 1.28 |
| `D0.r128.HKir.mc.d5.s1` | 1.28; 1.52 | 1.30; 1.66 | 1.31; 1.45 | 1.36; 1.65 | 1.31; 1.45 | 1.37; 1.25 |
| `D0.r128.HX.mc.d5.s1` | 1.11; 1.01 | 1.12; 1.01 | 1.12; 0.99 | 1.18; 1.05 | 1.15; 1.03 | 1.18; 1.01 |
| `D0.tr.H+dau.mc.d5.s1` | 1.02; 0.92 | 1.04; 0.96 | 1.04; 0.96 | 1.10; 0.99 | 1.09; 0.99 | 1.13; 0.99 |
| `D0.tr.H+dau.mc.d5.s2` | 1.02; 0.92 | 1.04; 0.95 | 1.04; 0.96 | 1.08; 0.97 | 1.09; 1.00 | 1.19; 0.98 |
| `D0.tr.H+dau.mc.d5.s3` | 1.03; 0.93 | 1.04; 0.95 | 1.03; 0.97 | 1.09; 1.00 | 1.06; 0.96 | 1.15; 0.99 |
| `D0.tr.H+e04.mc.d5.s1` | 1.39; 1.49 | 1.57; 1.79 | 1.55; 1.73 | 1.52; 1.66 | 1.70; 1.73 | 1.85; 1.51 |
| `D0.tr.H+e04.mc.d5.s2` | 1.35; 1.48 | 1.52; 1.73 | 1.62; 1.68 | 1.49; 1.57 | 1.74; 1.77 | 1.77; 1.47 |
| `D0.tr.H+e04.mc.d5.s3` | 1.34; 1.45 | 1.50; 1.69 | 1.49; 1.62 | 1.47; 1.59 | 1.73; 1.69 | 1.61; 1.49 |
| `D0.tr.H+e04.mc.n.s1` | 1.47; 1.61 | 1.64; 1.79 | 1.56; 1.53 | 1.61; 1.73 | 1.79; 1.66 | 1.61; 1.43 |
| `D0.tr.H+e04.sbmc.d5.s1` | 1.30; 1.60 | 1.41; 1.86 | 1.38; 1.75 | 1.37; 1.62 | 1.42; 1.94 | 1.26; 1.67 |
| `D0.tr.H+e10.mc.d5.s1` | 1.51; 1.80 | 1.63; 1.96 | 1.60; 1.72 | 1.66; 2.06 | 1.75; 1.86 | 1.94; 1.64 |
| `D0.xgb.H.mc.n.s1` | 1.02; 0.97 | 1.04; 0.98 | 1.05; 0.93 | 1.05; 0.97 | 1.06; 0.98 | 1.01; 0.83 |
| `D0.xgb.H.mc.n.s2` | 1.02; 0.97 | 1.04; 0.98 | 1.06; 0.94 | 1.05; 0.97 | 1.06; 0.98 | 1.07; 0.83 |
| `D0.xgb.H.mc.n.s3` | 1.02; 0.97 | 1.05; 0.98 | 1.06; 0.93 | 1.05; 0.97 | 1.07; 0.97 | 1.02; 0.84 |
| `D0.xgb.H.sb.n.s1` | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 |
| `D0.xgb.H.sb.n.s2` | 1.00; 1.00 | 1.00; 1.00 | 0.99; 1.00 | 1.01; 1.00 | 0.99; 1.00 | 0.99; 0.99 |
| `D0.xgb.H.sb.n.s3` | 1.00; 1.00 | 1.00; 1.00 | 0.99; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 0.96; 0.97 |
| `D0.xgb.HK.mc.n.s1` | 1.12; 1.10 | 1.13; 1.07 | 1.13; 0.97 | 1.17; 1.13 | 1.13; 1.05 | 1.13; 0.84 |
| `D0.xgb.HK.sb.n.s1` | 1.02; 1.46 | 0.99; 1.25 | 0.95; 1.12 | 1.06; 1.56 | 1.06; 1.23 | 1.01; 1.03 |
| `D0.xgb.HKP.mc.n.s1` | 1.13; 1.10 | 1.13; 1.06 | 1.14; 0.99 | 1.18; 1.13 | 1.14; 1.05 | 1.12; 0.85 |
| `D0.xgb.HX.mc.n.s1` | 1.11; 1.04 | 1.12; 1.04 | 1.13; 0.98 | 1.16; 1.06 | 1.12; 1.02 | 1.13; 0.86 |
| `D0.xgb.HX.sb.n.s1` | 1.08; 1.08 | 1.07; 1.07 | 1.06; 1.04 | 1.13; 1.09 | 1.09; 1.07 | 1.00; 1.02 |

#### Rozrzut seedów (odchylenie std G_mc WP70 między seedami)

| rodzina | seedy | 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|
| `D0.m64.H.mc.d5` | 3 | 0.980 ± 0.004 | 1.012 ± 0.003 | 1.039 ± 0.018 |
| `D0.r128.HK.mc.d25` | 3 | 0.998 ± 0.015 | 1.023 ± 0.013 | 0.999 ± 0.008 |
| `D0.r128.HK.mc.d5` | 3 | 1.159 ± 0.002 | 1.167 ± 0.004 | 1.142 ± 0.008 |
| `D0.xgb.H.mc.n` | 3 | 1.034 ± 0.002 | 1.057 ± 0.002 | 1.072 ± 0.008 |
| `D0.xgb.H.sb.n` | 3 | 0.998 ± 0.002 | 1.000 ± 0.001 | 0.987 ± 0.012 |
| `D0.tr.H+dau.mc.d5` | 3 | 1.042 ± 0.007 | 1.061 ± 0.008 | 1.050 ± 0.015 |
| `D0.tr.H+e04.mc.d5` | 3 | 1.450 ± 0.024 | 1.615 ± 0.032 | 1.546 ± 0.008 |

#### SF = ε_S(dane)/ε_S(MC) z tagu D* (fit Δm pass/fail) i C2ST_1D sygnału

| konfiguracja | run:klasa | bin | SF WP80 | SF WP70 | SF WP50 | C2ST_1D |
|---|---|---|---|---|---|---|
| `D0.ds.H+dau.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.34 ± 0.40 (n) | 1.41 ± 0.51 (n) | 1.34 ± 0.47 | 0.586 |
| `D0.ds.H+dau.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.05 ± 0.12 | 1.04 ± 0.13 | 0.92 ± 0.13 | 0.520 |
| `D0.ds.H+dau.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.95 ± 0.08 | 0.95 ± 0.08 | 0.98 ± 0.12 | 0.522 |
| `D0.ds.H+dau.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.72 (n) | 2.13 ± 2.38 (n) | 2.28 ± 2.72 (n) | – |
| `D0.ds.H+dau.mc.d5.s1` | 428770:J12only | 5-8 | 0.88 ± 0.16 | 0.95 ± 0.18 | 1.04 ± 0.21 | 0.546 |
| `D0.ds.H+dau.mc.d5.s1` | 428770:J12only | 8-30 | 0.99 ± 0.05 | 1.00 ± 0.05 | 1.04 ± 0.06 | 0.509 |
| `D0.ds.H+dau.mc.d5.s1` | 428770:all | 3-5 | 1.83 ± 2.09 (n) | 1.78 ± 0.71 (n) | 1.81 ± 0.87 (n) | – |
| `D0.ds.H+dau.mc.d5.s1` | 428770:all | 5-8 | 0.83 ± 0.11 | 0.88 ± 0.12 | 0.94 ± 0.13 | 0.567 |
| `D0.ds.H+dau.mc.d5.s1` | 428770:all | 8-30 | 0.97 ± 0.04 | 0.97 ± 0.04 | 1.00 ± 0.05 | 0.511 |
| `D0.ds.H+e04.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.33 ± 0.38 (n) | 1.36 ± 0.42 | 1.29 ± 0.45 (n) | 0.640 |
| `D0.ds.H+e04.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.03 ± 0.12 | 1.04 ± 0.13 | 1.05 ± 0.14 | 0.516 |
| `D0.ds.H+e04.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.96 ± 0.09 | 0.95 ± 0.09 | 0.98 ± 0.12 | 0.520 |
| `D0.ds.H+e04.mc.d5.s1` | 428770:J12only | 3-5 | 1.62 ± 1.84 (n) | 1.53 ± 1.26 (n) | 2.05 ± 2.43 (n) | – |
| `D0.ds.H+e04.mc.d5.s1` | 428770:J12only | 5-8 | 0.89 ± 0.17 | 0.89 ± 0.17 | 0.94 ± 0.19 | 0.546 |
| `D0.ds.H+e04.mc.d5.s1` | 428770:J12only | 8-30 | 0.94 ± 0.05 | 0.94 ± 0.05 | 0.96 ± 0.05 | 0.527 |
| `D0.ds.H+e04.mc.d5.s1` | 428770:all | 3-5 | 1.61 ± 1.14 (n) | 1.52 ± 0.62 (n) | 1.69 ± 0.80 (n) | – |
| `D0.ds.H+e04.mc.d5.s1` | 428770:all | 5-8 | 0.82 ± 0.11 | 0.81 ± 0.11 | 0.84 ± 0.12 | 0.576 |
| `D0.ds.H+e04.mc.d5.s1` | 428770:all | 8-30 | 0.93 ± 0.04 | 0.93 ± 0.04 | 0.94 ± 0.05 | 0.529 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.35 ± 0.38 (n) | 1.34 ± 0.47 (n) | 1.39 ± 0.48 | 0.598 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.06 ± 0.13 | 0.95 ± 0.13 | 0.528 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.92 ± 0.08 | 0.93 ± 0.08 | 1.03 ± 0.12 | 0.531 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 2.05 (n) | 2.13 ± 2.41 (n) | 2.16 ± 2.57 (n) | – |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 428770:J12only | 5-8 | 0.92 ± 0.17 | 0.95 ± 0.18 | 1.00 ± 0.20 | 0.534 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 1.01 ± 0.05 | 0.98 ± 0.06 | 0.506 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 428770:all | 3-5 | 1.73 ± 1.73 (n) | 1.76 ± 0.72 (n) | 1.82 ± 0.85 (n) | – |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 428770:all | 5-8 | 0.85 ± 0.11 | 0.89 ± 0.12 | 0.91 ± 0.13 | 0.560 |
| `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 0.98 ± 0.04 | 0.94 ± 0.05 | 0.515 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.35 ± 0.38 (n) | 1.38 ± 0.46 | 1.40 ± 0.49 | 0.633 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.05 ± 0.13 | 1.06 ± 0.13 | 0.99 ± 0.14 | 0.522 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.08 | 0.93 ± 0.08 | 1.05 ± 0.13 | 0.528 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.64 (n) | 2.13 ± 2.39 (n) | 2.13 ± 2.55 (n) | – |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 428770:J12only | 5-8 | 0.94 ± 0.18 | 0.96 ± 0.18 | 0.97 ± 0.19 | 0.524 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 428770:J12only | 8-30 | 1.03 ± 0.05 | 1.05 ± 0.05 | 1.08 ± 0.06 | 0.520 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 428770:all | 3-5 | 1.76 ± 1.81 (n) | 1.80 ± 1.35 (n) | 1.79 ± 0.83 (n) | – |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 428770:all | 5-8 | 0.86 ± 0.12 | 0.89 ± 0.12 | 0.90 ± 0.13 | 0.554 |
| `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1` | 428770:all | 8-30 | 1.01 ± 0.04 | 1.01 ± 0.04 | 1.04 ± 0.05 | 0.510 |
| `D0.m64.H.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.52 ± 0.44 (n) | 1.38 ± 0.47 (n) | 1.38 ± 0.48 | 0.595 |
| `D0.m64.H.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.05 ± 0.12 | 1.03 ± 0.13 | 0.94 ± 0.13 | 0.518 |
| `D0.m64.H.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.94 ± 0.08 | 0.93 ± 0.09 | 0.93 ± 0.11 | 0.523 |
| `D0.m64.H.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.49 (n) | 2.13 ± 2.42 (n) | 2.30 ± 2.73 (n) | – |
| `D0.m64.H.mc.d5.s1` | 428770:J12only | 5-8 | 0.92 ± 0.17 | 0.95 ± 0.18 | 1.04 ± 0.21 | 0.532 |
| `D0.m64.H.mc.d5.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 1.00 ± 0.05 | 1.02 ± 0.06 | 0.504 |
| `D0.m64.H.mc.d5.s1` | 428770:all | 3-5 | 1.87 ± 2.11 (n) | 1.84 ± 1.82 (n) | 1.92 ± 0.92 (n) | – |
| `D0.m64.H.mc.d5.s1` | 428770:all | 5-8 | 0.87 ± 0.12 | 0.89 ± 0.12 | 0.94 ± 0.13 | 0.553 |
| `D0.m64.H.mc.d5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 0.98 ± 0.04 | 0.98 ± 0.05 | 0.509 |
| `D0.m64.H.mc.n.s1` | 427929:MB427929 | 3-5 | 1.46 ± 0.42 (n) | 1.29 ± 0.44 | 1.36 ± 0.47 | 0.601 |
| `D0.m64.H.mc.n.s1` | 427929:MB427929 | 5-8 | 1.03 ± 0.12 | 1.05 ± 0.13 | 0.96 ± 0.13 | 0.518 |
| `D0.m64.H.mc.n.s1` | 427929:MB427929 | 8-30 | 0.91 ± 0.07 | 0.94 ± 0.08 | 0.95 ± 0.10 | 0.538 |
| `D0.m64.H.mc.n.s1` | 428770:J12only | 3-5 | 1.87 ± 1.61 (n) | 2.13 ± 2.34 (n) | 2.38 ± 2.79 (n) | – |
| `D0.m64.H.mc.n.s1` | 428770:J12only | 5-8 | 0.93 ± 0.17 | 0.96 ± 0.19 | 1.02 ± 0.20 | 0.527 |
| `D0.m64.H.mc.n.s1` | 428770:J12only | 8-30 | 1.02 ± 0.05 | 1.02 ± 0.05 | 1.03 ± 0.06 | 0.509 |
| `D0.m64.H.mc.n.s1` | 428770:all | 3-5 | 1.87 ± 2.08 (n) | 1.90 ± 1.92 (n) | 1.90 ± 0.89 (n) | – |
| `D0.m64.H.mc.n.s1` | 428770:all | 5-8 | 0.88 ± 0.12 | 0.90 ± 0.12 | 0.93 ± 0.13 | 0.550 |
| `D0.m64.H.mc.n.s1` | 428770:all | 8-30 | 1.00 ± 0.04 | 1.00 ± 0.04 | 1.00 ± 0.05 | 0.501 |
| `D0.m64.H.sb.d5.s1` | 427929:MB427929 | 3-5 | 1.40 ± 0.38 (n) | 1.37 ± 0.45 | 1.34 ± 0.45 | 0.630 |
| `D0.m64.H.sb.d5.s1` | 427929:MB427929 | 5-8 | 1.05 ± 0.12 | 1.01 ± 0.13 | 0.93 ± 0.13 | 0.519 |
| `D0.m64.H.sb.d5.s1` | 427929:MB427929 | 8-30 | 0.97 ± 0.09 | 0.96 ± 0.09 | 0.93 ± 0.11 | 0.518 |
| `D0.m64.H.sb.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.66 (n) | 2.13 ± 2.27 (n) | 2.23 ± 2.67 (n) | – |
| `D0.m64.H.sb.d5.s1` | 428770:J12only | 5-8 | 0.95 ± 0.18 | 0.96 ± 0.18 | 0.99 ± 0.19 | 0.522 |
| `D0.m64.H.sb.d5.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 0.99 ± 0.05 | 0.96 ± 0.05 | 0.511 |
| `D0.m64.H.sb.d5.s1` | 428770:all | 3-5 | 1.87 ± 1.94 (n) | 1.86 ± 1.88 (n) | 1.83 ± 0.87 (n) | – |
| `D0.m64.H.sb.d5.s1` | 428770:all | 5-8 | 0.87 ± 0.12 | 0.91 ± 0.12 | 0.90 ± 0.13 | 0.551 |
| `D0.m64.H.sb.d5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 0.96 ± 0.04 | 0.94 ± 0.05 | 0.516 |
| `D0.m64.HK.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.34 ± 0.40 | 1.37 ± 0.47 (n) | 1.31 ± 0.45 | 0.637 |
| `D0.m64.HK.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.05 ± 0.12 | 1.09 ± 0.14 | 1.00 ± 0.14 | 0.531 |
| `D0.m64.HK.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.95 ± 0.09 | 0.96 ± 0.09 | 1.00 ± 0.12 | 0.522 |
| `D0.m64.HK.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.74 (n) | 2.14 ± 2.10 (n) | 2.23 ± 2.65 (n) | – |
| `D0.m64.HK.mc.d5.s1` | 428770:J12only | 5-8 | 0.93 ± 0.17 | 0.95 ± 0.18 | 1.02 ± 0.20 | 0.527 |
| `D0.m64.HK.mc.d5.s1` | 428770:J12only | 8-30 | 1.01 ± 0.05 | 1.01 ± 0.05 | 1.04 ± 0.06 | 0.511 |
| `D0.m64.HK.mc.d5.s1` | 428770:all | 3-5 | 1.76 ± 1.36 (n) | 1.90 ± 1.91 (n) | 1.83 ± 0.86 (n) | – |
| `D0.m64.HK.mc.d5.s1` | 428770:all | 5-8 | 0.86 ± 0.11 | 0.90 ± 0.12 | 0.92 ± 0.13 | 0.558 |
| `D0.m64.HK.mc.d5.s1` | 428770:all | 8-30 | 0.99 ± 0.04 | 0.98 ± 0.04 | 1.01 ± 0.05 | 0.506 |
| `D0.r128.H.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.36 ± 0.40 (n) | 1.38 ± 0.45 | 1.36 ± 0.47 | 0.632 |
| `D0.r128.H.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.01 ± 0.13 | 1.00 ± 0.14 | 0.528 |
| `D0.r128.H.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.92 ± 0.08 | 0.95 ± 0.09 | 0.94 ± 0.11 | 0.531 |
| `D0.r128.H.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.47 (n) | 2.13 ± 2.12 (n) | 2.11 ± 2.52 (n) | – |
| `D0.r128.H.mc.d5.s1` | 428770:J12only | 5-8 | 0.95 ± 0.18 | 0.95 ± 0.18 | 0.96 ± 0.19 | 0.519 |
| `D0.r128.H.mc.d5.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 1.00 ± 0.05 | 0.99 ± 0.06 | 0.502 |
| `D0.r128.H.mc.d5.s1` | 428770:all | 3-5 | 1.87 ± 1.81 (n) | 1.98 ± 2.07 (n) | 1.80 ± 0.84 (n) | – |
| `D0.r128.H.mc.d5.s1` | 428770:all | 5-8 | 0.88 ± 0.12 | 0.88 ± 0.12 | 0.89 ± 0.13 | 0.547 |
| `D0.r128.H.mc.d5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 0.98 ± 0.04 | 0.96 ± 0.05 | 0.510 |
| `D0.r128.HK.mc.d1.s1` | 427929:MB427929 | 3-5 | 1.33 ± 0.37 (n) | 1.32 ± 0.45 (n) | 1.31 ± 0.45 | 0.579 |
| `D0.r128.HK.mc.d1.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.07 ± 0.13 | 1.05 ± 0.14 | 0.528 |
| `D0.r128.HK.mc.d1.s1` | 427929:MB427929 | 8-30 | 0.92 ± 0.08 | 0.90 ± 0.07 | 0.91 ± 0.11 | 0.537 |
| `D0.r128.HK.mc.d1.s1` | 428770:J12only | 3-5 | 1.87 ± 1.81 (n) | 2.13 ± 1.22 (n) | 2.26 ± 2.64 (n) | – |
| `D0.r128.HK.mc.d1.s1` | 428770:J12only | 5-8 | 1.00 ± 0.18 | 0.94 ± 0.18 | 0.99 ± 0.19 | 0.521 |
| `D0.r128.HK.mc.d1.s1` | 428770:J12only | 8-30 | 1.02 ± 0.05 | 1.02 ± 0.05 | 0.97 ± 0.06 | 0.509 |
| `D0.r128.HK.mc.d1.s1` | 428770:all | 3-5 | 1.78 ± 1.48 (n) | 1.75 ± 0.71 (n) | 1.86 ± 0.86 (n) | – |
| `D0.r128.HK.mc.d1.s1` | 428770:all | 5-8 | 0.91 ± 0.12 | 0.88 ± 0.12 | 0.92 ± 0.13 | 0.541 |
| `D0.r128.HK.mc.d1.s1` | 428770:all | 8-30 | 1.00 ± 0.04 | 0.99 ± 0.04 | 0.94 ± 0.05 | 0.515 |
| `D0.r128.HK.mc.d100.s1` | 427929:MB427929 | 3-5 | 1.26 ± 0.33 (n) | 1.25 ± 0.43 (n) | 1.25 ± 0.44 (n) | – |
| `D0.r128.HK.mc.d100.s1` | 427929:MB427929 | 5-8 | 0.99 ± 0.09 | 0.94 ± 0.10 | 1.00 ± 0.13 | 0.520 |
| `D0.r128.HK.mc.d100.s1` | 427929:MB427929 | 8-30 | 1.05 ± 0.07 | 1.03 ± 0.09 | 1.07 ± 0.12 | 0.521 |
| `D0.r128.HK.mc.d100.s1` | 428770:J12only | 3-5 | 1.46 ± 0.67 (n) | 2.14 ± 2.16 (n) | 2.25 ± 2.64 (n) | – |
| `D0.r128.HK.mc.d100.s1` | 428770:J12only | 5-8 | 0.85 ± 0.12 (n) | 0.94 ± 0.16 | 0.98 ± 0.20 | 0.522 |
| `D0.r128.HK.mc.d100.s1` | 428770:J12only | 8-30 | 0.99 ± 0.03 | 1.03 ± 0.04 | 1.05 ± 0.06 | 0.513 |
| `D0.r128.HK.mc.d100.s1` | 428770:all | 3-5 | 1.61 ± 0.56 (n) | 1.96 ± 1.95 (n) | 1.89 ± 0.88 (n) | – |
| `D0.r128.HK.mc.d100.s1` | 428770:all | 5-8 | 0.86 ± 0.09 | 0.88 ± 0.11 | 0.86 ± 0.12 | 0.555 |
| `D0.r128.HK.mc.d100.s1` | 428770:all | 8-30 | 0.98 ± 0.03 | 1.01 ± 0.04 | 1.01 ± 0.05 | 0.508 |
| `D0.r128.HK.mc.d25.s1` | 427929:MB427929 | 3-5 | 1.23 ± 0.41 (n) | 1.42 ± 0.51 (n) | 1.37 ± 0.48 | 0.593 |
| `D0.r128.HK.mc.d25.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.05 ± 0.13 | 0.99 ± 0.14 | 0.526 |
| `D0.r128.HK.mc.d25.s1` | 427929:MB427929 | 8-30 | 0.96 ± 0.09 | 0.92 ± 0.09 | 0.92 ± 0.11 | 0.528 |
| `D0.r128.HK.mc.d25.s1` | 428770:J12only | 3-5 | 1.87 ± 1.71 (n) | 2.14 ± 2.28 (n) | 2.17 ± 2.59 (n) | – |
| `D0.r128.HK.mc.d25.s1` | 428770:J12only | 5-8 | 0.88 ± 0.16 | 0.97 ± 0.19 | 1.00 ± 0.20 | 0.549 |
| `D0.r128.HK.mc.d25.s1` | 428770:J12only | 8-30 | 0.97 ± 0.05 | 1.00 ± 0.05 | 1.03 ± 0.06 | 0.512 |
| `D0.r128.HK.mc.d25.s1` | 428770:all | 3-5 | 1.80 ± 1.74 (n) | 1.97 ± 2.03 (n) | 1.73 ± 0.82 (n) | – |
| `D0.r128.HK.mc.d25.s1` | 428770:all | 5-8 | 0.81 ± 0.11 | 0.88 ± 0.12 | 0.91 ± 0.13 | 0.575 |
| `D0.r128.HK.mc.d25.s1` | 428770:all | 8-30 | 0.96 ± 0.04 | 0.97 ± 0.04 | 0.99 ± 0.05 | 0.517 |
| `D0.r128.HK.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.35 ± 0.37 (n) | 1.34 ± 0.44 | 1.35 ± 0.47 | 0.620 |
| `D0.r128.HK.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.06 ± 0.12 | 1.06 ± 0.13 | 1.02 ± 0.14 | 0.523 |
| `D0.r128.HK.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.08 | 0.93 ± 0.08 | 1.02 ± 0.12 | 0.528 |
| `D0.r128.HK.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.62 (n) | 2.13 ± 1.39 (n) | 2.24 ± 2.65 (n) | – |
| `D0.r128.HK.mc.d5.s1` | 428770:J12only | 5-8 | 0.97 ± 0.18 | 0.95 ± 0.18 | 0.96 ± 0.19 | 0.519 |
| `D0.r128.HK.mc.d5.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 1.02 ± 0.05 | 0.96 ± 0.06 | 0.509 |
| `D0.r128.HK.mc.d5.s1` | 428770:all | 3-5 | 1.80 ± 1.97 (n) | 1.77 ± 0.71 (n) | 1.85 ± 0.86 (n) | – |
| `D0.r128.HK.mc.d5.s1` | 428770:all | 5-8 | 0.88 ± 0.12 | 0.88 ± 0.12 | 0.89 ± 0.12 | 0.550 |
| `D0.r128.HK.mc.d5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 1.00 ± 0.04 | 0.94 ± 0.05 | 0.516 |
| `D0.r128.HK.mc.m5.s1` | 427929:MB427929 | 3-5 | 1.31 ± 0.38 (n) | 1.38 ± 0.48 (n) | 1.40 ± 0.48 | 0.600 |
| `D0.r128.HK.mc.m5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.06 ± 0.13 | 1.05 ± 0.14 | 0.527 |
| `D0.r128.HK.mc.m5.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.08 | 0.90 ± 0.08 | 0.93 ± 0.11 | 0.535 |
| `D0.r128.HK.mc.m5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.58 (n) | 2.12 ± 1.95 (n) | 2.32 ± 2.66 (n) | – |
| `D0.r128.HK.mc.m5.s1` | 428770:J12only | 5-8 | 0.97 ± 0.18 | 0.93 ± 0.18 | 0.94 ± 0.18 | 0.523 |
| `D0.r128.HK.mc.m5.s1` | 428770:J12only | 8-30 | 1.01 ± 0.05 | 1.01 ± 0.05 | 0.98 ± 0.06 | 0.505 |
| `D0.r128.HK.mc.m5.s1` | 428770:all | 3-5 | 1.87 ± 1.55 (n) | 1.78 ± 0.70 (n) | 1.93 ± 0.88 (n) | – |
| `D0.r128.HK.mc.m5.s1` | 428770:all | 5-8 | 0.91 ± 0.12 | 0.86 ± 0.12 | 0.87 ± 0.12 | 0.550 |
| `D0.r128.HK.mc.m5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 0.99 ± 0.04 | 0.96 ± 0.05 | 0.511 |
| `D0.r128.HK.mc.n.s1` | 427929:MB427929 | 3-5 | 1.41 ± 0.42 (n) | 1.33 ± 0.46 (n) | 1.31 ± 0.45 | 0.579 |
| `D0.r128.HK.mc.n.s1` | 427929:MB427929 | 5-8 | 1.11 ± 0.13 | 1.06 ± 0.13 | 1.03 ± 0.14 | 0.543 |
| `D0.r128.HK.mc.n.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.08 | 0.90 ± 0.08 | 0.97 ± 0.11 | 0.534 |
| `D0.r128.HK.mc.n.s1` | 428770:J12only | 3-5 | 1.87 ± 2.05 (n) | 2.13 ± 2.37 (n) | 2.26 ± 2.63 (n) | – |
| `D0.r128.HK.mc.n.s1` | 428770:J12only | 5-8 | 0.96 ± 0.18 | 0.95 ± 0.18 | 1.00 ± 0.20 | 0.517 |
| `D0.r128.HK.mc.n.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 1.01 ± 0.05 | 0.98 ± 0.06 | 0.504 |
| `D0.r128.HK.mc.n.s1` | 428770:all | 3-5 | 1.68 ± 1.66 (n) | 1.83 ± 1.78 (n) | 1.89 ± 0.87 (n) | – |
| `D0.r128.HK.mc.n.s1` | 428770:all | 5-8 | 0.89 ± 0.12 | 0.89 ± 0.12 | 0.91 ± 0.13 | 0.546 |
| `D0.r128.HK.mc.n.s1` | 428770:all | 8-30 | 0.99 ± 0.04 | 0.98 ± 0.04 | 0.96 ± 0.05 | 0.511 |
| `D0.r128.HK.sb.d5.s1` | 427929:MB427929 | 3-5 | 1.32 ± 0.39 | 1.33 ± 0.45 (n) | 1.33 ± 0.46 | 0.628 |
| `D0.r128.HK.sb.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.03 ± 0.13 | 0.94 ± 0.13 | 0.529 |
| `D0.r128.HK.sb.d5.s1` | 427929:MB427929 | 8-30 | 0.94 ± 0.08 | 0.95 ± 0.09 | 0.97 ± 0.12 | 0.526 |
| `D0.r128.HK.sb.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.90 (n) | 2.00 ± 2.34 (n) | 2.11 ± 2.54 (n) | – |
| `D0.r128.HK.sb.d5.s1` | 428770:J12only | 5-8 | 0.99 ± 0.19 | 0.95 ± 0.18 | 0.95 ± 0.19 | 0.519 |
| `D0.r128.HK.sb.d5.s1` | 428770:J12only | 8-30 | 0.99 ± 0.05 | 0.97 ± 0.05 | 0.88 ± 0.05 | 0.530 |
| `D0.r128.HK.sb.d5.s1` | 428770:all | 3-5 | 1.73 ± 1.22 (n) | 1.69 ± 0.70 (n) | 1.78 ± 0.83 (n) | – |
| `D0.r128.HK.sb.d5.s1` | 428770:all | 5-8 | 0.90 ± 0.12 | 0.88 ± 0.12 | 0.88 ± 0.12 | 0.543 |
| `D0.r128.HK.sb.d5.s1` | 428770:all | 8-30 | 0.97 ± 0.04 | 0.95 ± 0.04 | 0.86 ± 0.05 | 0.535 |
| `D0.r128.HK.sbmc.d5.s1` | 427929:MB427929 | 3-5 | 1.31 ± 0.42 (n) | 1.30 ± 0.41 | 1.39 ± 0.48 | 0.605 |
| `D0.r128.HK.sbmc.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.05 ± 0.13 | 0.98 ± 0.13 | 0.528 |
| `D0.r128.HK.sbmc.d5.s1` | 427929:MB427929 | 8-30 | 0.90 ± 0.08 | 0.91 ± 0.08 | 0.95 ± 0.11 | 0.539 |
| `D0.r128.HK.sbmc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.58 (n) | 2.03 ± 2.36 (n) | 2.10 ± 2.12 (n) | – |
| `D0.r128.HK.sbmc.d5.s1` | 428770:J12only | 5-8 | 0.97 ± 0.18 | 0.96 ± 0.18 | 0.93 ± 0.18 | 0.519 |
| `D0.r128.HK.sbmc.d5.s1` | 428770:J12only | 8-30 | 0.98 ± 0.05 | 0.94 ± 0.05 | 0.76 ± 0.05 | 0.562 |
| `D0.r128.HK.sbmc.d5.s1` | 428770:all | 3-5 | 1.67 ± 1.61 (n) | 1.77 ± 0.70 (n) | 1.77 ± 0.84 (n) | – |
| `D0.r128.HK.sbmc.d5.s1` | 428770:all | 5-8 | 0.91 ± 0.12 | 0.88 ± 0.12 | 0.86 ± 0.12 | 0.541 |
| `D0.r128.HK.sbmc.d5.s1` | 428770:all | 8-30 | 0.96 ± 0.04 | 0.92 ± 0.04 | 0.75 ± 0.04 | 0.563 |
| `D0.r128.HKP.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.36 ± 0.38 (n) | 1.37 ± 0.46 (n) | 1.32 ± 0.46 | 0.582 |
| `D0.r128.HKP.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.05 ± 0.13 | 0.99 ± 0.13 | 0.527 |
| `D0.r128.HKP.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.08 | 0.95 ± 0.09 | 1.06 ± 0.13 | 0.527 |
| `D0.r128.HKP.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.59 (n) | 2.13 ± 1.38 (n) | 2.01 ± 1.25 (n) | – |
| `D0.r128.HKP.mc.d5.s1` | 428770:J12only | 5-8 | 0.96 ± 0.18 | 0.93 ± 0.18 | 0.98 ± 0.19 | 0.525 |
| `D0.r128.HKP.mc.d5.s1` | 428770:J12only | 8-30 | 1.01 ± 0.05 | 1.03 ± 0.05 | 0.96 ± 0.05 | 0.511 |
| `D0.r128.HKP.mc.d5.s1` | 428770:all | 3-5 | 1.81 ± 2.01 (n) | 1.73 ± 0.68 (n) | 1.76 ± 0.82 (n) | – |
| `D0.r128.HKP.mc.d5.s1` | 428770:all | 5-8 | 0.88 ± 0.12 | 0.88 ± 0.12 | 0.90 ± 0.13 | 0.550 |
| `D0.r128.HKP.mc.d5.s1` | 428770:all | 8-30 | 0.99 ± 0.04 | 1.01 ± 0.04 | 0.94 ± 0.05 | 0.516 |
| `D0.r128.HKia.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.39 ± 0.40 (n) | 1.42 ± 0.48 | 1.35 ± 0.48 (n) | 0.650 |
| `D0.r128.HKia.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.04 ± 0.13 | 1.05 ± 0.14 | 0.527 |
| `D0.r128.HKia.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.94 ± 0.08 | 0.92 ± 0.09 | 0.97 ± 0.11 | 0.530 |
| `D0.r128.HKia.mc.d5.s1` | 428770:J12only | 3-5 | 1.73 ± 2.00 (n) | 1.89 ± 2.19 (n) | 1.84 ± 1.16 (n) | – |
| `D0.r128.HKia.mc.d5.s1` | 428770:J12only | 5-8 | 0.88 ± 0.17 | 0.88 ± 0.17 | 0.84 ± 0.17 | 0.550 |
| `D0.r128.HKia.mc.d5.s1` | 428770:J12only | 8-30 | 0.96 ± 0.05 | 0.96 ± 0.05 | 0.87 ± 0.05 | 0.533 |
| `D0.r128.HKia.mc.d5.s1` | 428770:all | 3-5 | 1.56 ± 0.60 (n) | 1.60 ± 0.69 (n) | 1.61 ± 0.76 (n) | – |
| `D0.r128.HKia.mc.d5.s1` | 428770:all | 5-8 | 0.80 ± 0.11 | 0.81 ± 0.11 | 0.77 ± 0.11 | 0.579 |
| `D0.r128.HKia.mc.d5.s1` | 428770:all | 8-30 | 0.95 ± 0.04 | 0.94 ± 0.04 | 0.85 ± 0.04 | 0.538 |
| `D0.r128.HKir.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.40 ± 0.40 (n) | 1.43 ± 0.48 | 1.38 ± 0.50 (n) | 0.652 |
| `D0.r128.HKir.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.05 ± 0.13 | 1.05 ± 0.13 | 1.01 ± 0.14 | 0.521 |
| `D0.r128.HKir.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.95 ± 0.08 | 0.91 ± 0.09 | 0.98 ± 0.12 | 0.532 |
| `D0.r128.HKir.mc.d5.s1` | 428770:J12only | 3-5 | 1.75 ± 1.77 (n) | 1.95 ± 2.26 (n) | 1.83 ± 1.17 (n) | – |
| `D0.r128.HKir.mc.d5.s1` | 428770:J12only | 5-8 | 0.87 ± 0.17 | 0.87 ± 0.17 | 0.85 ± 0.17 | 0.553 |
| `D0.r128.HKir.mc.d5.s1` | 428770:J12only | 8-30 | 0.97 ± 0.05 | 0.95 ± 0.05 | 0.84 ± 0.05 | 0.541 |
| `D0.r128.HKir.mc.d5.s1` | 428770:all | 3-5 | 1.56 ± 0.61 (n) | 1.62 ± 0.71 (n) | 1.60 ± 0.75 (n) | – |
| `D0.r128.HKir.mc.d5.s1` | 428770:all | 5-8 | 0.80 ± 0.11 | 0.80 ± 0.11 | 0.77 ± 0.11 | 0.581 |
| `D0.r128.HKir.mc.d5.s1` | 428770:all | 8-30 | 0.95 ± 0.04 | 0.94 ± 0.04 | 0.83 ± 0.04 | 0.544 |
| `D0.r128.HX.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.39 ± 0.39 (n) | 1.34 ± 0.43 | 1.30 ± 0.44 | 0.620 |
| `D0.r128.HX.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.09 ± 0.13 | 1.05 ± 0.13 | 1.01 ± 0.14 | 0.537 |
| `D0.r128.HX.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.90 ± 0.07 | 0.93 ± 0.08 | 1.06 ± 0.13 | 0.540 |
| `D0.r128.HX.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.68 (n) | 2.13 ± 2.26 (n) | 2.13 ± 2.57 (n) | – |
| `D0.r128.HX.mc.d5.s1` | 428770:J12only | 5-8 | 0.92 ± 0.17 | 0.94 ± 0.18 | 0.96 ± 0.19 | 0.531 |
| `D0.r128.HX.mc.d5.s1` | 428770:J12only | 8-30 | 1.01 ± 0.05 | 1.03 ± 0.05 | 0.99 ± 0.06 | 0.511 |
| `D0.r128.HX.mc.d5.s1` | 428770:all | 3-5 | 1.75 ± 1.72 (n) | 1.84 ± 1.77 (n) | 1.75 ± 0.81 (n) | – |
| `D0.r128.HX.mc.d5.s1` | 428770:all | 5-8 | 0.85 ± 0.11 | 0.88 ± 0.12 | 0.88 ± 0.12 | 0.561 |
| `D0.r128.HX.mc.d5.s1` | 428770:all | 8-30 | 0.99 ± 0.04 | 1.01 ± 0.04 | 0.96 ± 0.05 | 0.509 |
| `D0.tr.H+dau.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.47 ± 0.41 (n) | 1.36 ± 0.46 | 1.37 ± 0.49 (n) | 0.628 |
| `D0.tr.H+dau.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.07 ± 0.13 | 1.04 ± 0.13 | 1.00 ± 0.14 | 0.529 |
| `D0.tr.H+dau.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.92 ± 0.08 | 0.93 ± 0.08 | 0.97 ± 0.11 | 0.533 |
| `D0.tr.H+dau.mc.d5.s1` | 428770:J12only | 3-5 | 1.87 ± 1.62 (n) | 2.13 ± 2.07 (n) | 2.27 ± 2.71 (n) | – |
| `D0.tr.H+dau.mc.d5.s1` | 428770:J12only | 5-8 | 0.91 ± 0.17 | 0.96 ± 0.19 | 1.02 ± 0.20 | 0.536 |
| `D0.tr.H+dau.mc.d5.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 1.01 ± 0.05 | 1.06 ± 0.06 | 0.515 |
| `D0.tr.H+dau.mc.d5.s1` | 428770:all | 3-5 | 1.81 ± 1.83 (n) | 1.90 ± 1.91 (n) | 1.83 ± 0.87 (n) | – |
| `D0.tr.H+dau.mc.d5.s1` | 428770:all | 5-8 | 0.87 ± 0.11 | 0.89 ± 0.12 | 0.93 ± 0.13 | 0.553 |
| `D0.tr.H+dau.mc.d5.s1` | 428770:all | 8-30 | 0.98 ± 0.04 | 0.99 ± 0.04 | 1.02 ± 0.05 | 0.509 |
| `D0.tr.H+e04.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.43 ± 0.33 (n) | 1.62 ± 1.38 (n) | 1.44 ± 0.49 | 0.649 |
| `D0.tr.H+e04.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.03 ± 0.11 | 1.02 ± 0.11 | 1.07 ± 0.13 | 0.525 |
| `D0.tr.H+e04.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.92 ± 0.06 | 0.97 ± 0.06 | 0.90 ± 0.08 | 0.538 |
| `D0.tr.H+e04.mc.d5.s1` | 428770:J12only | 3-5 | 1.61 ± 1.64 (n) | 1.48 ± 1.52 (n) | 1.40 ± 0.75 (n) | – |
| `D0.tr.H+e04.mc.d5.s1` | 428770:J12only | 5-8 | 0.97 ± 0.18 | 0.93 ± 0.17 | 0.81 ± 0.16 | 0.568 |
| `D0.tr.H+e04.mc.d5.s1` | 428770:J12only | 8-30 | 0.93 ± 0.04 | 0.90 ± 0.04 | 0.83 ± 0.04 | 0.563 |
| `D0.tr.H+e04.mc.d5.s1` | 428770:all | 3-5 | 1.61 ± 1.41 (n) | 1.52 ± 0.51 (n) | 1.36 ± 0.57 (n) | – |
| `D0.tr.H+e04.mc.d5.s1` | 428770:all | 5-8 | 0.88 ± 0.11 | 0.84 ± 0.11 | 0.75 ± 0.10 | 0.588 |
| `D0.tr.H+e04.mc.d5.s1` | 428770:all | 8-30 | 0.93 ± 0.04 | 0.90 ± 0.04 | 0.83 ± 0.04 | 0.563 |
| `D0.tr.H+e04.mc.n.s1` | 427929:MB427929 | 3-5 | 1.61 ± 1.29 (n) | 1.73 ± 1.76 (n) | 1.57 ± 0.46 (n) | – |
| `D0.tr.H+e04.mc.n.s1` | 427929:MB427929 | 5-8 | 1.14 ± 0.12 | 1.08 ± 0.12 | 1.09 ± 0.13 | 0.566 |
| `D0.tr.H+e04.mc.n.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.06 | 0.92 ± 0.06 | 0.93 ± 0.08 | 0.537 |
| `D0.tr.H+e04.mc.n.s1` | 428770:J12only | 3-5 | 1.61 ± 0.35 (n) | 1.73 ± 1.69 (n) | 1.46 ± 0.77 (n) | – |
| `D0.tr.H+e04.mc.n.s1` | 428770:J12only | 5-8 | 0.96 ± 0.17 | 0.90 ± 0.17 | 0.88 ± 0.17 | 0.544 |
| `D0.tr.H+e04.mc.n.s1` | 428770:J12only | 8-30 | 0.94 ± 0.04 | 0.93 ± 0.04 | 0.90 ± 0.05 | 0.539 |
| `D0.tr.H+e04.mc.n.s1` | 428770:all | 3-5 | 1.61 ± 0.23 (n) | 1.73 ± 1.29 (n) | 1.44 ± 0.59 (n) | – |
| `D0.tr.H+e04.mc.n.s1` | 428770:all | 5-8 | 0.87 ± 0.11 | 0.83 ± 0.11 | 0.81 ± 0.11 | 0.576 |
| `D0.tr.H+e04.mc.n.s1` | 428770:all | 8-30 | 0.93 ± 0.04 | 0.93 ± 0.04 | 0.90 ± 0.04 | 0.537 |
| `D0.tr.H+e04.sbmc.d5.s1` | 427929:MB427929 | 3-5 | 1.35 ± 0.34 (n) | 1.50 ± 0.40 (n) | 1.38 ± 0.49 (n) | – |
| `D0.tr.H+e04.sbmc.d5.s1` | 427929:MB427929 | 5-8 | 1.09 ± 0.12 | 1.03 ± 0.12 | 0.99 ± 0.12 | 0.543 |
| `D0.tr.H+e04.sbmc.d5.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.07 | 0.93 ± 0.07 | 0.91 ± 0.08 | 0.531 |
| `D0.tr.H+e04.sbmc.d5.s1` | 428770:J12only | 3-5 | 1.34 ± 1.07 (n) | 1.28 ± 0.68 (n) | 1.33 ± 0.77 (n) | – |
| `D0.tr.H+e04.sbmc.d5.s1` | 428770:J12only | 5-8 | 0.88 ± 0.16 | 0.84 ± 0.16 | 0.76 ± 0.15 | 0.584 |
| `D0.tr.H+e04.sbmc.d5.s1` | 428770:J12only | 8-30 | 0.85 ± 0.04 | 0.77 ± 0.04 | 0.62 ± 0.04 | 0.629 |
| `D0.tr.H+e04.sbmc.d5.s1` | 428770:all | 3-5 | 1.39 ± 0.49 (n) | 1.27 ± 0.51 (n) | 1.32 ± 0.57 (n) | – |
| `D0.tr.H+e04.sbmc.d5.s1` | 428770:all | 5-8 | 0.82 ± 0.11 | 0.76 ± 0.10 | 0.70 ± 0.10 | 0.605 |
| `D0.tr.H+e04.sbmc.d5.s1` | 428770:all | 8-30 | 0.86 ± 0.04 | 0.79 ± 0.04 | 0.65 ± 0.03 | 0.618 |
| `D0.tr.H+e10.mc.d5.s1` | 427929:MB427929 | 3-5 | 1.33 ± 0.34 (n) | 1.38 ± 0.39 (n) | 1.46 ± 0.52 (n) | – |
| `D0.tr.H+e10.mc.d5.s1` | 427929:MB427929 | 5-8 | 1.03 ± 0.11 | 0.99 ± 0.11 | 1.04 ± 0.13 | 0.515 |
| `D0.tr.H+e10.mc.d5.s1` | 427929:MB427929 | 8-30 | 0.91 ± 0.06 | 0.95 ± 0.07 | 0.97 ± 0.10 | 0.541 |
| `D0.tr.H+e10.mc.d5.s1` | 428770:J12only | 3-5 | 1.60 ± 1.09 (n) | 1.42 ± 1.48 (n) | 1.77 ± 1.99 (n) | – |
| `D0.tr.H+e10.mc.d5.s1` | 428770:J12only | 5-8 | 0.90 ± 0.17 | 0.85 ± 0.16 | 0.82 ± 0.16 | 0.568 |
| `D0.tr.H+e10.mc.d5.s1` | 428770:J12only | 8-30 | 0.92 ± 0.04 | 0.89 ± 0.04 | 0.79 ± 0.04 | 0.578 |
| `D0.tr.H+e10.mc.d5.s1` | 428770:all | 3-5 | 1.52 ± 1.50 (n) | 1.37 ± 0.52 (n) | 1.55 ± 0.68 (n) | – |
| `D0.tr.H+e10.mc.d5.s1` | 428770:all | 5-8 | 0.82 ± 0.11 | 0.76 ± 0.10 | 0.75 ± 0.10 | 0.607 |
| `D0.tr.H+e10.mc.d5.s1` | 428770:all | 8-30 | 0.92 ± 0.04 | 0.89 ± 0.04 | 0.79 ± 0.04 | 0.576 |
| `D0.xgb.H.mc.n.s1` | 427929:MB427929 | 3-5 | 1.37 ± 0.39 (n) | 1.28 ± 0.43 (n) | 1.27 ± 0.43 | 0.568 |
| `D0.xgb.H.mc.n.s1` | 427929:MB427929 | 5-8 | 1.06 ± 0.13 | 1.01 ± 0.13 | 0.99 ± 0.14 | 0.524 |
| `D0.xgb.H.mc.n.s1` | 427929:MB427929 | 8-30 | 0.95 ± 0.09 | 0.95 ± 0.09 | 0.96 ± 0.12 | 0.520 |
| `D0.xgb.H.mc.n.s1` | 428770:J12only | 3-5 | 1.87 ± 1.71 (n) | 2.13 ± 2.32 (n) | 2.45 ± 2.90 (n) | – |
| `D0.xgb.H.mc.n.s1` | 428770:J12only | 5-8 | 0.94 ± 0.17 | 1.01 ± 0.19 | 1.04 ± 0.21 | 0.525 |
| `D0.xgb.H.mc.n.s1` | 428770:J12only | 8-30 | 1.06 ± 0.05 | 1.07 ± 0.05 | 1.14 ± 0.06 | 0.536 |
| `D0.xgb.H.mc.n.s1` | 428770:all | 3-5 | 1.82 ± 2.07 (n) | 1.83 ± 1.82 (n) | 1.93 ± 0.91 (n) | – |
| `D0.xgb.H.mc.n.s1` | 428770:all | 5-8 | 0.87 ± 0.11 | 0.91 ± 0.12 | 0.97 ± 0.14 | 0.553 |
| `D0.xgb.H.mc.n.s1` | 428770:all | 8-30 | 1.03 ± 0.04 | 1.04 ± 0.05 | 1.10 ± 0.05 | 0.525 |
| `D0.xgb.H.sb.n.s1` | 427929:MB427929 | 3-5 | 1.40 ± 0.39 (n) | 1.31 ± 0.43 (n) | 1.28 ± 0.43 | 0.571 |
| `D0.xgb.H.sb.n.s1` | 427929:MB427929 | 5-8 | 1.08 ± 0.13 | 1.01 ± 0.13 | 0.95 ± 0.13 | 0.532 |
| `D0.xgb.H.sb.n.s1` | 427929:MB427929 | 8-30 | 0.94 ± 0.08 | 0.97 ± 0.10 | 0.97 ± 0.12 | 0.526 |
| `D0.xgb.H.sb.n.s1` | 428770:J12only | 3-5 | 1.87 ± 2.12 (n) | 2.13 ± 1.23 (n) | 2.47 ± 2.89 (n) | – |
| `D0.xgb.H.sb.n.s1` | 428770:J12only | 5-8 | 0.92 ± 0.17 | 0.93 ± 0.18 | 1.00 ± 0.20 | 0.534 |
| `D0.xgb.H.sb.n.s1` | 428770:J12only | 8-30 | 0.99 ± 0.05 | 0.99 ± 0.05 | 0.98 ± 0.06 | 0.505 |
| `D0.xgb.H.sb.n.s1` | 428770:all | 3-5 | 1.81 ± 2.02 (n) | 1.85 ± 1.30 (n) | 2.02 ± 0.95 (n) | – |
| `D0.xgb.H.sb.n.s1` | 428770:all | 5-8 | 0.85 ± 0.11 | 0.88 ± 0.12 | 0.92 ± 0.13 | 0.561 |
| `D0.xgb.H.sb.n.s1` | 428770:all | 8-30 | 0.97 ± 0.04 | 0.97 ± 0.04 | 0.95 ± 0.05 | 0.512 |
| `D0.xgb.HK.mc.n.s1` | 427929:MB427929 | 3-5 | 1.35 ± 0.39 (n) | 1.33 ± 0.44 | 1.36 ± 0.48 (n) | 0.617 |
| `D0.xgb.HK.mc.n.s1` | 427929:MB427929 | 5-8 | 1.02 ± 0.12 | 1.06 ± 0.13 | 0.99 ± 0.13 | 0.522 |
| `D0.xgb.HK.mc.n.s1` | 427929:MB427929 | 8-30 | 0.97 ± 0.09 | 0.92 ± 0.08 | 0.99 ± 0.12 | 0.527 |
| `D0.xgb.HK.mc.n.s1` | 428770:J12only | 3-5 | 1.87 ± 1.63 (n) | 2.13 ± 2.28 (n) | 2.06 ± 2.50 (n) | – |
| `D0.xgb.HK.mc.n.s1` | 428770:J12only | 5-8 | 0.96 ± 0.18 | 0.98 ± 0.19 | 0.99 ± 0.20 | 0.516 |
| `D0.xgb.HK.mc.n.s1` | 428770:J12only | 8-30 | 1.05 ± 0.05 | 1.09 ± 0.05 | 1.12 ± 0.06 | 0.533 |
| `D0.xgb.HK.mc.n.s1` | 428770:all | 3-5 | 1.85 ± 2.06 (n) | 1.91 ± 1.87 (n) | 1.75 ± 0.81 (n) | – |
| `D0.xgb.HK.mc.n.s1` | 428770:all | 5-8 | 0.88 ± 0.12 | 0.90 ± 0.12 | 0.92 ± 0.13 | 0.550 |
| `D0.xgb.HK.mc.n.s1` | 428770:all | 8-30 | 1.03 ± 0.04 | 1.06 ± 0.05 | 1.08 ± 0.05 | 0.521 |
| `D0.xgb.HK.sb.n.s1` | 427929:MB427929 | 3-5 | 1.37 ± 0.37 (n) | 1.29 ± 0.42 | 1.28 ± 0.44 (n) | 0.603 |
| `D0.xgb.HK.sb.n.s1` | 427929:MB427929 | 5-8 | 1.00 ± 0.12 | 1.06 ± 0.14 | 0.94 ± 0.13 | 0.523 |
| `D0.xgb.HK.sb.n.s1` | 427929:MB427929 | 8-30 | 0.94 ± 0.08 | 1.00 ± 0.10 | 1.05 ± 0.13 | 0.525 |
| `D0.xgb.HK.sb.n.s1` | 428770:J12only | 3-5 | 1.87 ± 1.77 (n) | 2.13 ± 2.10 (n) | 2.37 ± 1.89 (n) | – |
| `D0.xgb.HK.sb.n.s1` | 428770:J12only | 5-8 | 0.96 ± 0.18 | 0.96 ± 0.18 | 1.00 ± 0.20 | 0.516 |
| `D0.xgb.HK.sb.n.s1` | 428770:J12only | 8-30 | 1.00 ± 0.05 | 0.99 ± 0.05 | 0.96 ± 0.06 | 0.510 |
| `D0.xgb.HK.sb.n.s1` | 428770:all | 3-5 | 1.87 ± 1.20 (n) | 1.83 ± 1.78 (n) | 1.95 ± 0.92 (n) | – |
| `D0.xgb.HK.sb.n.s1` | 428770:all | 5-8 | 0.89 ± 0.12 | 0.89 ± 0.12 | 0.91 ± 0.13 | 0.546 |
| `D0.xgb.HK.sb.n.s1` | 428770:all | 8-30 | 0.99 ± 0.04 | 0.97 ± 0.04 | 0.93 ± 0.05 | 0.518 |
| `D0.xgb.HKP.mc.n.s1` | 427929:MB427929 | 3-5 | 1.34 ± 0.39 (n) | 1.36 ± 0.48 (n) | 1.31 ± 0.45 | 0.578 |
| `D0.xgb.HKP.mc.n.s1` | 427929:MB427929 | 5-8 | 1.02 ± 0.12 | 1.06 ± 0.13 | 1.01 ± 0.14 | 0.523 |
| `D0.xgb.HKP.mc.n.s1` | 427929:MB427929 | 8-30 | 0.95 ± 0.09 | 0.96 ± 0.08 | 1.02 ± 0.13 | 0.519 |
| `D0.xgb.HKP.mc.n.s1` | 428770:J12only | 3-5 | 1.87 ± 1.60 (n) | 2.13 ± 2.21 (n) | 2.11 ± 2.56 (n) | – |
| `D0.xgb.HKP.mc.n.s1` | 428770:J12only | 5-8 | 0.93 ± 0.17 | 0.98 ± 0.19 | 0.97 ± 0.19 | 0.526 |
| `D0.xgb.HKP.mc.n.s1` | 428770:J12only | 8-30 | 1.05 ± 0.05 | 1.08 ± 0.05 | 1.14 ± 0.06 | 0.536 |
| `D0.xgb.HKP.mc.n.s1` | 428770:all | 3-5 | 1.83 ± 2.09 (n) | 1.83 ± 1.80 (n) | 1.80 ± 0.85 (n) | – |
| `D0.xgb.HKP.mc.n.s1` | 428770:all | 5-8 | 0.86 ± 0.11 | 0.91 ± 0.12 | 0.90 ± 0.13 | 0.556 |
| `D0.xgb.HKP.mc.n.s1` | 428770:all | 8-30 | 1.03 ± 0.04 | 1.05 ± 0.05 | 1.10 ± 0.05 | 0.526 |
| `D0.xgb.HX.mc.n.s1` | 427929:MB427929 | 3-5 | 1.37 ± 0.39 (n) | 1.34 ± 0.44 | 1.33 ± 0.47 (n) | 0.620 |
| `D0.xgb.HX.mc.n.s1` | 427929:MB427929 | 5-8 | 1.04 ± 0.12 | 1.04 ± 0.13 | 1.01 ± 0.14 | 0.516 |
| `D0.xgb.HX.mc.n.s1` | 427929:MB427929 | 8-30 | 0.93 ± 0.08 | 0.96 ± 0.08 | 1.03 ± 0.13 | 0.527 |
| `D0.xgb.HX.mc.n.s1` | 428770:J12only | 3-5 | 1.87 ± 1.81 (n) | 2.14 ± 2.31 (n) | 2.17 ± 2.60 (n) | – |
| `D0.xgb.HX.mc.n.s1` | 428770:J12only | 5-8 | 0.93 ± 0.17 | 0.98 ± 0.19 | 0.99 ± 0.19 | 0.530 |
| `D0.xgb.HX.mc.n.s1` | 428770:J12only | 8-30 | 1.05 ± 0.05 | 1.08 ± 0.05 | 1.14 ± 0.06 | 0.535 |
| `D0.xgb.HX.mc.n.s1` | 428770:all | 3-5 | 1.77 ± 1.42 (n) | 1.91 ± 1.89 (n) | 1.82 ± 0.85 (n) | – |
| `D0.xgb.HX.mc.n.s1` | 428770:all | 5-8 | 0.86 ± 0.11 | 0.90 ± 0.12 | 0.92 ± 0.13 | 0.556 |
| `D0.xgb.HX.mc.n.s1` | 428770:all | 8-30 | 1.03 ± 0.04 | 1.04 ± 0.05 | 1.10 ± 0.05 | 0.525 |
| `D0.xgb.HX.sb.n.s1` | 427929:MB427929 | 3-5 | 1.31 ± 0.39 | 1.33 ± 0.43 | 1.31 ± 0.46 (n) | 0.623 |
| `D0.xgb.HX.sb.n.s1` | 427929:MB427929 | 5-8 | 1.12 ± 0.13 | 1.06 ± 0.13 | 0.95 ± 0.13 | 0.547 |
| `D0.xgb.HX.sb.n.s1` | 427929:MB427929 | 8-30 | 0.92 ± 0.08 | 0.98 ± 0.09 | 0.97 ± 0.12 | 0.531 |
| `D0.xgb.HX.sb.n.s1` | 428770:J12only | 3-5 | 1.87 ± 2.11 (n) | 2.13 ± 1.21 (n) | 2.23 ± 2.62 (n) | – |
| `D0.xgb.HX.sb.n.s1` | 428770:J12only | 5-8 | 0.92 ± 0.17 | 0.94 ± 0.18 | 0.98 ± 0.19 | 0.531 |
| `D0.xgb.HX.sb.n.s1` | 428770:J12only | 8-30 | 0.99 ± 0.05 | 0.99 ± 0.05 | 0.98 ± 0.06 | 0.505 |
| `D0.xgb.HX.sb.n.s1` | 428770:all | 3-5 | 1.72 ± 1.18 (n) | 1.78 ± 1.77 (n) | 1.85 ± 0.86 (n) | – |
| `D0.xgb.HX.sb.n.s1` | 428770:all | 5-8 | 0.85 ± 0.11 | 0.87 ± 0.12 | 0.91 ± 0.13 | 0.562 |
| `D0.xgb.HX.sb.n.s1` | 428770:all | 8-30 | 0.97 ± 0.04 | 0.97 ± 0.04 | 0.95 ± 0.05 | 0.514 |

#### Werdykt wg kryteriów §4

- `D0.ds.H+dau.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.ds.H+e04.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **warunkowo**
- `D0.ens.xgb.H.sb.n.s1+r128.HK.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.ens.xgb.HX.mc.n.s1+r128.HK.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.m64.H.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **nie** → **odrzucona**
- `D0.m64.H.mc.n.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.m64.H.sb.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.m64.HK.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.H.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HK.mc.d1.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **nie**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HK.mc.d100.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HK.mc.d25.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **nie** → **odrzucona**
- `D0.r128.HK.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **tak** → **odrzucona**
- `D0.r128.HK.mc.m5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **nie**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HK.mc.n.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **nie**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HK.sb.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (1/3); K2 rzeźbienie nie gorsze niż H: **nie**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HK.sbmc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **nie**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HKP.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.r128.HKia.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **warunkowo**
- `D0.r128.HKir.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **warunkowo**
- `D0.r128.HX.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.tr.H+dau.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **tak** → **odrzucona**
- `D0.tr.H+e04.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **nie**; K4 seedy: **tak** → **odrzucona**
- `D0.tr.H+e04.mc.n.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **tak**; K4 seedy: **brak (1 seed)** → **warunkowo**
- `D0.tr.H+e04.sbmc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **nie**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `D0.tr.H+e10.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **nie**; K4 seedy: **brak (1 seed)** → **odrzucona**

### Dplus — odniesienie `Dplus.xgb.H.sb.n.s1`; wierszy: sig 135929, bmc 556336, bmc_all 3158876, sb428770 953050, dall428770 1658241, sb427929 333939, dall427929 580222

#### Zproxy względem XGB H (G = Z/Z_H, bootstrap sparowany), WP70; R_sculpt i corr na tle MC

| konfiguracja | cech | G_mc 3–5 | G_mc 5–8 | G_mc 8–30 | G_tr(SB 428770) 3–5 | 5–8 | 8–30 | R_sculpt MC (WP70) | R dane (pasmo) | nachylenie /100 MeV | corr(s, m) | par. | czas [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `Dplus.xgb.H.mc.n.s1` | 18 | 0.984 ± 0.004 | 1.023 ± 0.013 | 0.921 ± 0.032 | 0.955 ± 0.004 | 0.979 ± 0.010 | 0.856 ± 0.014 | 0.913 ± 0.014 | 0.980 ± 0.022 | -0.23 | -0.160 | 1140 | 167.2 |
| `Dplus.xgb.H.sb.n.s1` | 18 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 0.904 ± 0.014 | 0.970 ± 0.022 | -0.24 | -0.174 | 1194 | 223.5 |
| `Dplus.xgb.HDLZ.mc.n.s1` | 26 | 1.274 ± 0.010 | 1.184 ± 0.021 | 1.120 ± 0.050 | 1.394 ± 0.014 | 1.266 ± 0.021 | 0.989 ± 0.019 | 0.977 ± 0.020 | 1.024 ± 0.029 | -0.28 | -0.221 | 1189 | 256.9 |
| `Dplus.xgb.HDLZ.sb.n.s1` | 26 | 1.253 ± 0.010 | 1.151 ± 0.019 | 1.150 ± 0.043 | 1.442 ± 0.015 | 1.309 ± 0.024 | 1.236 ± 0.027 | 1.001 ± 0.020 | 1.095 ± 0.032 | -0.30 | -0.263 | 1363 | 306.0 |
| `Dplus.ens.xgb.H.sb.n.s1+r128.HDLZ.mc.d5.s1` | – | 1.230 ± 0.008 | 1.161 ± 0.017 | 1.098 ± 0.033 | 1.245 ± 0.008 | 1.161 ± 0.017 | 1.089 ± 0.015 | 0.936 ± 0.018 | 1.013 ± 0.028 | -0.24 | -0.088 | – | 0 |
| `Dplus.ens.xgb.HDLZ.mc.n.s1+r128.HDLZ.mc.d5.s1` | – | 1.332 ± 0.011 | 1.254 ± 0.024 | 1.176 ± 0.045 | 1.436 ± 0.013 | 1.296 ± 0.023 | 1.104 ± 0.021 | 0.958 ± 0.020 | 1.050 ± 0.031 | -0.27 | -0.132 | – | 0 |
| `Dplus.m64.H.mc.d5.s1` | 18 | 0.912 ± 0.005 | 0.965 ± 0.015 | 1.000 ± 0.034 | 0.860 ± 0.006 | 0.941 ± 0.012 | 0.908 ± 0.014 | 0.912 ± 0.013 | 0.991 ± 0.021 | -0.21 | -0.019 | 9601 | 1017.0 |
| `Dplus.r128.HDLZ.mc.d5.s1` | 26 | 1.260 ± 0.011 | 1.196 ± 0.024 | 1.098 ± 0.046 | 1.240 ± 0.011 | 1.168 ± 0.021 | 1.035 ± 0.018 | 0.937 ± 0.019 | 1.023 ± 0.028 | -0.23 | 0.016 | 54145 | 3232.9 |
| `Dplus.r128.HDLZ.mc.d5.s2` | 26 | 1.245 ± 0.010 | 1.201 ± 0.023 | 1.098 ± 0.043 | 1.234 ± 0.011 | 1.149 ± 0.020 | 1.070 ± 0.018 | 0.921 ± 0.018 | 1.012 ± 0.028 | -0.24 | 0.021 | 54145 | 1999.2 |
| `Dplus.r128.HDLZ.mc.d5.s3` | 26 | 1.263 ± 0.011 | 1.207 ± 0.023 | 1.126 ± 0.047 | 1.236 ± 0.010 | 1.156 ± 0.020 | 1.044 ± 0.018 | 0.920 ± 0.018 | 0.990 ± 0.027 | -0.24 | 0.013 | 54145 | 1580.2 |
| `Dplus.r128.HDLZ.mc.n.s1` | 26 | 1.285 ± 0.012 | 1.123 ± 0.025 | 0.988 ± 0.045 | 1.833 ± 0.025 | 1.597 ± 0.040 | 1.181 ± 0.026 | 1.447 ± 0.031 | 1.413 ± 0.046 | -0.29 | -0.233 | 54145 | 998.0 |
| `Dplus.r128.HDLZ.sbmc.d5.s1` | 26 | 1.238 ± 0.011 | 1.156 ± 0.021 | 1.104 ± 0.044 | 1.431 ± 0.015 | 1.299 ± 0.022 | 1.118 ± 0.022 | 1.047 ± 0.021 | 1.136 ± 0.033 | -0.28 | -0.014 | 54145 | 2033.1 |

#### WP80 / WP50 (G_mc; G_tr)

| konfiguracja | WP80 3–5 | 5–8 | 8–30 | WP50 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|---|---|
| `Dplus.ens.xgb.H.sb.n.s1+r128.HDLZ.mc.d5.s1` | 1.20; 1.21 | 1.14; 1.16 | 1.12; 1.10 | 1.24; 1.28 | 1.24; 1.24 | 1.24; 1.16 |
| `Dplus.ens.xgb.HDLZ.mc.n.s1+r128.HDLZ.mc.d5.s1` | 1.30; 1.40 | 1.21; 1.27 | 1.14; 1.12 | 1.34; 1.47 | 1.34; 1.31 | 1.38; 1.19 |
| `Dplus.m64.H.mc.d5.s1` | 0.91; 0.83 | 0.95; 0.93 | 0.98; 0.92 | 0.92; 0.89 | 1.04; 0.98 | 1.03; 0.93 |
| `Dplus.r128.HDLZ.mc.d5.s1` | 1.22; 1.16 | 1.13; 1.13 | 1.10; 1.05 | 1.26; 1.31 | 1.26; 1.24 | 1.35; 1.12 |
| `Dplus.r128.HDLZ.mc.d5.s2` | 1.21; 1.16 | 1.13; 1.13 | 1.12; 1.06 | 1.24; 1.29 | 1.33; 1.23 | 1.24; 1.09 |
| `Dplus.r128.HDLZ.mc.d5.s3` | 1.22; 1.16 | 1.15; 1.13 | 1.12; 1.05 | 1.28; 1.31 | 1.30; 1.22 | 1.35; 1.10 |
| `Dplus.r128.HDLZ.mc.n.s1` | 1.26; 1.79 | 1.11; 1.61 | 1.04; 1.25 | 1.30; 1.82 | 1.28; 1.58 | 1.22; 1.25 |
| `Dplus.r128.HDLZ.sbmc.d5.s1` | 1.23; 1.37 | 1.14; 1.28 | 1.08; 1.13 | 1.27; 1.48 | 1.22; 1.26 | 1.14; 1.12 |
| `Dplus.xgb.H.mc.n.s1` | 1.00; 0.96 | 1.02; 0.98 | 0.99; 0.89 | 0.98; 0.95 | 1.04; 0.97 | 0.93; 0.84 |
| `Dplus.xgb.H.sb.n.s1` | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 | 1.00; 1.00 |
| `Dplus.xgb.HDLZ.mc.n.s1` | 1.26; 1.38 | 1.19; 1.26 | 1.10; 1.06 | 1.26; 1.40 | 1.32; 1.27 | 1.15; 1.01 |
| `Dplus.xgb.HDLZ.sb.n.s1` | 1.23; 1.42 | 1.15; 1.31 | 1.09; 1.24 | 1.27; 1.44 | 1.27; 1.31 | 1.29; 1.31 |

#### Rozrzut seedów (odchylenie std G_mc WP70 między seedami)

| rodzina | seedy | 3–5 | 5–8 | 8–30 |
|---|---|---|---|---|
| `Dplus.r128.HDLZ.mc.d5` | 3 | 1.256 ± 0.010 | 1.201 ± 0.005 | 1.107 ± 0.016 |

#### Werdykt wg kryteriów §4

- `Dplus.ens.xgb.H.sb.n.s1+r128.HDLZ.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (2/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **brak pomiaru**; K4 seedy: **brak (1 seed)** → **warunkowo**
- `Dplus.ens.xgb.HDLZ.mc.n.s1+r128.HDLZ.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **brak pomiaru**; K4 seedy: **brak (1 seed)** → **warunkowo**
- `Dplus.m64.H.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **nie** (0/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **brak pomiaru**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `Dplus.r128.HDLZ.mc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (2/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **brak pomiaru**; K4 seedy: **tak** → **warunkowo**
- `Dplus.r128.HDLZ.mc.n.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (2/3); K2 rzeźbienie nie gorsze niż H: **nie**; K3 SF: **brak pomiaru**; K4 seedy: **brak (1 seed)** → **odrzucona**
- `Dplus.r128.HDLZ.sbmc.d5.s1`: K1 zysk ≥ 10 % w ≥ 2 binach (MC i transfer): **tak** (3/3); K2 rzeźbienie nie gorsze niż H: **tak**; K3 SF: **brak pomiaru**; K4 seedy: **brak (1 seed)** → **warunkowo**

<!-- WYNIKI:END -->

## 11. Ocena etapu 1 (17.11.2026 17:45; tabele w §10)
G = Zproxy/Zproxy(XGB H, tło z SB danych) przy WP70, biny pT 3–5 / 5–8 / 8–30 GeV; „MC” = tło MC truth ≥ 8, „dane” = ε_B z dalekich SB 428770.

| konfiguracja | G na MC | G na danych | R_sculpt MC | werdykt |
|---|---|---|---|---|
| D⁰ sieć r128, H + kinematyka, DisCo λ = 5 (3 ziarna) | 1,159 ± 0,002 / 1,167 ± 0,004 / 1,142 ± 0,008 | 1,03 / 1,05 / 1,04 | 0,99 | odrzucona (K1: zysk nie przenosi się na dane) |
| D⁰ sieć r128, te same 18 cech co H, λ = 5 | 1,04 / 1,08 / 1,07 | 0,96 / 0,99 / 0,99 | 0,95 | = drzewo H trenowane na MC (1,03 / 1,06 / 1,07; 0,97 / 0,98 / 0,93) |
| D⁰ sieć, H + kinematyka, bez dekorelacji | 1,25 / 1,22 / 1,12 | 4,1 / 2,4 / 1,4 | 10,9 | odrzucona (rzeźbi masę) |
| D⁰ sieć, λ = 1 / 25 / 100 | 1,19 / 0,98 / 0,56 (3–5) | 1,88 / 0,88 / 0,51 | 2,60 / 0,99 / 1,00 | λ = 1 rzeźbi, λ ≥ 25 zabija zysk |
| D⁰ drzewo, H + kinematyka, tło MC | 1,15 / 1,16 / 1,11 | 1,11 / 1,08 / 0,96 | 1,07 | rzeźbi (R 1,07 ± 0,01 wobec 0,94) |
| D⁰ drzewo, H + cechy dodatkowe (HX), tło MC | 1,14 / 1,14 / 1,14 | 1,04 / 1,05 / 0,99 | 0,97 | zysk < 10 % na danych |
| D⁰ sieć z izolacją (HKia / HKir), λ = 5 | 1,34 / 1,35 / 1,34 | 1,61 / 1,59 / 1,40 | 1,00 | warunkowo (1 ziarno; SF niżej) |
| D⁺ sieć r128, H+DLZ, λ = 5 (3 ziarna) | 1,256 ± 0,010 / 1,201 ± 0,005 / 1,107 ± 0,016 | 1,24 / 1,16 / 1,04 | 0,92–0,94 | warunkowo (brak pomiaru SF dla D⁺) |
| D⁺ drzewo H+DLZ, tło MC | 1,27 / 1,18 / 1,12 | 1,39 / 1,27 / 0,99 | 0,98 | odniesienie |

SF z tagu D* (D⁰, WP70, 5–8 / 8–30 GeV): sieć bez izolacji 1,06 ± 0,13 / 0,93 ± 0,08 (427929), 0,88 ± 0,12 / 1,00 ± 0,04 (428770 all) — jak H (1,01 ± 0,13 / 0,97 ± 0,10; 0,88 ± 0,12 / 0,97 ± 0,04). Sieć z izolacją: 1,04 ± 0,13 / 0,92 ± 0,09 (427929), 0,81 ± 0,11 / 0,94 ± 0,04 (428770 all), przy WP50 w 8–30 GeV 0,85 ± 0,04. Bin 3–5 GeV niemierzalny (błędy 40–200 %).

Wnioski:
1. Technologia nie ma znaczenia: przy tych samych cechach sieć i drzewo dają to samo (różnice ≤ 2 %). Decydują cechy i źródło tła treningowego.
2. Dla D⁰ bez izolacji żadna konfiguracja nie spełnia K1: zysk +14–17 % na MC kurczy się do +3–5 % na danych. Drzewo H zostaje.
3. Dekorelacja masy jest konieczna, gdy na wejściu jest kinematyka rozpadu: bez niej R_sculpt = 10,9. DisCo λ = 5 sprowadza R do 0,99 i corr(score, m) do 0,02 przy stracie ok. 9 pp zysku na MC (1,25 → 1,16 w 3–5 GeV); λ = 25 kasuje cały zysk.
4. Izolacja to jedyne wejście dające duży zysk (+34 % MC, +40–60 % na SB danych), ale kosztuje SF w runie jetowym (0,81 ± 0,11 w 5–8 GeV, 0,85 ± 0,04 przy WP50 w 8–30 GeV); w 427929 SF zgodny z 1. To potwierdza wcześniejszy wynik o izolacji, teraz z siecią zdekorelowaną.
5. D⁺: sieć z DisCo dorównuje drzewu H+DLZ na MC, a ma corr(score, m) 0,02 wobec −0,22 i R na danych 1,02 wobec 1,10 dla drzewa z tłem z SB. To kandydat na obejście niezdanego kryterium rzeźbienia H+DLZ (DLZ_dplus/RAPORT.md §9), ale kalibracja ε_S dla D⁺ pozostaje nierozwiązana.
6. Warunek etapu 2 (§4): (a) niespełniony, (b) spełniony co do zysku (G na danych 1,6 wobec 1,03), R_SF = 0,92–0,94 w 428770 i 0,98–0,99 w 427929 (próg 0,95). Etap 2 liczy się (przyjęto 17.11).

## 12. Ocena etapu 2 (18.11.2026 10:50; tabele w §10, wiersze `D0.ds.*` i `D0.tr.*`, liczone na wierszach wspólnych z odniesieniem — „podzbiór”)

Wejścia: 18 cech globalnych H + zbiór torów (córki; córki + otoczenie ΔR < 0,4 = `e04`, 14 tokenów; ΔR < 1,0 = `e10`, 22 tokeny), 11 cech na tor. DeepSets (`ds`) i mały transformer (`tr`), DisCo λ = 5 (`d5`), WP70, biny 3–5 / 5–8 / 8–30 GeV.

| konfiguracja | G na MC | G na danych (SB 428770) | R_sculpt MC | SF WP70 427929 (5–8 / 8–30) | SF WP70 428770 all (5–8 / 8–30) |
|---|---|---|---|---|---|
| XGB H (odniesienie) | 1 | 1 | 0,942 ± 0,006 | 1,01 ± 0,13 / 0,97 ± 0,10 | 0,88 ± 0,12 / 0,97 ± 0,04 |
| `tr` same córki (3 ziarna) | 1,042 ± 0,007 / 1,061 ± 0,008 / 1,050 ± 0,015 | 0,94 / 0,98 / 0,97 | 0,97 | 1,04 ± 0,13 / 0,93 ± 0,08 | 0,89 ± 0,12 / 0,99 ± 0,04 |
| `ds` + e04 | 1,20 / 1,26 / 1,26 | 1,29 / 1,47 / 1,41 | 0,974 ± 0,012 | 1,04 ± 0,13 / 0,95 ± 0,09 | 0,81 ± 0,11 / 0,93 ± 0,04 |
| `tr` + e04 (3 ziarna) | 1,450 ± 0,024 / 1,615 ± 0,032 / 1,546 ± 0,008 | 1,54–1,58 / 1,72–1,79 / 1,60–1,63 | 0,975 ± 0,015 | 1,02 ± 0,11 / 0,97 ± 0,06 | 0,84 ± 0,11 / 0,90 ± 0,04 |
| `tr` + e04 bez DisCo | 1,58 / 1,71 / 1,59 | 1,69 / 1,79 / 1,55 | 0,969 (corr −0,053, nachylenie −0,25) | 1,08 ± 0,12 / 0,92 ± 0,06 | 0,83 ± 0,11 / 0,93 ± 0,04 |
| `tr` + e04, tło SB+MC | 1,37 / 1,45 / 1,33 | 1,65 / 1,85 / 1,77 | 0,971 | 1,03 ± 0,12 / 0,93 ± 0,07 | 0,76 ± 0,10 / 0,79 ± 0,04 |
| `tr` + e10 | 1,62 / 1,71 / 1,66 | 1,94 / 1,89 / 1,75 | 0,971 ± 0,016 | 0,99 ± 0,11 / 0,95 ± 0,07 | 0,76 ± 0,10 / 0,89 ± 0,04 |

Wg kryteriów §4 dla `tr` + e04, λ = 5:
- **K1 (zysk): spełnione** z dużym zapasem na MC i na danych (+45–60 % na MC, +54–79 % w ε_S/√ε_B na pasmach bocznych danych).
- **K2 (rzeźbienie): spełnione** — R_sculpt 0,975 ± 0,015 (H: 0,942), R na paśmie danych 1,00–1,01 ± 0,04, corr(score, m) 0,003–0,005, nachylenie −0,10…−0,13 (H: −0,30).
- **K4 (stabilność): spełnione** — rozrzut między 3 ziarnami 0,024 / 0,032 przy zysku 0,45 / 0,62.
- **K3 (kalibracja): niespełnione w 428770, spełnione w 427929.** 428770 all 8–30 GeV: 0,90 ± 0,04 (H 0,97 ± 0,04), J12only 0,90 ± 0,04 (H 0,99 ± 0,05); przy WP50 0,83 ± 0,04. C2ST_1D sygnału 0,56–0,59 (H 0,51–0,56). W 427929 SF 1,02 ± 0,11 / 0,97 ± 0,06 — zgodne z 1 i z H.
- **Werdykt: nieprzyjęta jako nominał** (K3). W czystym minimum-bias (427929) spełnia wszystkie cztery kryteria w granicach błędu 6–11 % — kandydat warunkowy tam, gdzie próbka jest MBTS.

Wnioski:
1. Cały zysk pochodzi z torów otoczenia, nie z architektury ani z parametrów córek: transformer na samych córkach daje +4–6 % na MC i −2…−6 % na danych.
2. Otoczenie z torów bije dwa stożki izolacji z etapu 1 (`r128.HKia`: MC 1,34, dane 1,61 / 1,59 / 1,40) i transformer bije DeepSets; szersze otoczenie (ΔR < 1,0) daje jeszcze więcej na danych, ale gorszy SF w 428770 (0,76 ± 0,10 w 5–8 GeV).
3. Koszt jest ten sam co dla izolacji: efektywność sygnału w próbce jetowej 7–10 % (WP70) do 17 % (WP50) poniżej MC. Zysk netto pozostaje duży, ale wymaga SF per klasa triggera.
4. Trening z tłem z pasm bocznych danych (`sbmc`) uczy się różnicy dane/MC: najlepsze G na danych, najgorszy SF (0,79 ± 0,04; WP50 0,65 ± 0,03) → odrzucone.
5. DisCo λ = 5 kosztuje 5–8 % zysku na MC i usuwa resztkową korelację z masą; zostaje.
Zastrzeżenia: SF w 3–5 GeV nieokreślony; e10, DeepSets i warianty bez DisCo/sbmc na jednym ziarnie; G na danych z pasm bocznych 428770 (jetowy); 5 bloków liczone na CPU po awarii GPU (V100), wyniki nie zależą od urządzenia.

## 13. Ocena etapu 3 i zamknięcie linii sieci (18.11.2026 23:45; tabele: `runs/lxplus3/REPORT_etap3.md`, liczby: `runs/lxplus3/dnet3_results.json`)

Kryteria: `KRYTERIA_etap3.md` (zapisane przed wysłaniem). Jeden odczyt, 1 ziarno, WP70; próg E1 (|ΔSF| ≤ 0,04 przy błędzie SF ±0,04) ma małą moc — wynik jest wskazówką, nie rozstrzygnięciem.

| konfiguracja | SF 428770 all 8–30 | SF 428770 J12only 8–30 | SF 427929 5–8 / 8–30 | G pasma 428770 (3–5 / 5–8 / 8–30) | G pasma 427929 | G MC | R_sculpt MC |
|---|---|---|---|---|---|---|---|
| drzewo H | 0,968 ± 0,043 | 0,994 ± 0,050 | 1,005 ± 0,126 / 0,966 ± 0,097 | 1 / 1 / 1 | 1 / 1 / 1 | 1 / 1 / 1 | 0,942 ± 0,006 |
| `tr.H+e04` (etap 2) | 0,895 ± 0,038 | 0,898 ± 0,044 | 1,018 ± 0,115 / 0,971 ± 0,062 | 1,58 / 1,79 / 1,60 | 1,31 / 1,37 / 1,44 | 1,47 / 1,65 / 1,55 | 0,975 ± 0,015 |
| `e04g` (otoczenie bez pędów) | 0,940 ± 0,042 | 0,962 ± 0,049 | 1,093 ± 0,136 / 0,915 ± 0,091 | 1,14 / 1,26 / 1,28 | 1,02 / 1,01 / 1,17 | 1,13 / 1,17 / 1,22 | 0,965 ± 0,011 |
| `e04pt1` (tory > 1 GeV) | 0,953 ± 0,043 | 0,978 ± 0,050 | 1,067 ± 0,132 / 0,894 ± 0,086 | 1,21 / 1,35 / 1,40 | 1,01 / 1,04 / 1,19 | 1,15 / 1,21 / 1,25 | 0,966 ± 0,011 |
| `e04c` (stożek z 20 torów) | 0,912 ± 0,039 | 0,907 ± 0,044 | 1,111 ± 0,125 / 0,964 ± 0,061 | 1,57 / 1,73 / 1,66 | 1,32 / 1,29 / 1,51 | 1,49 / 1,66 / 1,61 | 0,972 ± 0,015 |
| `aS2` (dekor. aktywności, sygnał, λ 2) | 0,927 ± 0,039 | 0,924 ± 0,045 | 0,992 ± 0,115 / 0,945 ± 0,064 | 1,41 / 1,58 / 1,53 | 1,23 / 1,31 / 1,36 | 1,36 / 1,55 / 1,54 | 0,967 ± 0,013 |
| `aS10` (λ 10) | 0,912 ± 0,039 | 0,916 ± 0,045 | 0,992 ± 0,116 / 0,914 ± 0,076 | 1,21 / 1,49 / 1,38 | 1,10 / 1,26 / 1,26 | 1,20 / 1,46 / 1,46 | 0,967 ± 0,012 |
| `aB10` (dekor. aktywności, tło, λ 10) | 0,956 ± 0,042 | 0,974 ± 0,049 | 0,967 ± 0,114 / 0,914 ± 0,075 | 0,86 / 1,05 / 0,96 | 0,88 / 1,02 / 0,99 | 0,92 / 1,03 / 1,01 | 0,973 ± 0,009 |

Wszystkie fity SF określone. dCor²(score, aktywność) na sygnale MC (n torów; Σ pT): bez członu 0,005–0,010; 0,006–0,018 → `aS2` 0,0018; 0,0055 → `aS10` 0,0008; 0,0021 (dekorelacja zadziałała); na tle: bez członu 0,04–0,10 → `aB10` 0,0012; 0,0006.

**Werdykt wg kryteriów.**
- E1 (kalibracja w 428770): spełniają `e04g`, `e04pt1`, `aB10`; nie spełniają `e04c`, `aS2` (all 0,041; J12only 0,070), `aS10`.
- E2 (G ≥ 1,20 w ≥ 2 binach na pasmach 428770): `e04g` (2 z 3), `e04pt1` (3 z 3) spełniają; `aB10` nie (G ≈ 1).
- E3, E4: `e04g` i `e04pt1` spełniają.
- Formalnie `e04g` i `e04pt1` spełniają E1–E4 → wg zapisu byliby „kandydatami do potwierdzenia na 3 ziarnach”.

**Odczyt fizyczny (to jest wynik etapu).**
1. Odzyskanie kalibracji kosztuje większość zysku: `e04g`/`e04pt1` mają G 1,14–1,40 na pasmach 428770 wobec 1,58–1,79 pełnej sieci. Zysk i rozbieżność dane/MC siedzą w tej samej informacji: miękkich torach otoczenia i ich pędach.
2. `aB10`: po uniezależnieniu score tła od aktywności zysk znika (G ≈ 1) przy poprawnym SF → **zysk etapu 2 to rozdzielanie tła po aktywności otoczenia**.
3. `aS2`/`aS10`: dekorelacja score sygnału od aktywności działa (dCor² spada 5–10×), ale SF wraca tylko częściowo (0,91–0,93 wobec 0,895) → rozbieżność nie jest wyłącznie zależnością od liczby i sumy pT torów.
4. `e04c`: obcięcie do 12 torów nie było przyczyną (SF 0,91).
5. **G na pasmach czystego runu 427929 jest niższe niż na pasmach 428770**: pełna sieć 1,31 / 1,37 / 1,44 (nie 1,58 / 1,79 / 1,60), a warianty „skalibrowane” 1,01–1,04 poniżej 8 GeV i 1,17–1,19 w 8–30 GeV. Liczby „G na danych” z §12 (pasma 428770, na których trenowano drzewo odniesienia i które są próbką dżetową — skim, `SKIM_review_2026-11-18`) zawyżają zysk; właściwą liczbą dla próbki minimum-bias jest G z 427929.

**Decyzja (18.11).** Potwierdzenia na 3 ziarnach NIE robimy i linię sieci zamykamy:
- w 428770 warianty spełniające E1 dają na danych 1,14–1,40, ale ta próbka jest dżetowym skimem bez klasy minimum-bias, a różnica SF 0,94–0,95 wobec 0,895 to ~1σ na jednym ziarnie;
- w 427929, gdzie sieć miałaby sens, warianty „skalibrowane” nie dają zysku poniżej 8 GeV, a pełna `tr.H+e04` daje 1,31–1,44 przy SF 1,02 ± 0,12 / 0,97 ± 0,06, czyli nie wymaga poprawek.
Stan końcowy: nominałem zostaje drzewo H; `tr.H+e04` = kandydat warunkowy wyłącznie dla próbki minimum-bias (427929) z zyskiem na danych 1,3–1,4 (nie 1,5–1,8). Warunek powrotu: decyzja o przejściu analizy na 427929 i ≥ 3× większa statystyka D* do kalibracji.

## 14. Kierunek od 19.11 i ocena etapu 5a (19.11.2026 19:15)

**Uchylenie decyzji z §13.** Decyzja z 19.11: rozwijamy własną sieć, zgodność z cudzą siecią nie jest celem. Kryterium „SF sieci = SF drzewa” porzucone; SF ≠ 1 dopuszczalny, jeśli zmierzony i stabilny. Etap 4 (`KRYTERIA_etap4.md`) i 5a (`KRYTERIA_etap5.md`, plan `PLAN_etap5.md`) wysłane 19.11.

**Etap 5a — odczyt jednorazowy** (`runs/lxplus5/WERDYKT_etap5.md`, skrypt `verdict5.py`; komplet 7/7 wariantów, `fail/` puste; treningi na CPU po foldach, GPU tylko sonda).

| pozycja | wynik (3–5 / 5–8 / 8–30 GeV) | próg | werdykt |
|---|---|---|---|
| cechy par w uwadze (`trp`), Z/Z(tr) na MC | 1,007 ± 0,003 / 1,002 ± 0,012 / 1,030 ± 0,032; pasma 427929: 1,007 / 0,978 / 0,955 | ≥ 1,05 i 3σ; ≥ 1,03 | porzucić |
| krzywa uczenia, Z(100 %)/Z(50 %) na MC | 1,023 ± 0,003 / 1,028 ± 0,011 / 1,014 ± 0,024; Z(50)/Z(25): 1,024 / 1,006 / 0,979 | zamknięte < 0,02; otwiera ≥ 0,03 | nierozstrzygnięte wg zapisu |
| tagger otoczenia E po cięciu H, D⁰: G_E na pasmach 427929 (G_Htight) | 1,22 (2,03) / 1,38 (1,51) / 1,31 (1,00) | G_E ≥ 1,15 i G_E − G_Htight ≥ 0,10 w ≥ 2 binach | niespełnione (1/3) |
| tagger E, D⁺ | 1,12 (1,50) / 1,39 (1,35) / 2,04 (1,23) | jw. | niespełnione (1/3) |
| sonda GPU | kolejka 1047 s; 2,6 s/epokę vs 195 s CPU (×75) | ≤ 2 h i ≥ 10× | spełnione |

Decyzje: (i) `trp` zamknięte; (ii) krzywa uczenia: dodatkowych ziaren wariantu 50 % NIE zlecam (odstępstwo od zapisu §2, uzasadnienie: 1–3 % na podwojenie nie uzasadnia pretreningu ani wspólnego treningu niezależnie od rozstrzygnięcia) — pozycje #8 planu (pretrening, wspólny koder) zamknięte jako nieopłacalne; (iii) kaskada H × E zamknięta dla obu gatunków, pomiar SF_E (5b) nieuruchamiany.

**Obserwacja post hoc (poza kryteriami):** G_Htight = 2,03 (3–5 GeV) i 1,51 (5–8 GeV) na pasmach 427929 oznacza, że samo zaostrzenie drzewa H z WP70 do wydajności ≈ 0,35–0,42 podnosi Z/Z_H(WP70) bardziej niż sieć przy WP70. Dotychczasowa miara G (stała ε = 0,7) nie porównuje klasyfikatorów w ich optymalnych punktach pracy; przy niskim pT optimum leży przy ε ≪ 0,7 (zgodne z 2–3 GeV: optimum ε ≈ 0,2). Do odczytu etapu 4 dopisać porównanie max_ε Z(sieć) / max_ε Z(H) jako post hoc.

Lekcje infrastruktury: GPU 60–75× szybsze na epokę, kolejka kart nieprzewidywalna (17 min rano, > 3 h potem); priorytet konta wyczerpany ~400 jobami `sdvar` → etapy wysyłać wg ważności; `dn4_train` inicjalizuje fold 0 przed ustawieniem ziarna (niepowtarzalny bitowo).

## 15. Ocena etapu 4 (19.11.2026 20:20; `runs/lxplus4/WERDYKT_etap4.md`, skrypt `verdict4.py`; odczyt jednorazowy wg `KRYTERIA_etap4.md`)

Komplet: 9/9 modeli, app 100 % plików obu runów i MC, `fail/` puste. Brak: `dnet4_d23.json` (węzeł `ana23` jeszcze w kolejce o 20:07) — bin 2–3 GeV nierozstrzygnięty.

**Werdykt: sieć NIE zostaje nominałem; drzewo H zostaje.** Trzy niezależne powody:

1. **Wolny pion D\* (P1, MC, WP70):** wydajność sieci dla D⁰ z D\* ε_a = 0,868 / 0,890 / 0,887 (3–5 / 5–8 / 8–30 GeV), dla pozostałych D⁰ ε_b = 0,630 / 0,608 / 0,604, a dla D⁰←D\* z zamaskowanym pionem ε_c = 0,617 / 0,600 / 0,597. Podwójny stosunek sieć/drzewo 1,37 / 1,45 / 1,44 (drzewo: 1,00–1,02). Sieć rozpoznaje wolny pion w otoczeniu. Skutek: wszystkie dotychczasowe SF sieci z tagu D\* (§12–§13: 0,90 w 428770, 1,02 / 0,97 w 427929) kalibrują tylko D⁰ z D\*, nie próbkę inkluzywną; część „zysku” sieci na sygnale to preferowanie D⁰←D\*. Dotyczy wszystkich ziaren i wariantu sb2mc.
2. **Metoda R dla D⁰ nie osiąga precyzji już na MC:** R_fit − R_prawda = +0,16 ± 0,15 / +0,01 ± 0,06 / −0,03 ± 0,07; nieobciążona (średni pull 21 domknięć +0,04), ale σ > 0,05 we wszystkich binach → bramka zablokowała dane (192 punkty). Dodatkowo R_prawda ≠ 1 na MC: 0,957 / 0,932 / 0,905 (przenoszenie ε między próbkami MC). D⁰ nie ma więc żadnej kalibracji sieci na próbce inkluzywnej.
3. **D⁺ (domknięcie przechodzi w 5–8 i 8–30 GeV, tam dane otwarte):** R = 0,933 ± 0,043 / 0,940 ± 0,038 w 427929 i 0,874 ± 0,036 / 0,892 ± 0,027 w 428770 (ziarno 1; rozstęp ziaren ≤ 0,04 przy 8–30 GeV). **Zysk rzeczywisty z fitu piku G_fit = 1,11 / 1,03 w 427929 (próg 1,15 niespełniony) i 1,24 / 1,17 w 428770 (spełniony)**, podczas gdy proxy z pasm dawało 1,36 / 1,26 i więcej. W 428770 R zależy od klasy triggera przy 8–30 GeV: J12only 0,868 ± 0,019 vs reszta 1,030 ± 0,045 (3,3σ) → SF per klasa. Bin 3–5 GeV zablokowany (σ domknięcia 0,054).

Wnioski: (i) proxy Z na pasmach zawyża zysk — w czystym minimum-bias realny zysk w błędzie uzysku D⁺ to 3–11 %; (ii) zysk realny jest w próbce dżetowej (17–24 %), tam gdzie SF zależy od triggera; (iii) razem z §14 (kaskada, cechy par, krzywa uczenia) linia „sieć zamiast drzewa” zamknięta wynikiem negatywnym. Jedyny nietestowany wariant: sieć z usuniętym z wejścia torem zgodnym z wolnym pionem D\* (wtedy tag D\* kalibrowałby próbkę inkluzywną) — do decyzji, oczekiwany zysk w 427929 nie większy niż dla D⁺.
