# charmbkg — przegląd wyników i stanu (20.11.2026)

Projekt B. Pytanie wyjściowe: czy tło pod pikiem masy mezonu D da się odtworzyć **bez zaglądania do sygnału** — z obiektów różniących się od
sygnału jednym składnikiem (masą, hipotezą K↔π, pochodzeniem MC, torem) — i czy takie tło jako tło treningowe klasyfikatora podnosi
istotność. Od 17.11.2026 projekt obejmuje też kalibrację dane/MC, pomiar inkluzywny D⁰, D⁺, D* w runie 427929 i region z tagiem protonu
w AFP. Pakiet Pythona `charmbkg` + skrypty analizy (`analysis/`) + raporty z liczbami (`docs/results/`).


Osobny projekt od **double charm** (pary D–D, repo [charmtag](https://github.com/333963/charmtag), pakiet `charmtag`). charmtag korzysta
z `charmbkg` w przypiętej wersji (SALAD jako systematyka kształtu tła, HYP, źródła tła treningowego); charmbkg nie zależy od charmtag.
Wydzielony z `charmtag` 2.0.3 15.11.2026.

Kod bazowy analizy: `charmpy` z CERN GitLab (`berkeleylab/CharmPhysics/charmpy`, gałąź minimum-bias, **tylko odczyt**).
Własne poprawki do charmpy (k-fold, DisCo, tło treningowe z sidebandów, szablon refleksji, weta D⁺) są w prywatnym forku
[333963/charmpy](https://github.com/333963/charmpy) (gałęzie `e_fix` → `integration`) — to łatki na kod GitLaba, nie część tego pakietu.

## Stan (20.11.2026)

1. **Hipoteza wyjściowa — „lepsze tło treningowe → większe Z” — nie potwierdziła się** (stan z 15.11, tabela w sekcji
   [Tło treningowe](#tło-treningowe-stan-z-15112026)). Kolejne próby (sieć, nowe cechy, tło i sygnał „przez konstrukcję”) dały ten sam
   wniosek: wynik ogranicza statystyka, nie metoda. Drzewo XGB „H” przy WP70 zostaje nominałem.
2. **Pomiar inkluzywny** σ(D⁰), σ(D⁺), σ(D*±) w runie 427929 (MBTS, 0,4648 nb⁻¹) — policzony, wstępny, w większości post hoc.
3. **Region z tagiem protonu (AFP)**: tło przypadkowe z mieszania protonów, regiony trygera i świetlność na łańcuch, pierwszy odczyt
   D⁰ + proton (wskazanie na granicy 3σ, nie obserwacja), ułamek fiducjalny R^fid.

Indeks wszystkich raportów (jedno zdanie na temat): [docs/results/README.md](docs/results/README.md). Szkic noty metodycznej:
[docs/METHODS_NOTE_DRAFT.md](docs/METHODS_NOTE_DRAFT.md); dodatek o regionie z tagiem: [docs/TAGGED_ADDENDUM.md](docs/TAGGED_ADDENDUM.md).

### Wyniki z liczbami

Run 427929, WP70 drzewa H, D + anty-D, prompt + nieprompt, |η(D)| < 2,2; σ [µb] ± stat ± syst ± świetlność.
Źródło: [docs/results/xsec_final/v2/tabele/final.md](docs/results/xsec_final/v2/tabele/final.md) (D*: [xsec_atlas](docs/results/xsec_atlas/RAPORT.md)).

| | 3–5 GeV | 5–8 GeV | 8–30 GeV |
|---|---|---|---|
| σ(D⁰) | 1028 ± 119 ± 274 ± 23 | 278 ± 29 ± 49 ± 6 | 100 ± 10 ± 21 ± 2 |
| σ(D⁺) | 515 ± 44 ± 82 ± 11 | 165 ± 10 ± 25 ± 4 | 49,2 ± 4,0 ± 6,9 ± 1,1 |
| σ(D*±) | – | 134 ± 18 ± 23 ± 3 | 60,6 ± 7,1 ± 6,8 ± 1,3 |
| D⁺/D⁰ | 0,50 ± 0,07 ± 0,15 | 0,59 ± 0,07 ± 0,12 | 0,49 ± 0,06 ± 0,12 |
| uzysk D⁰ (model M8) | 4059 ± 468 | 1536 ± 161 | 700 ± 70 |
| uzysk D⁺ (model A2) | 3404 ± 292 | 1630 ± 102 | 651 ± 53 |
| C = A·ε (D⁰ / D⁺) | 22,3 / 16,1 % | 30,0 / 23,1 % | 39,6 / 30,8 % |
| SF klasyfikatora (D⁰ / D⁺) | 0,967 ± 0,033 / 0,949 ± 0,080 | 1,005 ± 0,135 / 0,993 ± 0,052 | 0,965 ± 0,109 / 0,990 ± 0,051 |

Pozostałe liczby (każda z plikiem źródłowym w repo):

| wielkość | wartość | status | źródło |
|---|---|---|---|
| świetlność, 427929 (nasza kopia plików) | 0,4648 ± 0,0103 nb⁻¹ (cały GRL: 0,4737) | skala bezwzględna to założenie 2,2 %; sprawdzian MBTS niezgodny (niżej) | [lumi](docs/results/lumi/RAPORT.md), `lumi.json` |
| świetlność, 428770, regiony trygera | R1 0,4687 ± 0,0104; R2 27,65 ± 0,61; R3 29,41 ± 0,68 nb⁻¹ | żaden region nie daje normalizacji bezwzględnej dla zdarzeń z D | [lumi](docs/results/lumi/RAPORT.md), [tagreg v2](docs/results/tagreg/v2/RAPORT.md), [reg428](docs/results/reg428/RAPORT.md) |
| ułamek tagów przypadkowych f_acc | 0,0087 ± 0,0044 (427929); 0,0477 ± 0,0135 (428770) | wstępny; tylko trygery nie-AFP (w 428770 to 4,4 % tagów) | [afp_mixing.md](docs/afp_mixing.md) §7, [afp_mix](docs/results/afp_mix/runs/lxplus/out_v2/RAPORT_DANE_v2.md) |
| odczyt z tagiem, D⁰, 427929, 3–30 GeV | N_tag = 54,7 −18,7 +19,4 (stat) ± 12,6 (kształt); Z = 3,05 lokalnie (zabawki p = 0,0015) | **wskazanie, nie obserwacja**; kryteria zapisane przed odczytem; wariant „okno − pasma” ≈ 30 | [tagged](docs/results/tagged/RAPORT.md) |
| odczyt z tagiem, D⁺, 427929 | N_tag = 29,1 −10,8 +11,5 ± 8,8 (Z = 2,9) | opisowo: nie spełnia 2 kryteriów kontrolnych | [tagged](docs/results/tagged/RAPORT.md) |
| odczyt z tagiem, 428770, region R2 (dżet + AFP) | D⁰ 278 ± 49 ± 48; D⁺ 185 ± 30 ± 24 | próbka silnie obciążona trygerem, bez normalizacji | [tagged](docs/results/tagged/RAPORT.md) §0 |
| R^fid = σ_tag^fid / σ_incl (D⁰) | 0,99 −0,34 +0,35 (stat) ± 0,29 (syst) %; wariant 0,54 ± 0,23 ± 0,09 % | **post hoc** (N_tag znane przed rachunkiem) | [TAGGED_ADDENDUM.md](docs/TAGGED_ADDENDUM.md) §4, [tagxsec](docs/results/tagxsec/RAPORT.md) |
| σ_tag^fid (D⁰) | 13,9 −4,8 +5,0 ± 4,6 µb (dodatek, wejścia końcowe); 13,6 −4,7 +4,8 ± 4,6 µb (raport `tagxsec`, wejścia wcześniejsze) | post hoc; dwie wersje — patrz „Znane problemy” | jw. |
| A·ε_AFP | 0,92 ± 0,05 (wcześniej zakładane 0,55) | akceptancja z szybkiej symulacji, ε z danych | afp_acceptance |
| ρ = C(tag)/C(ND) | 0,96 ± 0,12 | z bardzo małej próbki MC SD | [TAGGED_ADDENDUM.md](docs/TAGGED_ADDENDUM.md) §4 |
| Pythia 8, poziom cząstkowy: R(D⁰) | 0,384 ± 0,005 % (oficjalne A3), obwiednia modeli 0,154–1,53 % | przewidywanie, tylko błąd stat. MC | [sd_var](docs/results/sd_var/RAPORT.md), wcześniej [sd_pred](docs/results/sd_pred/RAPORT.md) |
| Pythia 8: przewidywany ułamek tagowany D⁰ | 0,322 % (0,129–1,28 %) przy A·ε_AFP = 0,92; R^fid 0,35 % (0,14–1,40 %) w dodatku | jw. | [sd_var, domknięcie](docs/results/sd_var/RAPORT_domkniecie.md), [TAGGED_ADDENDUM.md](docs/TAGGED_ADDENDUM.md) |
| SF izolacji D⁰ (ptcone40 < 2 GeV, tag D*) | 0,98 ± 0,04 (427929) wobec 0,51 ± 0,02 (428770) | deficyt idzie za trygerem dżetowym, nie za runem ani pile-upem | [b_trig](docs/results/b_trig/runs/lxplus/out/RAPORT_DANE.md), [ISO_origin](docs/results/b_trig/ISO_origin/RAPORT.md) |
| weto wierzchołków pile-up, ε dla sygnału D* | 93,95 ± 0,73 % (428770); 98,60 ± 0,69 % (427929); MC 99,25 ± 0,04 % | nie znosi się w stosunku tagowane/inkluzywne | [veto_pv](docs/results/veto_pv/RAPORT.md) |
| SF klasyfikatora dane/MC z tagu D* (D⁰, 5–30 GeV), SF_ref | 0,968 ± 0,029 (stat) ± 0,012 (model fitu) | post hoc; przenoszony na D⁺ i na D⁰ 3–5 GeV (warunkowo) | [sfabs](docs/results/sfabs/RAPORT.md), `akapit_do_noty.md` |
| współczynnik dane/MC arbitra hipotez K↔π (HYP; to nie jest SF klasyfikatora) | 0,977 ± 0,041 | zamknięte | [docs/hyp.md](docs/hyp.md) |
| pozycja i szerokość piku, dane − MC | −0,43 ± 0,51 MeV; stosunek szerokości 0,956 ± 0,034 | rozmiar systematyki, nie poprawka | [peak_pos](docs/results/peak_pos/RAPORT.md) |
| sieć (DNET, etap 6): zysk realny G_fit | 0,99 / 1,10 / 1,09 (427929); 1,13 / 1,25 / 1,17 (428770) | bez błędu, jedno ziarno; werdykt: SF nieokreślony, H zostaje | [dnet](docs/results/dnet/RAPORT.md) |

### Wyniki negatywne

Podejścia sprawdzone i odrzucone — każde z raportem (szczegóły i liczby w [indeksie](docs/results/README.md)):

- **Klasyfikator:** sieć zamiast drzewa (`dnet`), D⁺ z cechami Dalitza i kinematyką córek (`dlz_dplus`), nowe cechy D⁰ (`dlz_d0`),
  izolacja i dE/dx jako wejścia (`feat_charmnet`, `cn_clone`), sygnał z danych zamiast z MC (`sig_constr`).
- **Punkt pracy i fit:** ciasny punkt pracy z SF_rel z tego samego piku (`wptight`), fit w kategoriach score na danych (`catfit`),
  SF z kształtu score (`fit_scorebins/R3`), kalibracja ε_S skanem punktu pracy (`lowpt_calib`).
- **Tło:** pomiar ułamka MisMatch* z danych (`fmm_data`), szablon tła Δm z kombinacji D⁰ + tor poza 427929 8–30 GeV (`dstws`),
  przeważanie w zmiennych dżetów torowych i odjęcie piedestału jako wyjaśnienie deficytu izolacji (`b_trig/TRKJET`, `TRKJET2`, `iso_ped`);
  wcześniejsze: candrot, pseudo-WS, mixing toru ([docs/signal_like_background.md](docs/signal_like_background.md)).
- **Normalizacja:** żaden region runu 428770 nie daje drugiego punktu bezwzględnego (`reg428`).

### Co jest wstępne, a co post hoc

- **Wstępne jest wszystko**: nic nie przeszło recenzji wewnętrznej ATLAS, liczby pochodzą z własnej kopii ntupli v2.0.
- **Zapisane przed obejrzeniem wyniku**: odczyt z tagiem (`tagged`), bramki mieszania protonów
  (`afp_mix`), etapy sieci (`dnet`), test H+DLZ na danych (`dlz_dplus/DATA_2026-11-17`), `b_trig`, `iso_ped`, `lowpt_calib`.
- **Post hoc** (po poznaniu wcześniejszych liczb; raporty to oznaczają): model fitu D⁰ M8 i D⁺ A2, końcowy przekrój (`xsec_final` v2),
  SF bezwzględny (`sfabs`), skan punktu pracy (`wpscan`), skład tła (`bkgcat`), R^fid i σ_tag^fid (`tagxsec`), sprawdzian z tłem
  konstrukcyjnym (`konstr`), porównanie z Pythią. Wg oceny z 20.11 wszystko po 19.11 jest post hoc.

### Znane problemy

Z przeglądu kodu i stanu wyników (H2 §1.2, §3.2–3.4). Nic z tej listy nie zostało
naprawione w wynikach.

Wyniki:

1. **Model D⁰ (M8) jest wspólnym słabym punktem.** Domknięcie na MC przy 3–5 GeV: −22 ± 10 % (niepoprawione, niesione jako
   systematyka 24 %); dwa minima w skali szerokości przy 5–8 GeV. W ficie z tagiem z sygnałem związane są szerokie składowe
   (ok. 3,4 N_S), a rachunek „okno minus pasma” daje ok. 30 zamiast 55 zdarzeń.
2. **D⁺ (A2)**: kryterium domknięcia niespełnione w 5–8 GeV (+9,2 ± 6,2 %); żadne z 9 ustawień go nie spełnia.
3. **Świetlność**: skala bezwzględna niepotwierdzona; sprawdzian licznikiem MBTS dał σ_wid = 60,8 mb wobec oczekiwanych 69–72 mb
   (L zawyżona o 12–16 %?), wyjaśnienie tylko pośrednie.
4. **Stosunki wewnętrzne** D⁺/D⁰ i D*/D⁰ leżą powyżej pomiarów innych eksperymentów i sugerują σ(D⁰) niższe o 15–30 %;
   dane/Pythia i dane/CMS < 1 w większości binów.
5. **Dwie wersje σ_tag^fid**: 13,6 µb (`docs/results/tagxsec/`, `out/tagxsec_excl.json`) wobec 13,9 µb (dodatek, wejścia końcowe;
   R^fid syst 0,28 wobec 0,29). Rysunek T2 (`docs/results/presentation/`) czyta starszą wartość.
6. **Wycofana skala κ** (D⁰ × 1,364, D⁺ × 1,068; [cver](docs/results/cver/RAPORT.md) pkt 5): nadal jest wartością domyślną `--kscale` w
   `analysis/fitbias/skrypty/fb_var.py` i stałą w `analysis/konstr/skrypty/ka_data.py` — oznaczone w kodzie 20.11, nie zmienione,
   bo zapisane wyniki tych kroków z nią powstały ([analysis/fitbias/README.md](analysis/fitbias/README.md)). Uzyski przekroju końcowego
   są z fitów bez skali.
7. **f_acc**: werdykt „PASS” przy bramce G9 o statusie N/A; ekstrapolacja κ w ln ξ_cal traktowana jako systematyka (47 % błędu
   w 427929), nie poprawka; skok częstości poniżej ξ_cal ≈ 0,075 niewykrywany przez żadną bramkę; propozycja błędów asymetrycznych
   nieprzyjęta. Dotyczy tylko trygerów nie-AFP.
8. **ρ i akceptacja AFP** opierają się na małym MC dyfrakcyjnym (98 zrekonstruowanych D⁰, 25 z protonem w oknie).
9. Tory, wierzchołek 3-torowy i krotność — systematyki z literatury, niezmierzone. σ(D*±) bez własnego testu domknięcia w przeglądzie.
10. Ntuple runu 428770 pochodzą ze strumienia `physics_Main` (wybór zdarzeń, nie pełny strumień minimum-bias); łańcuchy AFP + tor
    i MBTS piszą do strumienia `MinBias` (potwierdzone w źródle menu trygera; istnienia AOD `physics_MinBias` nie sprawdzono):
    [docs/data_stream_428770.md](docs/data_stream_428770.md).
11. Przekroje w `sfabs/RAPORT.md` pochodzą sprzed korekty mianownika; obowiązuje `xsec_final/v2`.

Kod pakietu (niepoprawione; H2 §1.2):

- `charmbkg/hyp.py`, tryb `pertrack`: foldy treningu i predykcji liczone z różnych kluczy — gwarancja out-of-fold może nie działać
  (skutek dla wyników oceniony jako mały).
- `charmbkg/afp_mix.py`: `MixSum.total` / `MixSum.boot` po cichu ignorują `mask`; werdykt „PASS” mimo bramek N/A.
- `charmbkg/data/sf_rel_b2.json` + `sf_rel.py`: stan z 4.11, niezgodny z późniejszymi SF; moduł nieużywany. Nie brać SF z pakietu.
- `charmbkg/massfit.py` nie zawiera modeli M8 / A2 — silnik fitu wyników to `analysis/fitbias/skrypty/fb_fit.py` (druga kopia:
  `analysis/dnet/lxplus/dnet4/dn4_ana.py`).
- `charmbkg/bkggen/salad_v1.py`: kryterium `ok` spełnione z urzędu, gdy test pasma się nie wykona.

### Wynik → skrypt → raport

Skrypty w `analysis/` są kopiami z katalogu roboczego w stanie, w jakim liczyły; ścieżki wejść wskazują na lokalny obszar roboczy
i EOS, więc **z repo uruchamia się tylko to, co oznaczono „tak”**.

| wynik | skrypt | raport i liczby | z repo? |
|---|---|---|---|
| σ(D⁰), σ(D⁺), D⁺/D⁰ | `analysis/xsec_final/v2/xsec_final.py` ← `xsec_final/skrypty/kappa_fit.py` ← `bkgcat/skrypty/bc_m8.py` ← `fitbias/skrypty/{fb_fit,fb_var}.py`; D⁺: `fitbias_dplus/skrypty/fd_lib.py` | `docs/results/xsec_final/v2/tabele/final.md`, `out/final.json` | nie (histogramy poza repo) |
| σ(D*±) | `analysis/xsec_atlas/skrypty/xa_dstar.py` (+ `charmtag.dstar_refit`, kod projektu A) | `docs/results/xsec_atlas/RAPORT.md` | nie |
| C = A·ε | `analysis/xsec_atlas/skrypty/lxplus/xa_truth.py`, niezależnie `analysis/cver/skrypty/cver_*.py` | `docs/results/cver/RAPORT.md` | nie |
| SF klasyfikatora | tag D*: `analysis/dnet/lxplus/dnet4/dn4_ana.py`; kombinacja: `analysis/sfabs/skrypty/sfd_final.py` | `docs/results/sfabs/RAPORT.md`, `out/final.json` | nie |
| świetlność, regiony trygera | `analysis/lumi/skrypty/{a1…a4,l1…l4}`; bity i preskale: `analysis/tagtrig/` (+ `lxplus/tagtrig/tt_cool.py`); regiony: `analysis/tagreg/` | `docs/results/lumi/RAPORT.md`, `lumi.json`; `docs/results/tagreg/v2/RAPORT.md` | test `tagtrig`: tak; reszta nie |
| f_acc (tag przypadkowy) | `charmbkg/afp_mix.py` + `analysis/afp_mix/run_data.py` (ekstrakcja: `afp_extract.py`, `afp_merge.py`, `lxplus/afp_ext/`) | `docs/afp_mixing.md`, `docs/results/afp_mix/` | **tak** (toy i test; dane wymagają parquet z ekstrakcji) |
| odczyt z tagiem: N_tag, Z | `analysis/tagged/skrypty/{t1_prep,t2_unblind,tg_lib}.py` (silnik `fb_fit`) | `docs/results/tagged/RAPORT.md` | nie |
| R^fid, σ_tag^fid | `analysis/tagxsec/skrypty/tx_ana.py` (+ `lxplus/tx_mc.py`, `tx_excl.py`) → 13,6 µb; przeliczenie do 13,9 µb: `analysis/tagxsec/skrypty/przeliczenia.py` | `docs/results/tagxsec/RAPORT.md`, `docs/TAGGED_ADDENDUM.md` | nie |
| A·ε_AFP | `analysis/afp_acceptance/skrypty/{p1_extract,p3_analyze,a1_tables,a2_model}.py` | `docs/results/afp_acceptance/` | nie |
| przewidywania Pythii | `analysis/sd_var/` (`lxplus/sd_var/pythia_sd.cc`, `analyze_var.py`, `predict.py`), wcześniej `analysis/sd_pred/` | `docs/results/sd_var/RAPORT.md`, `RAPORT_domkniecie.md`, `docs/results/sd_pred/RAPORT.md` | nie (generacja na lxplus) |
| SF izolacji wg trygera | `analysis/b_trig/` (`btrig_*.py`, `ISO_origin/`, `TRKJET*/`) | `docs/results/b_trig/` | nie |
| testy sieci | `analysis/dnet/lxplus/dnet1…dnet6` | `docs/results/dnet/RAPORT.md` | nie |
| cechy klasyfikatora | `analysis/{feat_charmnet,dlz_dplus,dlz_d0,cn_clone,sig_constr}/` | `docs/results/<temat>/RAPORT.md` | nie (część plików pominięta — patrz README tematu) |
| tło treningowe, SALAD, HYP (stan z 15.11) | pakiet `charmbkg` (`bkggen/`, `s8salad`, `hyp*`); `analysis/{b1_salad_v1,salad_427929}/` | `docs/signal_like_background.md`, `docs/hyp.md`, `docs/results/b1_salad_v1/` | pakiet i testy: **tak**; liczby na danych: nie |
| rysunki | `analysis/presentation/skrypty/`, `docs/figures/` | `docs/results/presentation/PRZEGLAD.md` | nie |

### Jak odtworzyć

1. **Odczyt liczb**: każda liczba z tabel wyżej stoi w pliku wskazanym w kolumnie „źródło” / „raport i liczby”; małe pliki JSON z
   wynikami są obok raportów (`out/`, `tabele/`). Arytmetykę końcową przekroju, σ = N / (C·SF·BR·L), można przeliczyć z
   `docs/results/xsec_final/v2/out/final.json`.
2. **Z repo działa**: pakiet i testy (`pytest`), toy i test sterownika f_acc (`analysis/afp_mix/test_run_data.py`), test ekstrakcji
   bitów trygera (`python analysis/tagtrig/test_tt.py`, wymaga `uproot`, `awkward`).
3. **Fity i tabele** (`analysis/<temat>/`): wymagają histogramów i tabel z obszaru roboczego (`charm_work/…`, npz / h5), których w repo
   nie ma; skrypty importują się wzajemnie przez `sys.path` do sąsiednich katalogów roboczych. Kolejność i wejścia opisuje raport
   danego tematu. Jednego sterownika ani testu regresyjnego dla liczb nagłówkowych nie ma (zalecenie H2 §6.1 pkt 1 — otwarte).
4. **Wszystko przed histogramami** (score per kandydat, ekstrakcja z ntupli, tabele prawdy, generacja Pythii) liczono na lxplus / EOS
   pakietami z `analysis/<temat>/lxplus/`; ntupli w repo nie ma.

## Tło treningowe (stan z 15.11.2026)

**Główna hipoteza — „lepsze tło treningowe → większe Z” — nie potwierdziła się.** Z @WP70 (D⁰, MC closure): sidebandy 4,93,
SALAD 4,96, FETA 4,83, idealne tło 4,95; flow −7 %, kopuła −26 %. Klasyfikator jest nieczuły na mismodelling tła tej skali.

| element | wynik | status |
|---|---|---|
| SB daleki + weto refleksji | najlepsze tło treningowe | **nominał treningu** (używany w charmtag) |
| SALAD v1.0 (MC→dane, warunkowany m) | C2ST 0,507 (D⁰) / 0,510 (D⁺); ΔZ ≈ 0 | **systematyka kształtu tła** ±5/10/25 % @WP80/70/50 i estymator ε_B pod pikiem |
| FETA, flow p(x\|m) + CATHODE | closure 0,60–0,61 na CPU, nie bije SB | cross-check; FETA czeka na trening GPU |
| candrot, pseudo-WS, mixing toru, AD dla par | negatywne | odrzucone |
| strona sygnałowa (D⁰ z tagu D*, sWeights) | SF_rel 0,97 ± 0,05 (8–30 GeV) | zamknięte → systematyka ε_S (nowsze SF: tabela wyżej) |
| HYP: arbiter PID K↔π | SF arbitra 0,977 ± 0,041, nie rzeźbi; ≈ hash po poprawce r(m) | pary: hash, pojedyncze D⁰: keep_both — [docs/hyp.md](docs/hyp.md) |

Otwarte z tego etapu:
1. **Anomalia EV** (5–8 GeV, WP95): SF95 = 0,83 ± 0,10 w oknie ±55 MeV (fit pass/fail, nominalny;
   kontrola sb 0,74 ± 0,06; rdzeń ±36 MeV SF 0,85), płaski w WP; nie z heurystyki swapów ([docs/refl_truth.md](docs/refl_truth.md),
   [ev_dstar](docs/results/ev_dstar/EV_pdgid.md)). Run 427929 daje 1,05 ± 0,11 i nie rozróżnia „0,83” od „1”
   ([b_trig](docs/results/b_trig/runs/lxplus/out/RAPORT_DANE.md)). Raport źródłowy `charm_work/B2b_dstar/` nie jest w repo
   (zawiera odniesienia do projektu A). Dotyczy też systematyki okna masy w charmtag.
2. **Closure par z wyborem nogi hashem** — czy efekt r(m) występuje w tle kombinatorycznym fit2d (MC).
3. FETA na pełnej statystyce (GPU) — tylko jako cross-check SALAD.
4. Systematyka szablonu refleksji ±10 % → ~±5 % (R/S dane/MC @WP70 0,96 ± 0,04; `HYP_d0/t4/t4_results.md` §4–5, poza repo);
   kalibracja PID per tor z tagu D*.

Szczegóły i liczby: [docs/signal_like_background.md](docs/signal_like_background.md).

## Instalacja i szybki start
```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"                            # flow/FETA: .[flow]; wszystko: .[all]
ruff check .
pytest -q -m "not slow"                            # ok. 1–4 min
pytest -q analysis/afp_mix/test_run_data.py        # sterownik f_acc na toy (parquet)
```
Stan testów 20.11.2026 (Python 3.14, czyste środowisko z `.[dev]`): `ruff check .` czysty; `pytest -q -m "not slow"` → 85 passed, 8 skipped, 1 deselected (`slow`), 0 failed. Pominięte: 5 testów integracji z charmtag
(wymagają `charmtag`; z `charmtag` 2.1.0 przechodzą), 1 × `torch`/`zuko` (`.[flow]`), 2 × brak `s1_control.json`
(`scripts/studies/s1_control.py`). CI (`.github/workflows/tests.yml`): Python 3.10 i 3.12.

CLI: `charmbkg bkggen --method sideband|reflection|dstar_wrong|mix|generate|validate|candrot …` (= `charmbkg-bkggen …`),
`charmbkg s8salad …`, `charmbkg hyp-extract …`. API generatorów:
`from charmbkg.bkggen.api import make_generator` → `make_generator(name).fit(SB).sample(m)`; jeden test dwóch próbek: `charmbkg.bkggen.c2st`.

## Układ repo
- `charmbkg/` — pakiet: `bkggen/` (generatory i kontrole tła), `hyp*` + `dstar_tag` (HYP), `refldata`, `bkgsrc`, `s8salad`, `val_epsb`, `sf_rel`;
  rdzeń skopiowany z charmtag 2.0.3 (`config`, `io`, `samples`, `features`, `model`, `massfit`, …).
- `charmbkg/afp_mix*` — tło przypadkowego tagu AFP z mieszania protonów ([docs/afp_mix_method.md](docs/afp_mix_method.md)).
- `configs/` — SALAD v1.0, warianty tła; `scripts/studies/` (s1_control); `tests/`.
- `analysis/` — skrypty poza pakietem: sterujące z testami (`afp_mix/`, `tagtrig/`) i zarchiwizowane kroki analizy, jeden katalog
  na temat z krótkim `README.md` ([analysis/README.md](analysis/README.md)).
- `docs/results/<temat>/` — raporty, tabele, małe pliki JSON i rysunki każdego kroku (po polsku, tak jak powstały); indeks:
  [docs/results/README.md](docs/results/README.md).
- `docs/` — dokumenty metod (`signal_like_background.md`, `hyp.md`, `refl_truth.md`, `afp_mix_method.md`, `afp_mixing.md`), szkic noty
  (`METHODS_NOTE_DRAFT.md`), dodatek o tagu (`TAGGED_ADDENDUM.md`), ustalenie
  o strumieniu danych (`data_stream_428770.md`).

## Zasady
Bez ntupli, h5, modeli i danych binarnych; małe wyniki (JSON / CSV do 1 MB, rysunki) są w `docs/results/`.
Kod ATLAS na GitLabie CERN (`charmpy`, `charmpp`, `CharmAnalysis`) — tylko odczyt, bez kopii w tym repo.
Nic z projektu A (pary D–D) — ani kodu, ani wyników.
