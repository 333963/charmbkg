# B1 — SALAD v1.0: zamrożony generator tła pod pikiem + closure adwersarialny (2026-11-04)

Projekt B (obiekty podobne do sygnału), strona tła. Cel: (1) zamrozić SALAD jako artefakt produkcyjny v1.0 (estymator ε_B pod
pikiem + systematyka kształtu), (2) sprawdzić, czy closure S4/S8 nie jest kołowy — dotąd pseudo-dane zniekształcano **gładko** w m,
a SALAD zakłada właśnie gładkość w m. Kod/wyniki: `charm_work/B1_salad_v1/`, addytywnie w `charm_work/AA_merge/src`.
Liczone na maszynie roboczej na podpróbkach S6 (`S6_kand/sub/dp_*.npz`, D⁺ po wecie D*, zakres 1,70–2,20): MC tła mb_p0+p1 truth ≥ 8
2 116 565 kand. (SR 373 057), fd p0+p1 truth 0/1 47 594, dane data_p0 851 802 (SB 620 329). D0 nieliczone (czas) — kod gatunkowo
ogólny (`regions.D0` w yaml).

## 1. Artefakt v1.0 (co i gdzie)

| Plik (AA_merge/src) | Rola |
|---|---|
| `config/bkggen/salad_v1.yaml` | JEDNA konfiguracja: kind `hgb` (prod) + `sk_mlp` (cross-check), 18 cech H, t=(m−m0)/0,1, clip 0,05–20, seed 0, max_n 400 k, truth ≥ 8, kryteria akceptacji (n_eff/N(SR) ≥ 0,5, N_SB ≥ 20 k, C2ST pasma ≤ 0,53), regiony D0/D⁺ (SR/pasmo/SB) |
| `charmtag/bkggen/salad_v1.py` | `load_cfg`, `regions`, `SaladV1` (HGB/MLP z t), `fit_salad_v1(cfg, sp, data_sb, mc)`, `salad_weights(rw, mc)`, `run(...)` → `salad_w_<sp>_<sample>_v1.npz` (w float32, `row`, meta: sha256 kolumn MC i danych SB, wersja, kind, n_eff/N, AUC_SB, C2ST pasma, `ok`), `load_weights`; generator api `salad_v1` |
| `charmtag/bkggen/api.py` | hak: 4 linie na końcu (`try: from . import salad_v1`) — rejestruje `make_generator("salad_v1")`; nic innego nie zmienione |
| `tests/test_salad_v1.py` | 4 testy, ~55 s: konfig/regiony; closure gładki syntetyk (C2ST < 0,52, round-trip wag/sha256); closure adwersarialny syntetyk — liczby referencyjne (regresja): skład lokalny pred/true **0,615**, globalny 0,983; rejestracja api |

Bieg produkcyjny na danych (podpróbka data_p0, MC mb_p0+p1) — `B1_salad_v1/out/salad_w_Dplus_data_p0_{hgb,sk_mlp}_v1.npz` (10 MB każdy):

| kind | AUC_SB dane vs MC | n_eff/N (SR) | kwantyle w (SR) 1/50/99 % | C2ST pasma (dane vs MC ważone) | ok | S8 (ten sam test) |
|---|---|---|---|---|---|---|
| **hgb** | 0,643 | 0,67 | 0,35 / 0,84 / 3,78 | **0,507** | tak | 0,510, n_eff/N 0,68 |
| sk_mlp | 0,640 | 0,53 | — | 0,515 | tak | 0,531, 0,74 |

## 2. Closure: gładki (powtórka S8 tym samym kodem) vs adwersarialny

Pseudo-dane = MC B parzystość 1 (po zdarzeniach) z zniekształceniem, akceptacja-odrzucenie; „MC” = parzystość 0; prawda = pseudo SR.
ε_B: HGB fd(p0) vs pseudo-SB, progi WP z fd(p1) w binach pT 3–5–8–30. C2ST: HGB OOF 3-krotny, 60 k/stronę (podłoga 0,50 ± 0,003).
Tryby: **smooth** = S8 (w = exp(Σ(c0+c1·t)·z_j)); **cat_global** = skład per truth_category stały w m (MisMatch* 8–11 ×1,5, OtherNoB 12 ×0,8,
Other*B 13–15 ×1,5); **cat_local** = to samo ×g(m), g = exp(−(m−m0)²/2σ²), σ = 35 MeV (pod pikiem, niewidoczne w SB; f_MM w SR 8,9 → 13,6 %);
**cat_step** = MisMatch* ×1,5 tylko dla m < m0; **cuts_local** = Lxy > 0,3 ×1,5 i |IPsig| > 3 ×0,7 lokalnie (g(m)); **cat_local+smooth**.
Motywacja (V §2 i pomiar tutaj): ε_B@80 w SR per kategoria: OtherNoB 0,011, MisMatchNoB 0,057, MisMatchSomeB 0,098, MisMatchAllB 0,48,
Other*B 0,091 → po cięciu WP80 tło pod pikiem to ~42 % MisMatch* / 54 % OtherNoB mimo 9 / 90 % przed cięciem; AUC MisMatch vs OtherNoB (H, SR) 0,695.

| tryb | kind | f_MM SR pseudo/MC | C2ST raw SB | raw SR | SALAD SB | **SALAD SR** | pasmo | SB-kopia | podłoga | n_eff/N | ε_B pred/true @80 (SALAD / SB-kopia / MC raw) | @70 | @50 | σ_stat prawdy |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| smooth (S8) | hgb | 0,089/0,089 | 0,636 | 0,582 | 0,501 | **0,520** | 0,514 | 0,634 | 0,501 | 0,87 | **1,006** / 0,714 / 0,905 | **1,031** / 0,667 / 0,874 | 1,099 / 0,612 / 0,799 | 0,025/0,038/0,077 |
| smooth (S8) | sk_mlp | 0,089/0,089 | 0,636 | 0,582 | 0,513 | **0,515** | 0,511 | 0,634 | 0,501 | 0,88 | 0,887 / 0,714 / 0,905 | 0,889 / 0,667 / 0,874 | 0,936 / 0,612 / 0,799 | jw. |
| cat_global | hgb | 0,152/0,089 | 0,508 | 0,508 | 0,500 | **0,501** | 0,504 | 0,590 | 0,499 | 1,00 | 0,907 / 0,751 / 0,821 | 0,861 / 0,703 / 0,750 | 0,802 / 0,689 / 0,658 | 0,017/0,025/0,049 |
| **cat_local** | hgb | 0,136/0,089 | 0,502 | 0,508 | 0,503 | **0,507** | 0,503 | 0,598 | 0,500 | 1,00 | **0,845** / 0,647 / 0,846 | **0,809** / 0,584 / 0,810 | 0,701 / 0,508 / 0,702 | 0,017/0,025/0,049 |
| cat_local | sk_mlp | 0,136/0,089 | 0,502 | 0,508 | 0,504 | **0,507** | 0,505 | 0,598 | 0,500 | 0,99 | 0,839 / 0,647 / 0,846 | 0,804 / 0,584 / 0,810 | 0,708 / 0,508 / 0,702 | jw. |
| cat_step | hgb | 0,108/0,089 | 0,504 | 0,506 | 0,504 | **0,503** | 0,504 | 0,588 | 0,501 | 1,00 | 0,938 / 0,756 / 0,925 | 0,929 / 0,720 / 0,913 | 0,852 / 0,661 / 0,833 | jw. |
| cat_step | sk_mlp | 0,108/0,089 | 0,504 | 0,506 | 0,504 | **0,506** | 0,504 | 0,588 | 0,501 | 1,00 | 0,964 / 0,756 / 0,925 | 0,963 / 0,720 / 0,913 | 0,900 / 0,661 / 0,833 | jw. |
| cuts_local | hgb | 0,088/0,089 | 0,502 | 0,536 | 0,502 | **0,536** | 0,511 | 0,610 | 0,501 | 1,00 | 0,917 / 0,696 / 0,916 | 0,949 / 0,677 / 0,948 | 0,939 / 0,654 / 0,939 | 0,015/0,024/0,048 |
| cat_local+smooth | hgb | 0,137/0,089 | 0,637 | 0,583 | 0,503 | **0,508** | 0,509 | 0,632 | 0,500 | 0,88 | 0,792 / 0,576 / 0,750 | 0,754 / 0,503 / 0,675 | 0,656 / 0,389 / 0,516 | 0,030/0,043/0,080 |

Per bin pT (SALAD/prawda @80 / @70): smooth hgb 3–5: 1,02/1,04, 5–8: 0,96/1,00, 8–30: 0,69/0,72 (mało statystyki: pseudo SR 52 k);
cat_local hgb 3–5: 0,85/0,80, 5–8: 0,86/0,88, 8–30: 0,67/0,91. Pełne tabele: `B1_salad_v1/out/b1_closure_Dplus.md`, `closure_*.json`.

**Wnioski z closure**
1. Gładki (S8) odtworzony 1:1 na tym samym kodzie: hgb C2ST SR 0,520 (S8 0,542 na innej podpróbce), ε_B 1,01/1,03/1,10; sk_mlp C2ST lepszy
   (0,515) ale ε_B 0,89/0,89/0,94 (S8: 0,98/0,93/1,04) → **rozrzut sk_mlp między podpróbkami ~10 %**; hgb stabilniejszy w ε_B → hgb = prod.
2. **Closure nie jest kołowy w sensie trywialnym, ale C2ST jest ślepy na skład.** Zniekształcenie składu (Δf_MM = +4,7 pp pod pikiem, cat_local)
   daje C2ST 0,507 (≈ podłoga) w SR, paśmie i SB — niemierzalne — a ε_B pred/true **0,845 / 0,809 / 0,70** @80/70/50. Nawet skład **globalny**
   (ten sam w SB, cat_global) SALAD koryguje tylko w połowie (MC raw 0,82 → 0,91 @80): klasyfikator dane-vs-MC ma AUC 0,514 w SB i nie rozdziela
   4,7 pp domieszki kategorii o AUC 0,695 — a ε_B siedzi w ogonie, gdzie wag nie ma z czego nauczyć.
3. Czułość: **−3,3 % / −4,1 % / −6,4 % ε_B @80/70/50 na +1 pp MisMatch* pod pikiem** (cat_local: Δf = +4,7 pp). Próg w m (cat_step: lewy SB
   ×1,5) → 0,94/0,93 (hgb), 0,96/0,96 (sk_mlp, gładszy w m interpoluje lepiej). Cięcia fizyczne lokalne → 0,92/0,95/0,94 (SALAD = MC raw).
4. Co **jest** widoczne na danych: ε_B(SR) z fitu pol2 SB→okno po cięciu (S8 „data_eff”): SALAD/dane 1,02 ± 0,03 / 1,09 ± 0,04 / 1,25 ± 0,13
   (hgb), 0,97 / 1,01 / 1,14 (sk_mlp). To ogranicza nieznany skład w danych do |Δf_MM| ≲ 1,5 pp @80 (→ ±5 %) i jest jedynym testem, który
   składu nie przepuszcza. C2ST pasma (0,507) jest warunkiem koniecznym, nie wystarczającym.

## 3. Test wycieku sygnału (fd wstrzyknięte w 70–90 MeV od m0, obie strony, cechy sygnału; smooth, hgb)

| wyciek w 70–90 MeV | n_inj | C2ST SALAD SR | pasmo | n_eff/N | ε_B pred/true @80 / 70 / 50 | w̄ MC „sygnałopodobnego” (score ≥ WP80) SR / pasmo / SB | corr(w_SR, bez wycieku) / rel. rms Δw |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0,520 | 0,514 | 0,87 | 1,006 / 1,031 / 1,099 | 1,075 / 1,076 / 1,106 | 1 / 0 |
| 2 % | 478 | 0,527 | 0,522 | 0,86 | 1,030 / 1,070 / 1,170 | **1,156** / 1,158 / 1,104 | 0,949 / 0,136 |
| 5 % | 1232 | 0,515 | 0,512 | 0,88 | 1,011 / 1,061 / 1,232 | 1,117 / 1,115 / 1,135 | 0,877 / 0,197 |

SALAD **uczy się** sygnału z krawędzi SB i ekstrapoluje go do SR: 2 % domieszki → +2,4 / +3,8 / +6,5 % ε_B i +7,5 % wagi zdarzeń sygnałopodobnych
w SR; C2ST nie reaguje (0,527 vs 0,520). Realny wyciek D⁺: ogon fd w 70–90 MeV ≪ 0,1 % (S/B w SR ≈ 0,2 %), truth 2–7 w lewym SB 0,14 % →
efekt < 0,3 % ε_B @80. Wniosek: weto refleksji/D* przed fitem SALAD obowiązkowe (jak S8), a SB_min = 70 MeV nie zawężać.

## 4. Werdykt (kryteria BP §4 B1: adwersarialny C2ST ≤ 0,53 i |ε_B pred/true − 1| ≤ 10 % @80/70)

- C2ST ≤ 0,53: **spełnione trywialnie** (0,507) — kryterium puste dla zniekształceń składu; nie używać jako jedynego.
- |ε_B − 1| ≤ 10 %: **niespełnione** dla Δf_MM = +4,7 pp (15,5 % / 19 %); spełnione dla progu w m (6 % / 7 %) i cięć lokalnych (8 % / 5 %).
- Decyzja: **SALAD v1.0 zostaje estymatorem ε_B pod pikiem, bez nowego generatora**, z ograniczeniem zapisanym w yaml i tutaj:
  *SALAD zakłada gładkość składu tła w m; składnik zlokalizowany pod pikiem przechodzi do ε_B bez korekty: −3,3/−4,1/−6,4 % na +1 pp MisMatch*
  (+4,7 pp → −15/−19/−30 %), a C2ST tego nie wykrywa.* Systematyka kształtu tła: **±5 % @80, ±10 % @70 (S8, z fitu danych — bez zmian),
  @50 podnieść do ±25 %** (dane 1,25 ± 0,13; adwersarialny −30 %). Warunek użycia v1.0 w produkcji: przy każdym zastosowaniu kontrola
  ε_B(SR) z fitu pol2 danych (S8 `data_eff`) — zgodność w ±5/±10 % @80/70 jest kryterium akceptacji, C2ST pasma ≤ 0,53 tylko warunkiem wstępnym.
  Model (hgb vs sk_mlp) różni się o ≤ 10 % w closure gładkim → mieści się w ±5/10 % tylko dla hgb; sk_mlp = cross-check, nie nominał.

## 5. Otwarte
- ~~D0 closure adwersarialny~~ — **zrobione (etap B1-D0, §6)**: D0 spełnia oba kryteria; nowe otwarte: (i) `mc_truth_min: 2` dla D0 w `salad_v1.yaml`
  (per gatunek; zmiana konfiguracji → wersja 1.1) + powtórka `tests/test_salad_v1.py`; (ii) skład truth 2–7 w danych D0 (1 % tła) — normalizacja z MC mb
  (Pythia) nie jest kontrolowana; przy v1.2 (`hypo_swap`, pdgId córek) zmierzyć udział częściowo dopasowanego charmu w SB danych.
- Bin pT 8–30 niedoszacowany (0,67–0,69) we wszystkich trybach — statystyka pseudo SR; powtórzyć na pełnym MC na lxplus (S5_SALAD=1 z S8).
- Wejście v1.2 (`hypo_swap`, pdgId córek) pozwoli zmierzyć f_MM w danych wprost (kontrola składu) — wtedy czułość z pkt 3 zamienia się w systematykę liczoną.

## 6. D0 — closure adwersarialny (etap B1-D0, 2026-11-04)

Powtórka §1–3 dla D0 tym samym kodem SALAD v1.0 (`salad_v1.py`/`salad_v1.yaml`, `regions.D0`: SR |Δm| < 55, pasmo 55–75, SB > 75 MeV),
wejście `S_bkggen/sub/*.npz` (h5 v1.1, decayType 10) → cechy H `charmtag.features.compute`, zakres 1,75–2,00 GeV, pT 3–30, **weto refleksji
|m_swap−m0| > 25 MeV wszędzie** (S4), swapy 100/101 z `closest_mass`. Po wecie: MC mb_p0–p3 1 614 333 kand. (truth ≥ 8: 1 592 528, SR 706 627;
truth 2–7: 17 104 = **1,05 % tła w SR / 1,07 % w SB**; swapy 100/101: 652 = 0,025 % w SR; 0/1: 4 701), fd truth 0/1 106 173 (+ 14 573 swapów 100/101),
dane data_p0+p1 743 240 (SB 295 k). Kod/wyniki: `charm_work/B1_salad_v1/out_D0/` (`code/b1_{prep,study,data,dataeff}_d0.py`, `closure_*.json`,
`b1_closure_D0.md`, `data_eff_D0*.json`, `salad_w_D0_*_v1.npz`). Maszyna robocza, 2 CPU, ~2,5 min/bieg (27 biegów).

**Kategorie truth dla D0 — te same kody co D⁺** (ustalone z danych mb: udziały w SR / mediana Lxy): 8: 4,4 % / 0,44; 9: 0,56 % / 0,47; 10: 0,04 %; 11: 0,28 % / **0,89**
(B-podobna); 12 (OtherNoB): **92,4 %** / 0,44; 13: 0,58 % / 0,46; 14: 0,01 %; 15: 0,02 % / **0,75**. MisMatch* (8–11) = **5,4 %** w SR i 5,4 % w SB (D⁺: 8,9 %).
ε_B@80/70/50 w SR per kategoria (MC surowe, klasyfikator H fd vs pseudo-SB): tło 0,090/0,041/0,008; **MisMatch 0,163/0,096/0,034** (×1,8 tła; D⁺ ×5);
OtherNoB 0,085/0,038/0,006; Other*B 0,180/0,107/0,035; **swapy po wecie 0,91/0,82/0,64** (= sygnał); **truth 2–7 0,20/0,16/0,10** (częściowo dopasowany charm).

### 6a. Produkcja na danych (D0, mb_p0+p1 truth ≥ 8)

| kind | dane | AUC_SB | n_eff/N (SR) | kwantyle w (SR) 1/50/99 % | C2ST pasma | ok | S4 |
|---|---|---|---|---|---|---|---|
| **hgb** | data_p0 | 0,648 | 0,64 | 0,33 / 0,83 / 4,05 | **0,510** | tak | 0,515 |
| sk_mlp | data_p0 | 0,645 | 0,61 | 0,31 / 0,81 / 4,11 | 0,516 | tak | 0,523 |
| hgb | data_p0+p1 | 0,652 | 0,62 | — | 0,507 | tak | — |

Kontrola na danych (jak S8 `data_eff`, kod `b1_dataeff_d0.py`): klasyfikator H fd(p0) vs dane SB (parz. 0), progi z fd(p1), ocena na danych parz. 1;
ε_B(SR) z fitu pol2 SB → okno przed/po cięciu (sygnał poza fitem): **0,0733 ± 0,0021 / 0,0356 ± 0,0014 / 0,0088 ± 0,0007** @80/70/50.
Predykcja / dane: **SALAD hgb 1,03 / 0,99 / 0,84**, sk_mlp 0,98 / 0,94 / 0,78, SB-kopia 0,98 / 0,96 / 0,88, MC surowe 1,25 / 1,19 / 0,95 (MC truth ≥ 8);
**z MC truth ≥ 2** (patrz 6c real_t27): hgb **1,05 / 1,03 / 0,98**, sk_mlp 1,01 / 0,96 / 0,83. Jak S4 (SALAD mlp +6–8 % @80/70 vs fit danych; tu hgb w ±5 %).

### 6b. Closure gładki (smooth, jak S4/S8) — pseudo = MC parz. 1 zniekształcone, MC = parz. 0, prawda = pseudo SR

`w_smooth` z B1 (siła 0,8) daje dla D0 C2ST raw SB tylko **0,568** (D⁺ 0,636, S4 D0 0,631, dane vs MC 0,648) → policzone też **siła 1,6: raw SB 0,645 = skala danych**.

| siła | kind | n_pseudo SR | C2ST raw SB/SR | SALAD SB | **SALAD SR** | pasmo | SB-kopia | n_eff/N | ε_B pred/true @80 (SALAD / SB-kopia / MC raw) | @70 | @50 | σ_stat |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0,8 | **hgb** | 172 827 | 0,568 / 0,560 | 0,501 | **0,527** (seed 1: 0,520) | 0,527 | 0,549 | 0,92 | **1,028** / 0,914 / 1,065 | **1,040** / 0,873 / 1,082 | 1,003 / 0,776 / 1,032 | 0,008/0,012/0,027 |
| 0,8 | sk_mlp | 172 827 | jw. | 0,508 | **0,508** | 0,504 | 0,549 | 0,92 | 0,960 / 0,914 / 1,065 | 0,959 / 0,873 / 1,082 | 0,926 / 0,776 / 1,032 | jw. |
| 1,6 | **hgb** | 90 630 | 0,645 / 0,631 | 0,501 | **0,541** (seed 1: 0,548) | 0,540 | 0,604 | 0,74 | **1,029** / 0,855 / 1,272 (seed 1: 0,985) | **1,018** / 0,785 / 1,285 (0,959) | 0,996 / 0,657 / 1,246 (0,921) | 0,012/0,018/0,039 |
| 1,6 | sk_mlp | 90 630 | jw. | 0,517 | **0,528** | 0,521 | 0,604 | 0,69 | 0,960 / 0,855 / 1,272 | 0,932 / 0,785 / 1,285 | 0,887 / 0,657 / 1,246 | jw. |

S4 (mlp, inna podpróbka/zniekształcenie): C2ST 0,507, ε_B 0,98/0,94/0,88 — odtworzone w granicach rozrzutu: sk_mlp 0,96/0,96/0,93, hgb 1,03/1,04/1,00.
Rozrzut hgb↔sk_mlp **7 % @80, 8 % @70, 8–11 % @50**; rozrzut seedów hgb przy sile 1,6: 4/6/8 %. Per bin pT (hgb 0,8, @80/@70): 3–5 1,03/1,04, 5–8 1,10/1,12,
8–30 1,08/1,05 (D0 ma 3× więcej statystyki w 8–30 niż D⁺ w B1 — bin nie jest już niedoszacowany).

### 6c. Closure adwersarialny (hgb; sk_mlp gdzie podane; C2ST podłoga 0,50 ± 0,003)

| tryb | f_MM SR pseudo/MC | f_swap SR | C2ST raw SB / SR | SALAD SR | pasmo | SB-kopia | ε_B pred/true @80 (SALAD / SB-kopia / MC raw) | @70 | @50 |
|---|---|---|---|---|---|---|---|---|---|
| cat_global (MM ×1,5, OtherNoB ×0,8, Other*B ×1,5, stałe w m) | 0,096/0,054 | 0 | 0,500 / 0,505 | 0,501 | 0,497 | 0,526 | 0,973 / 0,928 / 0,965 | 0,953 / 0,888 / 0,939 | 0,862 / 0,803 / 0,827 |
| **cat_local** (jw. × g(m), σ 35 MeV; Δf_MM = +2,8 pp) | 0,082/0,054 | 0 | 0,496 / 0,501 | **0,501** | 0,499 | 0,532 | **0,972** / 0,900 / 0,972 (seed 1: 0,976) | **0,949** / 0,856 / 0,949 (0,958) | 0,850 / 0,729 / 0,849 (0,861) |
| cat_local sk_mlp | 0,082/0,054 | 0 | jw. | 0,503 | 0,502 | 0,532 | 0,960 / 0,900 / 0,972 | 0,936 / 0,856 / 0,949 | 0,840 / 0,729 / 0,849 |
| cat_local_MM (tylko MisMatch ×2 × g; Δf_MM = +3,5 pp) | 0,088/0,054 | 0 | 0,496 / 0,501 | 0,504 | 0,499 | 0,532 | 0,965 / 0,893 / 0,965 | 0,940 / 0,827 / 0,940 | 0,842 / 0,691 / 0,842 |
| cat_step (MM ×1,5 dla m < m0) | 0,067/0,054 | 0 | 0,499 / 0,501 | 0,502 | 0,497 | 0,525 | 0,986 / 0,933 / 0,991 (sk_mlp 0,993) | 0,970 / 0,903 / 0,976 (0,984) | 0,912 / 0,803 / 0,918 (0,948) |
| cuts_local (Lxy > 0,3 ×1,5, |IPsig| > 3 ×0,7 × g) | 0,053/0,054 | 0 | 0,499 / 0,537 | 0,537 | 0,505 | 0,555 | 0,970 / 0,906 / 0,971 | 0,956 / 0,865 / 0,957 | 0,917 / 0,803 / 0,919 |
| cat_local+smooth (0,8) | 0,081/0,054 | 0 | 0,567 / 0,563 | 0,517 | 0,511 | 0,550 | 0,956 / 0,856 / 1,037 | 0,926 / 0,796 / 1,021 | 0,832 / 0,611 / 0,919 |
| cat_local+smooth (1,6) | 0,081/0,054 | 0 | 0,646 / 0,629 | 0,541 | 0,533 | 0,610 | 0,939 / 0,772 / 1,237 | 0,898 / 0,713 / 1,234 | 0,802 / 0,551 / 1,110 |
| **swap_local ×1,5** (swapy 100/101 z mb+fd parz. 1, udział nominalny z mb × (1+0,5 g)) | 0,054/0,054 | 0,00018 | 0,499 / 0,498 | 0,498 | 0,499 | 0,534 | **1,000** / 0,950 / 1,001 (sk_mlp 1,020) | **1,003** / 0,927 / 1,004 (1,030) | 0,975 / 0,900 / 0,976 (1,027) |
| **swap_local ×0,5** | 0,054/0,054 | 0,00010 | 0,501 / 0,500 | 0,500 | 0,500 | 0,527 | **1,002** / 0,956 / 1,003 (seed 1: 1,000) | **0,998** / 0,932 / 1,000 (0,994) | 1,003 / 0,935 / 1,006 (0,983) |
| swap_plus (+1 pp·g(m) swapów; efektywnie **+0,30 pp** w SR — dziura |Δm| < 25 MeV po wecie) | 0,054/0,054 | 0,00299 | 0,498 / 0,498 | 0,498 | 0,501 | 0,530 | 0,971 / 0,929 / 0,971 | 0,947 / 0,881 / 0,946 | 0,802 / 0,776 / 0,801 |
| **real_t27** (pseudo = truth 2–7 + ≥ 8 + 100/101 parz. 1, bez zniekształcenia; MC truth ≥ 8) | 0,053/0,054 | 0,00025 | 0,500 / 0,502 | 0,501 | 0,495 | 0,526 | 0,984 / 0,961 / 0,984 | 0,971 / 0,951 / 0,969 | **0,879** / 0,950 / 0,865 |
| **real_t27, MC truth ≥ 2** (`mc_truth_min: 2`) | 0,053/0,053 | 0,00025 | 0,500 / 0,500 | 0,501 | 0,501 | 0,526 | **0,999** / 0,961 / 0,999 | **1,003** / 0,951 / 1,002 | **0,991** / 0,950 / 0,990 |
| swap_plus, MC truth ≥ 2 | 0,054/0,053 | 0,00299 | 0,502 / 0,503 | 0,500 | 0,499 | 0,530 | 0,985 / 0,929 / 0,986 | 0,977 / 0,881 / 0,980 | 0,911 / 0,776 / 0,920 |

σ_stat prawdy 0,005–0,010 / 0,008–0,015 / 0,02–0,03 @80/70/50. Pełne tabele (per bin pT, w̄ sygnałopodobnych) w `out_D0/b1_closure_D0.md`.

### 6d. Wyciek sygnału (fd truth 0/1 wstrzyknięte w 75–95 MeV od m0, obie strony)

| siła smooth | wyciek | n_inj | C2ST SALAD SR / pasmo | n_eff/N | ε_B pred/true @80 / 70 / 50 (bez wycieku) | w̄ MC sygnałopodobnego SR / pasmo / SB |
|---|---|---|---|---|---|---|
| 0,8 | 2 % | 1 296 | 0,513 / 0,507 | 0,94 | 1,025 / 1,042 / 1,045 (1,028 / 1,040 / 1,003) → **+0 / +0 / +4 %** | 0,985 / 0,985 / 0,920 (0,980 / 0,979 / 0,913) |
| 0,8 | 5 % | 3 341 | 0,530 / 0,523 | 0,92 | 1,108 / 1,184 / 1,437 → +8 / +14 / +43 % | 1,073 / 1,076 / 0,954 |
| 1,6 | 2 % | 732 | 0,589 / 0,589 | 0,67 | 1,136 / 1,157 / 1,189 (1,029 / 1,018 / 0,996) → **+10 / +14 / +19 %** | 0,901 / 0,905 / 0,778 (0,817 / 0,816 / 0,755) |

Przy zniekształceniu o skali danych (1,6) wyciek 2 % kosztuje +10/+14/+19 % — więcej niż D⁺ (+2,4/+3,8/+6,5 %), bo klasyfikator SALAD przy AUC 0,70 ma mniej
„zapasu” i n_eff/N spada do 0,67. **Realny wyciek D0 jest jednak pomijalny**: ogon fd truth 0/1 poza 75 MeV = 3,3 % sygnału, S/B(SR) ≈ 0,6 % (mb) → sygnał w 75–95 MeV
≈ 0,02 % tła ≪ 2 % → < 0,2 % ε_B. Prawdziwym „wyciekiem” D0 jest **truth 2–7 (częściowo dopasowany charm, 1,05 % tła w SR i tyle samo w SB)** — patrz real_t27.

### 6e. Czułość Δε_B/ε_B na +1 pp składnika pod pikiem (σ 35 MeV)

| składnik | @80 | @70 | @50 | źródło | D⁺ (B1) |
|---|---|---|---|---|---|
| MisMatch* | **−1,0 %** (analit. −0,8) | **−1,7 %** (−1,3) | **−4,6 %** (−3,3) | cat_local_MM (Δf +3,5 pp) / (ε_c−ε_B)/ε_B | −3,3 / −4,1 / −6,4 % |
| Other*B | −1,0 % | −1,6 % | −3,4 % | analit. | — |
| truth 2–7 | −1,2 % (zmierz. −1,6 przy 1,05 pp) | −2,8 % (−2,9) | −11 % (−12) | real_t27 / analit. | (0,14 % w SB, pominięte) |
| swapy K↔π po wecie | −9 % (zmierz. −10 przy 0,30 pp) | −19 % (−18) | −80 % (−66) | swap_plus / analit. | n/d |

Swapy są sygnałem (ε 0,91 @80), ale po wecie refleksji zostaje ich **0,025 % tła w SR** (mb; fd: 3,9 % sygnału w SR → 0,6 % × 3,9 % ≈ 0,02 %, zgodnie) i mają dziurę
|Δm| < 25 MeV → ×1,5 / ×0,5 = ±0,01 pp → **|Δε_B| < 0,3 % @80/70, < 2,5 % @50 (= σ_stat)**. Weto 25 MeV jest wystarczające; swap_local nie jest systematyką.

### 6f. Werdykt D0 (kryteria BP §4 B1) i porównanie z D⁺

- **C2ST ≤ 0,53**: spełnione dla wszystkich trybów adwersarialnych (0,50–0,54; cuts_local 0,537 i smooth 1,6 hgb 0,541/0,548 na granicy — widoczne, bo zniekształcenie
  dotyka cech, nie składu); jak dla D⁺ **ślepe na skład** (cat_local 0,501, real_t27 0,501, swap 0,498).
- **|ε_B pred/true − 1| ≤ 10 % @80/70**: **spełnione dla D0 we wszystkich trybach** (najgorzej cat_local+smooth 1,6: 0,939/0,898 → 6 %/10 %; cat_local 0,972/0,949;
  cat_global 0,973/0,953) — dla D⁺ było niespełnione (0,845/0,809).
- **D0 jest ~3× mniej czuły na skład niż D⁺**: MisMatch* to 5,4 % tła (D⁺ 8,9 %) i jest tylko ×1,8 bardziej sygnałopodobne od OtherNoB (D⁺ ×5, bo 3-ciałowy MisMatch
  z prawdziwym D w środku ma Lxy/pointing sygnału; dla 2 ciał „MisMatch” to głównie jeden prawdziwy tor + losowy) → −1,0 %/pp vs −3,3 %/pp @80; +2,8 pp daje −3 % zamiast −15 %.
- **Nowy składnik specyficzny dla D0: truth 2–7** (częściowo dopasowany charm: D0→Kππ⁰ bez π⁰, D*→D0π itd.), 1,05 % tła w SR, **niewidoczny dla C2ST (0,501) i
  nieskorygowany przez SALAD** (pred/true = MC raw 0,984/0,971/0,879), mimo że jest *gładki* w m (1,05 % SR / 1,07 % SB) — klasyfikator dane-vs-MC z AUC 0,502 nie uczy
  się 1 % domieszki (jak cat_global w B1). To samo widać **na danych**: SALAD/dane @50 0,84 (MC ≥ 8) → 0,98 (MC ≥ 2). **Naprawa = konfiguracja, nie nowy generator:
  dla D0 `mc_truth_min: 2`** (MC tła SALAD = wszystko poza truth 0/1, tj. 2–7, 8–15 i swapy 100/101) → closure 0,999/1,003/0,991; dane 1,05/1,03/0,98.
  Wymaga zmiany yaml (per gatunek) — tu NIE zrobione (etap analityczny); dla D⁺ truth 2–7 = 0,14 % w SB (B1) → bez znaczenia, można ujednolicić.
- **Systematyka SALAD dla D0 (rekomendacja)**: **±5 % @80, ±10 % @70** (jak D⁺/S8: rozrzut modelu 7–8 %, cat_local −3/−5 %, dane 1,03–1,05/0,99–1,03),
  **@50 ±15 %** z `mc_truth_min: 2` (closure 0,99, dane 0,98 ± 0,08, cat_local −15 %, rozrzut modelu 11 %) lub **±20 %** przy MC truth ≥ 8 (dane 0,84; real_t27 0,88) —
  wobec ±25 % dla D⁺. Kryterium akceptacji jak w §4: ε_B z fitu pol2 danych w ±5/10 % @80/70 (D0 spełnia: hgb 1,03/0,99 i 1,05/1,03); C2ST pasma ≤ 0,53 tylko warunek wstępny (0,507–0,516).
- hgb = nominał także dla D0 (sk_mlp 0,96/0,96/0,93 w closure gładkim, 0,98/0,94/0,78 na danych; przy sile 1,6 bin pT 8–30 0,77).
