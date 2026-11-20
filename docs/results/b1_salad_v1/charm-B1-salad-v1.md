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
- D0 closure adwersarialny (kod gotowy: `b1_prep.py`/`b1_study.py` z `--species`-owym yaml; dane `S_bkggen/sub`).
- Bin pT 8–30 niedoszacowany (0,67–0,69) we wszystkich trybach — statystyka pseudo SR; powtórzyć na pełnym MC na lxplus (S5_SALAD=1 z S8).
- Wejście v1.2 (`hypo_swap`, pdgId córek) pozwoli zmierzyć f_MM w danych wprost (kontrola składu) — wtedy czułość z pkt 3 zamienia się w systematykę liczoną.
