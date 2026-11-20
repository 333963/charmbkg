# H+DLZ dla D⁺→Kππ na danych — raport (17.11.2026, projekt B)

**Werdykt wg kryteriów zapisanych przed obejrzeniem piku (`KRYTERIA.md`): NIE PRZYJMOWAĆ H+DLZ (NO-GO).**
Zysk statystyczny jest (K1 PASS), ale nie przechodzi rzeźbienie (K2) ani względny SF (K3) — w obu runach,
także w rozstrzygającym 427929. Wariant zapasowy H+dalitz też nie przechodzi K3.
Pary D–D (projekt A) nie były dotykane.

Pełne tabele: `tables.md`; werdykt maszynowy: `verdict.json`; rysunki: `fig/`; dane wejściowe raportu: `runs/lxplus/`.

## 1. Co zrobiono (wg `DLZ_dplus/PROCEDURA.md` §2–§4)

- Liczone na lxplus, etap ORCH `dlz_data` (paczka `Q_lxplus/dlz_data`, DAG 54 węzły, 08:19–17:12 CET), wyniki na EOS
  `/eos/user/k/klyczak/charm_q/dlz_data`, pobrane do `runs/lxplus/`.
- Trening: 5 foldów po zdarzeniach, sygnał MC fd+mb (truth 0/1), tło = SB danych 428770 (|m−m0| > 70 MeV, po wecie D*), cap 3·10⁶,
  reważenie tła do sygnału w (pT, |η|); modele H (18 cech), H+DLZ (26), H+dalitz; dwa ziarna (12345, 23456) dla H i H+DLZ.
- Progi WP (95/90/80/70/50) z OOF sygnału MC per bin pT (3–5, 5–8, 8–30 GeV) — te same ε_S^MC dla H i H+DLZ.
- Dane h5 v2.0, wszystkie pliki: 428770 — 97 plików, 55 253 987 kandydatów po wecie D*; 427929 — 458 plików, 19 313 580 kandydatów,
  z czego 14 494 522 po `passGRL_MinBiasAFP ∧ HLT_mb_mbts_L1MBTS_1` (tylko te wchodzą do fitów 427929). L_xy ≥ 0,25 mm.
- Fit masy: 2 Gaussy o kształcie z MC (wolne N_S, dμ, skala σ) + exp(pol2), binowany Poisson ML, 5 MeV, 1,70–2,20 GeV.
  Błędy zysku i ΔSF_rel: bootstrap po zdarzeniach, 200 replik, sparowany H/H+DLZ.

## 2. Wynik: zysk σ(N_S)/N_S przy tych samych WP

Zysk = 1 − [σ/N]_DLZ / [σ/N]_H. W binie 3–5 GeV fit swobodny ma parametr kształtu na granicy przy WP90 (oba runy) →
zgodnie z regułą z `KRYTERIA.md` §5 dla K1/K3 obowiązuje fit o kształcie ustalonym na MC (kolumna „kształt ust.”).

**Run 427929 (czysty minbias), WP70**

| pT [GeV] | N_S (H) | N_S (H+DLZ) | S/B H → DLZ | Z H → DLZ | zysk, fit swobodny [%] | zysk, kształt ust. [%] |
|---|---|---|---|---|---|---|
| 3–5 | 3427 ± 298 | 3987 ± 280 | 0,157 → 0,291 | 19,0 → 24,5 | +19,4 ± 5,4 | **+20,7 ± 2,9** |
| 5–8 | 1565 ± 85 | 1717 ± 85 | 0,698 → 1,079 | 26,4 → 29,7 | **+9,1 ± 3,0** | +9,7 ± 2,3 |
| 8–30 | 631 ± 41 | 630 ± 40 | 1,336 → 1,720 | 22,2 → 23,4 | **+4,2 ± 2,9** | +4,4 ± 2,2 |

**Run 428770 (82 % L1_J12), WP70**

| pT [GeV] | N_S (H) | N_S (H+DLZ) | S/B H → DLZ | Z H → DLZ | zysk, fit swobodny [%] | zysk, kształt ust. [%] |
|---|---|---|---|---|---|---|
| 3–5 | 4206 ± 399 | 6426 ± 512 | 0,065 → 0,155 | 15,2 → 23,2 | +16,1 ± 7,5 | **+32,9 ± 3,5** |
| 5–8 | 3036 ± 179 | 3557 ± 187 | 0,217 → 0,352 | 24,2 → 29,3 | **+10,9 ± 3,1** | +16,6 ± 2,2 |
| 8–30 | 6145 ± 177 | 6304 ± 170 | 0,704 → 0,956 | 51,2 → 55,4 | **+6,4 ± 1,4** | +6,5 ± 1,1 |

(pogrubione = wartości użyte w K1). Średnia po binach: 427929 **+11,3 %**, 428770 **+16,7 %** → K1 PASS (próg: ≥ 0 w każdym binie, średnio ≥ 8 %).
Pozostałe WP (95/90/80/50) — `tables.md`; trend: zysk rośnie ku luźnym WP i niskiemu pT, ale tam jest też najmniej wiarygodny (§3).

**Zastrzeżenie do zysku (wynik):** przy tej samej ε_S^MC uzysk N_S powinien być równy, a H+DLZ daje systematycznie więcej:
N_DLZ/N_H przy WP70 = 1,163 ± 0,088 / 1,097 ± 0,040 / 0,997 ± 0,037 (427929) i 1,528 ± 0,135 / 1,172 ± 0,049 / 1,026 ± 0,019 (428770);
przy WP90–95 w 3–5 GeV aż 2–8×. Część „zysku” w 3–5 i 5–8 GeV to więc większy dopasowany N_S, a nie mniejszy błąd.

## 3. Rzeźbienie i stabilność piku

- **R_sculpt w danych** (ε_B w paśmie 45–70 MeV / ekstrapolacja liniowa z |Δm| > 70 MeV), WP70, H → H+DLZ:
  427929: 0,926 ± 0,010 → 1,025 ± 0,014 (3–5); 0,992 ± 0,031 → 1,112 ± 0,041 (5–8); 1,015 ± 0,067 → 1,106 ± 0,083 (8–30);
  428770: 0,929 ± 0,006 → 1,032 ± 0,008; 0,954 ± 0,012 → 1,041 ± 0,016; 0,974 ± 0,016 → 1,054 ± 0,019.
  H+DLZ podnosi R o ~+0,08–0,12 we wszystkich binach. W 3–5 GeV to „bliżej 1”, ale w 5–8 (427929) i 8–30 (428770) odchodzi od 1
  o > 2σ względem H → **K2 FAIL**.
- MC truth ≥ 8, okno sygnału: H 0,913 ± 0,008, H+DLZ 1,002 ± 0,011 (WP70) — ten sam kierunek i wielkość co w danych.
- Nachylenie ε_B(m): H −0,24, H+DLZ −0,30 /100 MeV (różnica 0,06, w dopuszczalnych ±0,10).
- χ²_SB fitu masy H+DLZ przy WP70: 0,68–1,31 (< 1,5, PASS). Ale przy luźnych WP w 3–5 GeV fit H+DLZ jest zły
  (χ²/ndf 1,78–2,36 w 427929, 2,76–6,35 w 428770), a N_S silnie zależy od modelu tła: 428770, 3–5 GeV, WP70:
  exp(pol2) 6426 ± 512, exp(pol3) 5029 ± 414, exp(pol1) 10985 ± 270; dla H rozrzut jest dużo mniejszy (4099–4230 bez pol1).
- **Closure fitu na MC mb** (ε_fit(WP70)/ε_fit(WP90) względem prawdy): H 1,20 ± 0,14 / 0,92 ± 0,08 / 0,98 ± 0,10;
  H+DLZ **0,59 ± 0,07** / 0,81 ± 0,07 / 0,93 ± 0,09. Fit po cięciu H+DLZ nie odtwarza prawdy już na MC w 3–5 GeV.
- **Pik (K5)**: Δμ(DLZ−H) przy WP70 w granicach (|Δμ| ≤ 0,96 MeV tam, gdzie określone); szerokość DLZ/H w 5–8 GeV
  1,098 ± 0,041 (427929) i 1,128 ± 0,050 (428770) — poza 2σ → **K5 FAIL (ostrzeżenie)**; 8–30: 0,990 ± 0,046 / 1,042 ± 0,022.
  3–5 GeV: nieokreślone (reguła granicy).

## 4. SF wydajności (PROCEDURA §3)

ΔSF_rel = SF_rel(DLZ)/SF_rel(H), SF_rel(WP) = [N(WP)/N(WP90)]_dane / [ε(WP)/ε(WP90)]_MC, przy WP70:

| pT [GeV] | 427929 | 428770 | 427929, kształt ust. | 428770, kształt ust. |
|---|---|---|---|---|
| 3–5 | 0,508 ± 0,116 | 0,394 ± 0,086 | **0,728 ± 0,068** | **0,530 ± 0,058** |
| 5–8 | **0,880 ± 0,066** | **0,737 ± 0,087** | 0,962 ± 0,043 | 0,880 ± 0,050 |
| 8–30 | **0,860 ± 0,062** | **0,880 ± 0,030** | 0,936 ± 0,043 | 0,927 ± 0,020 |

Próg: |ΔSF_rel − 1| < 2σ i < 0,10 w każdym binie → **K3 FAIL** w obu runach (427929: 3–5 −4,0σ, 8–30 −2,3σ; 5–8 mieści się w 2σ,
ale nie w 0,10). SF_rel samego H jest zgodny z 1 (427929, WP70: 1,10 ± 0,22 / 1,10 ± 0,09 / 0,97 ± 0,08).
Kontrole: referencja WP95 zamiast WP90 i połowy danych dają ten sam znak (`tables.md`).
SF warunkowy z sPlot (427929, przy H WP90): H WP70 1,08 ± 0,11 / 1,05 ± 0,05; H+DLZ WP70 1,42 ± 0,14 / 1,14 ± 0,06 —
ε „danych” > 1 przy WP80 dla H+DLZ w 3–5 GeV (1,82 ± 0,17) jest niefizyczne i pokazuje to samo co N_DLZ/N_H > 1.

K4 (C2ST sygnału, sWeights): ΔC2ST(H+DLZ − H) = −0,032 ± 0,006 (427929), −0,028 ± 0,003 (428770) → PASS (≤ 0,02).

## 5. Kryteria — zestawienie

| kryterium | 427929 | 428770 |
|---|---|---|
| K1 zysk σ/N (WP70) | PASS (+20,7 / +9,1 / +4,2 %, śr. 11,3 %) | PASS (+32,9 / +10,9 / +6,4 %, śr. 16,7 %) |
| K2 rzeźbienie | **FAIL** (R 5–8 GeV: 1,112 ± 0,041 vs 0,992 ± 0,031) | **FAIL** (R 8–30: 1,054 ± 0,019 vs 0,974 ± 0,016) |
| K3 ΔSF_rel | **FAIL** (0,73 / 0,88 / 0,86) | **FAIL** (0,53 / 0,74 / 0,88) — w tym runie tylko opisowe |
| K4 C2ST sygnału | PASS | PASS |
| K5 pik (ostrzeżenie) | FAIL (szerokość 5–8 GeV) | FAIL (szerokość 5–8 GeV) |
| K6 foldy | PASS: 0 naruszeń; AUC trening − OOF = +0,006…+0,008 | (ten sam model) |
| K7 ziarna, ΔN_S(WP70) H+DLZ | PASS (+0,97σ / −0,22σ / +0,40σ) | PASS (−0,82σ / +0,41σ / +0,63σ) |

Uwaga do K7: sam H w 428770, 3–5 GeV różni się między ziarnami o 3,1σ (4206 vs 5438) — fit w tym binie jest niestabilny niezależnie od DLZ.
Część K7 „ΔZproxy < 3 %” nie była liczona w tym biegu (tylko ΔN_S).

Wariant zapasowy **H+dalitz** (WP70): zysk +12,7 ± 5,7 / +8,4 ± 2,8 / +4,1 ± 2,6 % (427929), ΔSF_rel 0,54 ± 0,11 / 0,92 ± 0,07 / 0,86 ± 0,06 → też FAIL K3.

## 6. Interpretacja (hipotezy, nie wynik)

1. Niedobór ΔSF_rel i nadmiar N_DLZ/N_H mają najpewniej wspólną przyczynę w **ficie, nie w wydajności sygnału**: cięcie H+DLZ
   zmienia kształt tła pod pikiem (R rośnie z ~0,93 do ~1,03; tło przestaje być opisywane przez exp(pol2) przy luźnych WP), więc
   fit przypisuje część tła sygnałowi przy WP90/95 (mianownik SF_rel). Przemawia za tym closure na MC (0,59 w 3–5 GeV) oraz to,
   że ΔSF_rel podzielone przez stosunek closure DLZ/H wychodzi ~1,04 / 1,00 / 0,90 (427929; błędy ~15–20 %, nie liczone ściśle).
2. Alternatywa nieodróżniona w tym badaniu: rzeczywista składowa pikująca wzmacniana przez cechy Dalitza/pT córek
   (refleksje, np. D_s⁺→KKπ, D⁺→πππ, D*-feed). Wymagałoby szablonów refleksji z MC.
3. W 8–30 GeV, gdzie tło jest małe, odchylenie jest najmniejsze (0,86–0,94), ale zysk też (+4–6 %) — H+DLZ niewiele tam daje.

## 7. Odstępstwa i ograniczenia

- Model trenowany na SB 428770 stosowany do 427929 (zdarzenia rozłączne) — zapisane z góry.
- Bin 3–5 GeV: fit swobodny na granicy przy WP90 → K1/K3 z kształtu ustalonego na MC, K5 nieokreślone (reguła z §5 KRYTERIÓW,
  dopisana po teście technicznym na 5 plikach, przed pełnym biegiem).
- G2 (C2ST tła) nie powtarzane (było +0,0056). G5 (closure) opisowe — ale to ono najmocniej tłumaczy FAIL.
- Zysk „rozrzut replik” (kolumna w `tables.md`) bywa inny niż z hesjanu (np. 428770 5–8 GeV WP70: −5,1 % vs +10,9 %) —
  błędy hesjanu w fitach z wolnym kształtem nie są w pełni wiarygodne; to dodatkowy powód ostrożności wobec K1.

## 8. Co dalej (propozycje, do decyzji)

- Nie używać H+DLZ w analizie w obecnej postaci. Jeśli zysk w 3–8 GeV jest wart zachodu: (a) fit z modelem tła odpornym na
  rzeźbienie (szablon tła z SB po cięciu / pol3 + składowe refleksji), (b) powtórzyć K3 z closure ≈ 1 na MC jako warunkiem wstępnym,
  (c) dekorelacja score od masy (np. cechy Dalitza unormowane do m_D są już użyte — sprawdzić, które z 8 cech niosą R: przesiew pojedynczo).

## 9. Pliki

- Kryteria: `KRYTERIA.md`. Skrypty: `scripts/{dlzlib,apply,fitlib,analyze,mc_checks,report}.py` (kopie w paczce `Q_lxplus/dlz_data`).
- Rysunki: `fig/fig1_massfit_WP70_{427929,428770}.png`, `fig2_gain.png`, `fig3_sculpt_*.png`, `fig4_dSF.png`.
- Etap ORCH: `ORCH/stages/dlz_data.sh` (stan: done).
- Pozostałości: `h5/data_427929/` (~4 GB, niepełne pobranie v2.0, już niepotrzebne), `logs/`, puste `models/`, `oof/`;
  na lxplus `/tmp/klyczak_dlzt` (lxplus993, scratch testów).
