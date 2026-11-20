# SD_pred: przewidywania Pythii 8 dla D + proton w oknie AFP — paczka i wynik testu (17.11.2026)

Zadanie #8 planu prac z 17.11. **Wszystko na poziomie cząstek, przed wydajnościami detektora.** Bez danych ATLAS.
Stan: paczka gotowa i przetestowana w `/tmp` na lxplus (kilka–kilkadziesiąt tysięcy zdarzeń na model). **Joby produkcyjne nie zostały wysłane**;
Liczby w §4 pochodzą z testu i mają duże błędy statystyczne.

## 1. Co już było (nie duplikowane)

| miejsce | co zawiera | czego brak |
|---|---|---|
| `DIFF_model/RAPORT.md` | kształty pT(D) per proces, pseudo-tag z kodu procesu (SD + CD), rozrzut tune'ów A3/Monash/MPI±, lista wariantów do wygenerowania (§6), łatka ξ/t (§5) | ξ, t, M_X w próbkach; okno 0,03–0,10; ND/DD z protonem wiodącym; warianty strumienia; twarda dyfrakcja |
| `AFP_ddpred/` (`gen/A3` 500 × 1 M, `gen/syst_*`) | pary D–D w dyfrakcji wobec ND (poza zakresem tego zadania) | jak wyżej — pliki `gen_*.tsv.gz` nie mają ξ |
| `TUNE_pythia/`, `GEN_mcstat/`, `Q_lxplus/gen_norm` | systematyka tune'u dla μ par, wagi w2 (projekt A) | nie dotyczą tagu protonu; `gen_norm` nie był ruszany |
| `Q_lxplus/truth/pythia_dd.cc`, `Q_lxplus/syst` (etap `tune`) | wzorzec generacji Pythii na lxplus (LCG_110_swan, Pythia 8.317) | — |

`AFP_ddpred/leadscan.py` nie istnieje (ani lokalnie, ani w innych katalogach projektu). Żadna istniejąca próbka nie ma zapisanego protonu,
więc z gotowych plików nie dało się policzyć niczego z listy planu; `DIFF_model` stwierdza wprost, że ND i DD „nie dają nietkniętego protonu” —
test poniżej pokazuje, że dają (proton wiodący z fragmentacji), i to na poziomie ~10 % próbki z tagiem.

## 2. Definicja użyta w paczce

- D: D⁰, D⁺, D*⁺ z antycząstkami, ostatnia kopia, prompt i z b razem (flaga `fromB` zapisana), **2 ≤ pT < 30 GeV, |η| < 2,2**
  (region fitu głównego łańcucha (przegląd kodu charmpp `GITLAB_2026-11-17/RAPORT.md` §(b) i l. 127); biny 2–3, 3–5, 5–8, 8–30). D⁰ z D*⁺ liczone także jako D⁰. BR naturalne.
- Tag: proton końcowy z 0,030 ≤ ξ ≤ 0,100 po stronie A (+z) lub C, ξ = 1 − |p_z|/p_wiązki. „Dokładnie jeden proton w oknie” jest na poziomie cząstek
  spełnione automatycznie (dwa protony z ξ ≤ 0,1 po jednej stronie przekraczają energię wiązki). Strony sumowane jak region `tag` = A lub C;
  D w zdarzeniu z tagiem po obu stronach liczone dwa razy (histogram `pt_tagAC` daje ich udział).
- **Okno t:** w źródłach nie ma okna t akceptancji AFP (tylko „zależna od pT protonu, optyki i odległości od wiązki”, `RECON_program` D2).
  Nominał jest więc bez cięcia t; zapisany jest rozkład |t| (krok 0,05 GeV²) i wariant |t| < 1 GeV² (CDF, CMS-TOTEM), więc dowolne okno da się nałożyć po fakcie.
- ξ z cząstek centralnych: ξ_cal = Σ(E ± p_z)/√s po widocznych cząstkach |η| < 4,9 (+ dla strony A); ξ_trk: naładowane |η| < 2,5, pT > 0,5 GeV; ξ_tk1: pT > 0,1 GeV.
- Bez etykiety procesu w definicji; proces służy tylko do rozbicia na składowe.

## 3. Paczka

`charm_work/Q_lxplus/sdpred/` (lxplus: `~/private/charm_q/Q_lxplus/sdpred`, EOS `/eos/user/k/klyczak/charm_q/sdpred`):

| plik | rola |
|---|---|
| `pythia_sd.cc` | generator (na bazie `pythia_dd.cc`): ustawienia wyłącznie z plików `models/*.cmnd`; TSV tylko dla zdarzeń z D (pT > 2, \|η\| < 2,5): wiersz `E` (proces, flaga twardej dyfrakcji, n_ch, ξ i t protonu wiodącego A/C, ξ_cal, ξ_trk, ξ_tk1 A/C) + wiersze `D`; JSON: σ_gen i liczniki protonów w oknie dla wszystkich zdarzeń per proces |
| `models/*.cmnd`, `models.txt` | 19 modeli: nazwa, liczba jobów, zdarzenia/job, lista cmnd |
| `sd_reduce.py` | w jobie: TSV → histogramy (npz) per składowa × gatunek |
| `sd_merge.py` | węzeł końcowy: suma i suma kwadratów po jobach → `merged/<model>.json` (błąd stat. z rozrzutu między jobami) |
| `sd_job.sh`, `sd.sub`, `make_dag_sd.sh`, `submit_sd.sh` | węzeł, DAG (batch `charm_sd`), `check\|submit\|status\|fetch`; `check` kompiluje `bin/pythia_sd`, robi test doHard 200 zdarzeń, sprawdza `condor_q`, DONE i kwotę AFS; submit w czystym env |

`ORCH/stages/sdpred.sh` — etap (sync, check, submit, done po `DONE` nowszym niż `.orch_start`, fetch tylko `merged/`).
`SD_pred/analyze.py` — lokalnie, tylko json + numpy (~1 s): `SD_pred/out/tables.md`, `results.json`.

### Modele (`models.txt`)

| model | proces | różnica względem bazy A3 | joby × zdarzenia |
|---|---|---|---|
| `nd_A3` | SoftQCD:nonDiffractive | — (proton wiodący z ND; też mianownik R) | 100 × 1 M |
| `nd_Monash` | ND | tune Monash (inna fragmentacja resztek wiązki) | 40 × 1 M |
| `dif_A3` | SD + DD + CD | strumień DL, ε = 0,07, α′ = 0,25 (nominał ATLAS A3) | 20 × 1 M |
| `dif_SaS` | SD + DD + CD | `PomFlux = 1` (Schuler–Sjöstrand) | 10 × 1 M |
| `dif_DL085`, `dif_DL10`, `dif_DL04` | SD + DD + CD | DL z ε = 0,085 / 0,10 / 0,04 | 10 × 1 M każdy |
| `dif_MBR` | SD + DD + CD | `PomFlux = 5`, `SigmaTotal:mode = 2`, `SigmaDiffractive:mode = 2` (σ także z MBR) | 10 × 1 M |
| `dif_H1A`, `dif_H1B` | SD + DD + CD | `PomFlux = 6` / `7` (H1 2006 Fit A / B) | 10 × 1 M każdy |
| `dif_sigPp25` | SD + DD + CD | `Diffraction:sigmaRefPomP = 25` mb (mniej MPI w układzie X) | 10 × 1 M |
| `dif_Monash` | SD + DD + CD | tune Monash (strumień SaS) | 10 × 1 M |
| `hq_incl` | HardQCD cc̄ + bb̄, pTHatMin = 0 | bez doHard (mianownik twardy, proton wiodący bez dyfrakcji) | 8 × 0,5 M |
| `hq_t3`, `hq_t3_SaS`, `hq_t3_H1B` | jak wyżej + `Diffraction:doHard`, `sampleType = 3` | tylko dyfrakcyjne, **bez** tłumienia MPI; strumień DL-A3 / SaS / H1 Fit B | 8 × 0,25 M każdy |
| `hq_t4`, `hq_t4_SaS`, `hq_t4_H1B` | `sampleType = 4` | tylko dyfrakcyjne, **z** dynamiczną przerwą (weto MPI) | 10 × 0,2 M każdy |

Razem 312 jobów, 2,66·10⁸ zdarzeń. Poprawka do `DIFF_model` §6: MBR to `SigmaTotal:mode = 2` (nie 3); 3 = ABMST (sprawdzone w xmldoc 8.317).

## 4. Wynik testu (lxplus `/tmp`, jeden plik na model, ziarno 7; pełne tabele: `test/out/tables.md`)

Test techniczny: kompilacja z Pythią 8.317 (LCG_110_swan), `Diffraction:doHard` typ 3 i 4 **działają**, generator → reduce → merge → `analyze.py` przechodzi.

| model | zdarzenia | czas/zdarzenie | σ_gen [mb] | D fid. z tagiem (A + C) | P(proton w oknie) na stronę |
|---|---|---|---|---|---|
| `nd_A3` | 40 000 | 2,5 ms | 55,59 | 3 | ND 0,095 % |
| `dif_A3` | 40 000 | 0,6 ms | 23,00 | 77 | SD 4,8 %, DD 0,9 %, CD 14,6 % |
| `dif_MBR` | 20 000 | 0,7 ms | 21,50 | 22 | SD 3,1 %, DD 0,65 %, CD 12,8 % |
| `dif_H1B` | 20 000 | 0,9 ms | 23,00 | 36 | SD 4,7 %, DD 0,87 %, CD 15 % |
| `hq_incl` | 20 000 | 1,3 ms | 4,98 | 16 | 0,13 % |
| `hq_t3` | 20 000 | 1,3 ms | 0,750 | 1848 | 11 % |
| `hq_t4` | 10 000 | 1,6 ms | 0,331 | 476 | 10,7 % |

Czas/zdarzenie = czas całkowity/N, zawiera inicjalizację (~8 s), więc dla krótkich biegów jest zawyżony; ND (101 s na 40 000) jest miarodajne.

Liczby fizyczne z testu (błędy stat. 15–60 %, **orientacyjne**):

- **Miękka dyfrakcja + ND (A3):** σ(D⁰ + p, obie strony) = 27,5 ± 4,3 µb, R(D⁰) = 0,81 ± 0,13 %; D⁺ 12,1 ± 2,6 µb (0,72 ± 0,16 %); D*⁺ 8,9 ± 2,5 µb (0,70 ± 0,20 %).
  MBR: R(D⁰) = 0,47 ± 0,13 %; H1 Fit B: 0,83 ± 0,17 %. R jest tu liczone z oknem ξ i na obie strony; literaturowe 1,6–2 % dotyczy ξ < 0,1 bez dolnej granicy.
- **Twarda dyfrakcja (cc̄ + bb̄, strumień DL-A3):** typ 3: R ≈ 4,3–4,7 %; typ 4 (dynamiczna przerwa): R(D⁰) = 1,12 ± 0,10 %, D⁺ 1,19 ± 0,13 %, D*⁺ 1,43 ± 0,21 %.
  **S² = σ(typ 4)/σ(typ 3) w oknie = 0,26 ± 0,02 (D⁰)**, a nie ~0,07 jak dla dijetów: przy skali powabu weto na dodatkowe MPI jest dużo słabsze.
  Cały przekrój (bez okna): 0,331/0,750 = 0,44; udział dyfrakcyjny typu 3 w σ(cc̄ + bb̄) = 15 %.
- **Proton wiodący bez dyfrakcji:** ND daje ~4 µb D + p (3 zdarzenia), czyli rząd 10 % próbki z tagiem w scenariuszu miękkim; w `hq_incl` 4,0 ± 1,0 µb,
  co przy tłumieniu typu 4 stanowi 18 ± 7 % (D⁰). Udział rośnie z ξ (0 poniżej 0,05, ~20–35 % w 0,075–0,10) i protony te mają twardsze t
  (tylko 33–44 % z |t| < 1 GeV², dyfrakcyjne 94–100 %). DD w teście: 0 zdarzeń (za mała statystyka).
- **Asymetria:** SD miękka A_η = −0,25 ± 0,11 (D częściej po stronie przeciwnej do protonu, zgodnie z 60:40 z `DIFF_model`); twarda dyfrakcja ≈ 0 (+0,01 ± 0,02 typ 3, +0,08 ± 0,05 typ 4).
- **ξ_cal:** na poziomie cząstek ξ_cal ≤ ξ_p zawsze (zachowanie E + p_z), więc **cięcie ξ_cal ≤ ξ_AFP + 0,02 ma wydajność 100 % we wszystkich składowych**;
  straty mogą pochodzić wyłącznie z rozdzielczości, szumu i pile-upu. Przesunięcie: mediana ξ_cal − ξ_p = −0,011 (SD miękka; ξ_cal/ξ_p = 0,79), −0,019 (ND),
  −0,023 (typ 3), −0,032 (typ 4; ξ_cal/ξ_p = 0,39). ξ z torów (|η| < 2,5) to tylko ~5–10 % ξ_p — nie nadaje się jako zamiennik bez silnej korekty.
  Zdarzenia ND z D bez wymogu protonu: mediana ξ_cal = 0,11; P(ξ_cal ≤ 0,12) = 0,55 — tyle tła przypadkowego przepuszcza cięcie przy ξ_AFP = 0,10.
- **Oczekiwane liczby (poziom cząstek, przed jakąkolwiek wydajnością, obie strony):** na 1 nb⁻¹ ≈ 27 500 D⁰ + p (A3 miękkie) lub ≈ 9 700 (typ 4);
  po BR kanału: D⁰→Kπ ≈ 1 090 / 380, D⁺→Kππ ≈ 1 130 / 500, D*⁺→D⁰π→Kππ ≈ 240 / 130 na nb⁻¹.
- **Rozrzut modeli (test):** R(D⁰) 0,47–1,12 % (×2,4), σ(D⁰ + p) 9,7–28 µb (×2,9) między strumieniami miękkimi a twardą dyfrakcją typu 4; pełna wartość po produkcji.

## 5. Statystyka i czas

Z testu: D fid. z tagiem na zdarzenie: ND 7,5·10⁻⁵, dyfrakcja miękka 1,9·10⁻³, typ 3: 0,09, typ 4: 0,05, `hq_incl` 8·10⁻⁴.
Plan daje (wszystkie gatunki): ND ~7 500 (D*⁺ ~1 500 → 2,6 %; biny pT 5–8 %), `dif_A3` ~38 000, warianty ~19 000, twarde 100 000–180 000, `hq_incl` ~3 200.
Błąd stat. R na model < 3 % dla każdego gatunku; ogranicza składowa ND (wspólna dla wariantów strumienia).
Koszt: ND 140 × ~42 min, dyfrakcja 110 × ~10–15 min, twarde 62 × 5–12 min → **≈ 130 CPU·h**, ściana 1–4 h przy `-maxjobs 150` (zależnie od kolejki).
EOS: TSV ~4 MB na 1 M zdarzeń ND, razem < 2 GB; pobierane tylko `merged/*.json` (< 30 MB). AFS: tylko logi.

## 6. Uruchomienie

Wpis do `ORCH/pipeline.yaml`:
```
  - name: sdpred        # SDPRED: Pythia 8 D + proton w oknie AFP (poziom cząstek; ND/SD/DD/CD + twarda dyfrakcja) → SD_pred/runs/lxplus; analiza SD_pred/analyze.py (projekt B)
    max_hours: 12
```
Bez `needs` (etap niezależny). Po fetch: `.venv_h5/bin/python charm_work/SD_pred/analyze.py` (domyślnie `SD_pred/runs/lxplus/merged` → `SD_pred/out/`).
Test częściowy: `SD_ONLY=hq_t4,dif_A3 SD_NJOBS=2 SD_NEV=20000 ./submit_sd.sh submit`.

## 7. Zastrzeżenia

1. `SoftQCD:nonDiffractive` + `Diffraction:doHard` (typ 3) nie wygenerowało żadnego zdarzenia w 100 s — ta kombinacja jest nieużywalna. Twarda dyfrakcja
   jest więc tylko dla `HardQCD:hardccbar/hardbbbar` (LO 2 → 2, bez rozszczepienia gluonu i wzbudzenia zapachu z MPI/kaskady), stąd σ_incl(D⁰) = 0,87 mb
   wobec 3,4 mb w SoftQCD. R twarde i miękkie mają różne mianowniki; porównywać R, nie σ.
2. Brak okna t i akceptancji AFP (φ, pT protonu); brak Reggeonu (Pythia go nie ma); POMWIG nie robiony.
3. Proton wiodący w ND zależy od modelu resztek wiązki — jedyny wariant to Monash; rzeczywista niepewność jest większa niż rozrzut dwóch tune'ów.
4. σ_SD w A3 jest przeszacowane względem ALFA (×0,64 w 8 TeV) — normalizacja miękka jest górna; R skaluje się tak samo.
5. ξ_cal na poziomie cząstek nie zawiera progów klastrów ani szumu; wynik „ε = 100 %” mówi tylko, że nieefektywność cięcia jest czysto detektorowa.
6. W scenariuszach twardych składowa `HQnd` (proton wiodący z resztek wiązki) jest brana z `hq_incl` (bez doHard) i dodawana do dyfrakcyjnej z typu 3/4.
   To przybliżenie: w spójnej próbce inkluzywnej (`sampleType = 2`) część tych zdarzeń zostałaby zastąpiona dyfrakcyjnymi (różnica rzędu udziału dyfrakcji, 7–15 %).
7. Liczby w §4 pochodzą z testu jednego pliku na model; błędy przy jednym jobie są poissonowskie.

## 8. Wyniki produkcji (17.11.2026, 11:47) — POZIOM CZĄSTEK

Pełne tabele: `out/tables.md`, `out/results.json` (26 modeli, 18 scenariuszy). Definicja: D z antycząstkami, 2 ≤ pT < 30 GeV, |η| < 2,2; proton 0,030 ≤ ξ ≤ 0,100 po dowolnej stronie; R = σ(D + p)/σ_incl(D); błędy tylko stat. MC.

**Błąd konfiguracji pierwszego biegu (naprawiony etapem `sdpred2`).** W Pythii 8.3 `SigmaDiffractive:PomFlux` (z ε i α′) działa na miękką dyfrakcję tylko przy `SigmaDiffractive:mode = 0` (klasa `SigmaTotOwn`; sprawdzone w nagłówkach 8.317). Modele `dif_SaS`, `dif_DL04/085/10`, `dif_H1A/B` były więc identyczne z `dif_A3`, a `dif_A3` ma kształt SaS, nie DL z tune'u A3. Analiza pomija te sześć duplikatów. Poprawione modele `dif_own_*` (`models/own.cmnd`: mode = 0, σ_SD/DD/CD ustawione na wartości SaS, razem 23,0 mb) zmieniają tylko kształt w (ξ, t). Twarda dyfrakcja (`doHard`) strumień uwzględniała od początku.

| scenariusz miękki (ND z A3 + SD + DD + CD) | R(D⁰) [%] | σ(D⁰ + p) [µb] | udział nie-SD [%] |
|---|---|---|---|
| kształt SaS, mode = 1 (pierwotny „A3”) | 0,750 ± 0,006 | 24,5 ± 0,2 | 10,2 ± 0,2 |
| strumień SaS, mode = 0 | 0,816 ± 0,008 | 26,9 ± 0,3 | 9,3 ± 0,3 |
| strumień DL ε = 0,04 | 0,548 ± 0,007 | 17,9 ± 0,2 | 12,5 ± 0,4 |
| **strumień DL ε = 0,07, α′ = 0,25 (parametry A3)** | **0,356 ± 0,005** | 11,5 ± 0,2 | 16,6 ± 0,5 |
| strumień DL ε = 0,085 | 0,287 ± 0,006 | 9,2 ± 0,2 | 19,7 ± 0,6 |
| strumień DL ε = 0,10 | 0,216 ± 0,004 | 6,9 ± 0,1 | 25,5 ± 0,8 |
| strumień H1 Fit A | 0,111 ± 0,003 | 3,5 ± 0,1 | 43,3 ± 1,5 |
| strumień H1 Fit B | 0,122 ± 0,003 | 3,9 ± 0,1 | 40,5 ± 1,3 |
| MBR (kształt i σ) | 0,431 ± 0,005 | 13,8 ± 0,2 | 13,9 ± 0,4 |
| σ_ref(Pomeron–p) = 25 mb, kształt SaS | 0,319 ± 0,005 | 10,2 ± 0,2 | 22,5 ± 0,7 |

| scenariusz twardy (HardQCD cc̄ + bb̄) | R(D⁰) [%] typ 3 (bez tłumienia) | typ 4 (dynamiczna przerwa) | S² = typ 4/typ 3 |
|---|---|---|---|
| strumień DL-A3 | 4,17 ± 0,02 | 1,06 ± 0,01 | 0,255 ± 0,002 |
| strumień SaS | 6,13 ± 0,02 | 1,49 ± 0,01 | 0,243 ± 0,002 |
| strumień H1 Fit B | 9,52 ± 0,04 | 2,27 ± 0,01 | 0,239 ± 0,002 |

Wnioski:
1. R miękkie zależy od strumienia 7,4× (0,11–0,82 %); razem z twardą typu 4 rozrzut R wynosi 20×. D⁺ i D*⁺ dają te same R co D⁰ w granicach 3 %.
2. Zamknięcie mode 0 wobec mode 1 dla SaS: 0,816 wobec 0,750 (9 %) — tryby nie są identyczne nawet przy tym samym strumieniu; przyczyna niezbadana.
3. R spada z pT(D): DL-A3 0,41 → 0,24 %, SaS 0,87 → 0,54 %, twarda typ 4 1,41 → 0,29 % (biny 2–3 → 8–30 GeV).
4. Tło fizyczne (proton wiodący z ND, plus DD i CD) rośnie z ξ: dla DL-A3 udział ND 6,2 % przy ξ 0,03–0,04 i 17,6 % przy 0,075–0,10; dla H1 B 20 → 49 %. Im stromszy strumień, tym mniej czysta próbka tagowana.
5. Asymetria η względem strony protonu: SD −0,08…−0,12 (D częściej po stronie przeciwnej do protonu), twarda +0,01…+0,03, ND −0,08.
6. Cięcie ξ_cal ≤ ξ_p + 0,02 ma wydajność 100 % dla wszystkich składowych — wynika z zachowania pędu, więc każda nieefektywność w danych jest detektorowa.
7. Punkty 1–7 z §7 pozostają w mocy (brak akceptacji AFP w t i φ, brak Reggeonu, różne mianowniki R miękkiego i twardego).
