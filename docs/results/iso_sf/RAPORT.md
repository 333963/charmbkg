# ISO_sf — współczynnik skali izolacji D⁰ z tagu D* (16.11.2026)

**Pytanie.** Jak duży jest SF_iso = ε_dane/ε_MC dla cięć izolacji (ptcone40/100, bezwzględnie i względem pT(D)) dla prawdziwego D⁰?
Czy przeważenie MC do krotności danych domyka SF? (pozycja M3 CharmNet: „data/MC of the inputs”)

## Krótko
- **Statystyka danych jest nadal mała.** Kolumny ptcone są tylko w 2 lokalnych plikach h5 v1.2 (`data_428770_000002`, `_000025`).
  Duża próbka z tagiem D* (`HYP_d0/runs/lxplus_dstar`, 39 plików, 2,2 M wierszy) **nie ma kolumn izolacji ani n_trk**.
  Sygnał D* w 2 plikach (z fitu): **176 ± 26** (preselekcja score_H > 0,6) / 264 ± 68 (bez preselekcji).
  Wszystkie 39 plików dałyby ok. 20× więcej, czyli ok. 3,5–5 tys. D*, σ(SF) ≈ 0,03–0,05. Do tego potrzebne są ptcone z h5 v1.2 na EOS
  (`charm_h5_v12/parts/data_428770`). Skrypty przyjmują dowolną liczbę plików.
- **SF (pT 3–30, preselekcja A).** ptcone40 < 1: **0,45 ± 0,12**; < 2: 0,60 ± 0,10; < 4: 0,71 ± 0,09.
  ptcone40/pT < 0,25: 0,56 ± 0,11. ptcone100 < 4: 0,22 ± 0,09.
  Potwierdza to wcześniejsze 0,19 vs 0,49 (FEAT_charmnet), z inną metodą: izolacja wiersza D⁰, jednoczesny fit pass/fail.
- **Krotność wyjaśnia mniej więcej połowę efektu.**
  - Sygnał D* w danych siedzi w znacznie aktywniejszych zdarzeniach. Udział n_trkid ≥ 15: dane 0,71 ± 0,06, MC 0,23.
  - Po przeważeniu MC do rozkładu n_trkid sygnału w danych: SF(ptcone40 < 1) 0,45 → **0,69 ± 0,19**; ptcone40 < 2: 0,80 ± 0,14;
    ptcone40/pT < 0,25: 0,80 ± 0,16; ptcone40/pT < 0,5: 0,87 ± 0,13.
  - Niezależna mapa w1(n_trk) (1D, etap1) daje to samo: 0,69 / 0,77 / 0,76.
  - Reszta SF ≈ 0,7–0,85 pozostaje (około 1–2σ od 1).
  - ptcone100 nie domyka się (SF_rw 0,4–0,5) → **NO-GO**.
- **Rekomendacja.**
  - Izolacja w CharmNet tylko pod dwoma warunkami: (1) MC przeważone w krotności (w2(n_trk, lead_pt)) oraz (2) rezydualny SF z danych.
  - Dziś SF(ptcone40-typu) = 0,7–0,85 przy błędzie ±0,15–0,20 (20–30 %). To za dużo na nominał.
  - Tylko zmienne ptcone40 (lepiej względne ptcone40/pT, luźne progi). Bez ptcone100/80/60.
  - Systematyka do czasu pełnej statystyki: |1 − SF_rw| ⊕ σ, czyli ≈ 25 % na ε_S części izolacyjnej.
  - Pełna próbka (39 plików) jest warunkiem decyzji GO.

## 1. Dane i metoda
- **Dane.** Wiersze D⁰→Kπ (decayType 10) z tagiem D* z `HYP_d0/runs/lxplus_dstar/data_428770.parquet`, złączone z h5 po (file, row).
  - Uwaga: `row` jest globalny dla pliku ≠ 0.
  - Izolacja pochodzi z **wiersza D⁰** (tak jak w CharmNet), a nie z wiersza D*. Pion miękki leży w stożku D⁰ w danych i w MC jednakowo.
  - |m(Kπ) − m_D⁰| < 40 MeV, pT 3–30 GeV.
- **MC.** fd_nd, 8 lokalnych plików, to samo złączenie; truth ∈ {0,1}, Δm 143,5–147,5 MeV; N = 25 075 (A: ok. 21,8 tys.).
- **ε w danych.** `charmtag.dstar_refit.fit_passfail`: jednoczesny fit Δm pass/fail, kształt sygnału z MC truth.
  - Błąd z profilu NLL; kolumna „det” = `determined`.
  - Skala σ sygnału zatrzymuje się na dolnej granicy 0,8: w danych szczyt jest węższy niż w MC. To wspólne dla pass/fail i w stosunku się skraca.
- **Preselekcja A: score_H > 0,6** (XGB H, bez izolacji) — najlepsza istotność (Z = 6,9 vs 4,0 bez cięcia).
  - Uwaga: ε(score > 0,6) w danych ≈ 0,65 vs MC 0,87 — osobny problem modelowania score.
  - **B** = bez preselekcji (kontrola).
- **Proxy krotności.** Brak n_trk w h5, więc użyto **n_trkid** = liczby różnych torów córek wszystkich kandydatów w zdarzeniu.
  - W MC v9 korelacja z nTracksPt500Eta25 wynosi 0,77.
  - Frakcje n_trkid sygnału w danych pochodzą z fitów pass/fail „w binie/poza binem” (5 binów).
  - ε_MC przeważono do tych frakcji; błąd wag z toy.
- **MC v9 z n_trk.** Lokalne ntuple produkcji v9 (`mc23_minbiasND_Dpt3to150GeV_v9`, 41 plików): 1145 truth D*, 10,7 tys. truth D⁰ inkluzywnie.
  - Mają n_trk, tory i ptcone. ε zgodne z h5 fd (ptcone40 < 1: 0,48 vs 0,455; ptcone100 < 4: 0,38 vs 0,38).
  - Wariant **w1(n_trk)** = `etap1_definicje_masy_v20260925_01/weights_ntrk.csv`.
  - Mapa 2D w2(n_trk, lead_pt) (`w2d.npy`) jest tylko na lxplus (gen_norm), lokalnie brak.
  - Uwaga: w1 bez lead_pt przeszacowuje przesunięcie dla zdarzeń z charmem (walidacja etap1: tempo D* 3,0× vs 1,7× dla w2).

## 2. Tabela SF (preselekcja A, pT 3–30; pełne tabele w `out/sf_table.md`)

N_sygn = 176 ± 26.

| cięcie | ε dane | ε MC | **SF** | SF rw(n_trkid) | SF w1(n_trk) |
|---|---|---|---|---|---|
| ptcone40 < 0,5 | 0,096 ± 0,037 | 0,265 | **0,36 ± 0,14** | 0,59 ± 0,23 | 0,61 ± 0,25 |
| ptcone40 < 1 | 0,203 ± 0,053 | 0,451 | **0,45 ± 0,12** | 0,69 ± 0,19 | 0,69 ± 0,19 |
| ptcone40 < 2 | 0,397 ± 0,069 | 0,664 | **0,60 ± 0,10** | 0,80 ± 0,14 | 0,77 ± 0,14 |
| ptcone40 < 4 | 0,605 ± 0,076 | 0,858 | **0,71 ± 0,09** | 0,83 ± 0,11 | 0,80 ± 0,10 |
| ptcone100 < 4 | 0,083 ± 0,034 | 0,374 | **0,22 ± 0,09** | 0,47 ± 0,20 | 0,62 ± 0,27 |
| ptcone100 < 8 | 0,179 ± 0,052 | 0,711 | **0,25 ± 0,07** | 0,40 ± 0,12 | 0,41 ± 0,12 |
| ptcone40/pT < 0,1 | 0,101 ± 0,036 | 0,300 | **0,34 ± 0,12** | 0,55 ± 0,20 | 0,56 ± 0,21 |
| ptcone40/pT < 0,25 | 0,296 ± 0,059 | 0,531 | **0,56 ± 0,11** | 0,80 ± 0,16 | 0,76 ± 0,16 |
| ptcone40/pT < 0,5 | 0,537 ± 0,075 | 0,758 | **0,71 ± 0,10** | 0,87 ± 0,13 | 0,84 ± 0,12 |
| ptcone40/pT < 1 | 0,823 ± 0,082 | 0,918 | **0,90 ± 0,09** | 0,99 ± 0,10 | 0,95 ± 0,10 |
| ptcone100/pT < 1 | 0,248 ± 0,051 | 0,501 | **0,50 ± 0,10** | 0,86 ± 0,19 | 0,91 ± 0,20 |
| ptcone100/pT < 2 | 0,547 ± 0,079 | 0,814 | **0,67 ± 0,10** | 0,90 ± 0,14 | 0,86 ± 0,13 |

- ptcone100 < 1/2 GeV: ε w danych ≈ 0 (fit przy granicy) — progi bezużyteczne.
- Preselekcja B (bez cięcia score, S = 264 ± 68) daje zgodne wartości, ale z błędami 2× większymi.
  Na przykład SF(ptcone40 < 1) = 0,32 ± 0,19 → rw 0,43 ± 0,26; ptcone40/pT < 0,5: 0,72 → 0,83.

**Biny pT (A).** Statystyka na granicy użyteczności: N = 81 ± 21 / 36 ± 11 / ok. 59 (8–30: fit niestabilny, N_err = 0 — nie używać).

| bin | ptcone40 < 1 SF → rw | ptcone40 < 2 SF → rw | ptcone40/pT < 0,25 SF → rw | ptcone40/pT < 0,5 SF → rw |
|---|---|---|---|---|
| 3–5 | 0,70 ± 0,21 → 1,01 ± 0,32 | 0,76 ± 0,18 → 0,96 ± 0,23 | 0,66 ± 0,21 → 0,95 ± 0,30 | 0,71 ± 0,18 → 0,90 ± 0,23 |
| 5–8 | 0,45 ± 0,31 → 0,69 ± 0,48 (N) | 0,89 ± 0,24 → 1,20 ± 0,33 | 0,63 ± 0,27 → 0,89 ± 0,39 | 0,97 ± 0,27 → 1,18 ± 0,34 |
| 8–30 | 0,02 ± 0,10 → 0,03 (?) | 0,22 ± 0,11 → 0,31 ± 0,16 | 0,39 ± 0,11 → 0,55 ± 0,16 | 0,54 ± 0,09 → 0,65 ± 0,11 |

- W 3–5 GeV po przeważeniu SF jest zgodny z 1.
- W 8–30 GeV deficyt izolacji w danych jest największy (SF_rw 0,3–0,65). Wymaga potwierdzenia: fit 8–30 jest niestabilny.

## 3. Krotność
- **Rozkład n_trkid sygnału D*** (A):

| bin n_trkid | f dane | f MC |
|---|---|---|
| [3,15) | 0,29 ± 0,06 | 0,77 |
| [15,25) | 0,25 ± 0,06 | 0,16 |
| [25,35) | 0,13 ± 0,06 | 0,05 |
| [35,50) | 0,24 ± 0,06 | 0,02 |
| ≥ 50 | 0,08 ± 0,06 | 0,005 |

  Sygnał w danych żyje w dużo aktywniejszych zdarzeniach niż w MC, zgodnie z niedoszacowaniem krotności w Pythii (w1/w2).
  Uwaga: wagi w ogonie sięgają ok. 16×.
- **Zależność ε_MC od n_trk jest silna** (MC v9, truth D*, ptcone40 < 1): 0,80 / 0,67 / 0,47 / 0,50 / 0,43 / 0,20
  dla n_trk 0–10 / 10–20 / 20–30 / 30–40 / 40–60 / > 60. Przeważenie w1: ε 0,48 → 0,31.
- **SF w podziale na krotność** (A, pT 3–30):
  - n_trkid ≥ 25: SF = 1,33 ± 0,44 (ptcone40 < 1), 1,12 ± 0,30 (< 2), 1,13 ± 0,34 (rel40 < 0,25) — **zgodne z 1**.
  - n_trkid < 25: SF = 0,33 ± 0,15 / 0,55 ± 0,13 / 0,54 ± 0,14.
  - Część deficytu w „lo” wynika z rozkładu wewnątrz binu: dane są przesunięte ku 15–25.
  - Po pełnym przeważeniu 5-binowym zostaje SF ≈ 0,7–0,85 z błędem ok. 0,15–0,2.
- **Wniosek.**
  - Aktywność zdarzenia jest głównym, mierzalnym źródłem różnicy. Dwie niezależne metody (n_trkid z danych, w1(n_trk)) dają tę samą poprawkę.
  - Domyka ona połowę lub więcej różnicy dla ptcone40, ale nie dla ptcone100.
  - Rezydualny SF ≈ 0,7–0,85 jest na poziomie 1–2σ. Nie da się rozstrzygnąć, czy to efekt fizyczny (np. UE wokół charmu, fragmentacja), czy fluktuacja.
  - n_trkid jest tylko proxy; pełne w2(n_trk, lead_pt) wymaga n_trk w danych (h5 v2.x / lxplus).

## 4. Ograniczenia
- Tylko 2/39 plików danych ma ptcone lokalnie. To główne ograniczenie, a nie metoda.
- Część fitów ma ε < 0 lub nieokreślone: bardzo ciasne progi ptcone100 oraz bin 8–30. Nie używać.
- Skala σ sygnału w danych jest na granicy 0,8. Model tła: |ε_expo − ε| typowo < 0,02–0,04 (`out/sf_fits.json`).
- Wagi n_trkid są jednowymiarowe (bez lead_pt). Proxy zależy od rekonstrukcji kandydatów, więc nie musi skalować się jak n_trk.
- SF mierzone na D⁰ z D* (z pionem miękkim w stożku) przenosi się na D⁰ inkluzywne tylko przy założeniu, że efekt pionu jest dobrze modelowany.

## Pliki
- `scripts/extract_data.py` (dane/fd: złączenie parquet ↔ h5, n_cand, n_trkid), `scripts/evtools.py`
- `scripts/extract_root_mc.py` (MC v9: n_trk, lead_pt, ptcone wiersza D⁰), `scripts/mc_ntrk.py`
- `scripts/sf_fit.py` (106 fitów pass/fail), `scripts/summarize.py`
- `out/sf_table.md` (wszystkie tabele), `out/sf_fits.json`, `out/sf_summary.json`, `out/mc_ntrk.json`, parquety pośrednie
- Logi w `logs/`.

## Pełna próbka — lxplus (ISO-FULL, 16.11, przygotowane)
- Paczka `Q_lxplus/iso_full` (README tam), etap ORCH `iso_full` → `runs/lxplus/{data_428770,fd_nd,mb_nd}_tag.parquet` (te same kolumny co `out/*_tag.parquet` + ptcone60/80, run/event, eta/phi).
- Analiza: `scripts/run_lxplus.sh [NPROC]` → `runs/lxplus/out/` (sf_fit.py/summarize.py czytają `ISO_IN`/`ISO_OUT`; domyślnie `out/` jak dotąd, wyniki 2-plikowe nietknięte; kopie `*.bak_1116`).
- Test: ekstraktor na 2 + 8 lokalnych plikach daje identyczne wiersze i wartości jak `out/*_tag.parquet`.
