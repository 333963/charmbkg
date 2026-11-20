# D23 — kandydaci D z pT 2–3 GeV dla projektu B: skąd wziąć, paczka, opłacalność (18.11.2026)

Tylko odczyt lokalnych materiałów.

## 1. Gdzie są pliki
| miejsce | co | dostęp |
|---|---|---|
| grid (rucio, scope `user.malewick`, SCRATCHDISK) | ROOT CharmAnalysis `charm-prod-v2.0`: `user.malewick.<próbka>.CharmAnalysis.Dpt2to3GeV.v2.0.Oct28_EXT0`; pliki `…charm_minbias_dpt2to3.root` | **tak, bez nowych uprawnień**, o ile proxy `~/private/x509up_charm` żyje (wystawione 16.11 ~21:20 na 96 h → do ~20.11 21:20) i repliki nie wygasły |
| NERSC `/global/cfs/cdirs/atlas/wcharm/minbias_run3/{grid_v2.0, merged_v2.0_02-11-26, merged_v2.1_04-11-26}` | ROOT + h5 charmpp `<próbka>_d23.h5` | brak dostępu (NERSC) |
| nasz EOS `/eos/user/k/klyczak/charm_raw_v20` | tylko `Dpt3to150GeV` (`rucio_v20/datasets*.txt`) | d23 nie ma (do potwierdzenia `ls`) |

Próbki d23 wg `charmpp:scripts/minbias_v2/samples_v2.0.txt` i `docs/CHARMPP_v2.0_PRODUCTION.md` §1 (mirror, stan 05.11):
data_428770 (97 plików, zadanie 53347060, 206,4 M kandydatów), fd_nd (`mc23.forcedDecay.ND2`, 24, 26,0 M), mb_nd (`mc23.minbiasND`, 80, 72,1 M),
oraz fd_sd/dd/cd, mb_sd/dd/cd. **data_427929 d23: niepewne** — wiersz zakomentowany, zadanie 1.11 było w toku (≈ 450/583 plików); wersja d150 istniała 16.11 (458 plików), więc d23 prawdopodobnie też, ale tego nie widać lokalnie.
Selekcja d23 różni się od d150 wyłącznie `dPtMinGeV/dPtMaxGeV` (2/3 zamiast 3/150; CharmAnalysis 6007843).

Do rozstrzygnięcia przez `bin/lx` (tylko odczyt):
1. `bin/lx 'voms-proxy-info -file ~/private/x509up_charm -timeleft; ls /eos/user/k/klyczak/ | grep -iE "raw|h5"; find /eos/user/k/klyczak/charm_raw_v20 -name "*dpt2to3*" | wc -l'`
2. `bin/lx 'export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase; source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh -q >/dev/null; lsetup rucio >/dev/null; export X509_USER_PROXY=~/private/x509up_charm RUCIO_ACCOUNT=klyczak; for d in data.minbias.428770 data.minbias.427929 mc23.minbiasND mc23.forcedDecay.ND2; do echo "== $d"; rucio list-files user.malewick:user.malewick.$d.CharmAnalysis.Dpt2to3GeV.v2.0.Oct28_EXT0 | tail -3; done'` (liczba plików i rozmiar; pusty wynik = kontenera nie ma / wygasł)
3. `bin/lx 'eos root://eosuser.cern.ch quota /eos/user/k/klyczak | tail -4; fs lq ~; condor_q -batch | tail -3'`

Jeśli proxy < 10 h: odnowić na lxplus `voms-proxy-init -voms atlas -valid 96:00 -out ~/private/x509up_charm` (hasło certyfikatu wpisywane ręcznie).
Jeśli kontenery wygasły: potrzebne przedłużenie repliki (`rucio add-rule`) albo kopia `grid_v2.0/*Dpt2to3GeV*` / `merged_v2.1_04-11-26/*_d23*.h5` na EOS.

## 2. Co zbudowano (nie wysłano)
- `Q_lxplus/h5_d23/` — kopia `h5_v20` (konwerter i meta identyczne bajt w bajt) + pobranie wzorowane na `rucio_v20`; opis w `README.md`.
  Wejście `/eos/user/k/klyczak/charm_raw_v20_d23/<próbka>` (osobny katalog), wyjście `/eos/user/k/klyczak/charm_h5_v20_d23/parts/<próbka>`.
  Zmienne ścieżek: `d23_common.sh` (`PTPROD`, `RAW20`, `H5V20`), `submit_rv.sh` (`RAWD`), `rv_job.sh` (`DESTROOT`), listy `datasets_d23.txt`, `samples_d23.txt`.
- `ORCH/stages/rucio_d23.sh`, `ORCH/stages/h5_d23.sh` (dwa etapy, bo plików nie ma na EOS — odstępstwo od planu „jeden etap”).
- Wpis YAML (NIE dopisany):
```yaml
  - name: rucio_d23     # RUCIO_D23 (18.11, projekt B): pobranie ROOT v2.0 dpt2to3 (428770, 427929, mb_nd, fd_nd) przez rucio w condor → /eos/user/k/klyczak/charm_raw_v20_d23 (stages/rucio_d23.sh)
    max_hours: 30
  - name: h5_d23        # H5_D23: ROOT v2.0 dpt2to3 → h5 (schemat h5_v20) → charm_h5_v20_d23/parts (stages/h5_d23.sh); projekt B, osobno od dpt3to150
    needs: [rucio_d23]
    max_hours: 12
```
- Testy: `bash -n` (7 skryptów paczki, 2 etapy), `py_compile`, suchy DAG w piaskownicy z atrapami plików (28 węzłów, 659 procesów + finish; plik `dpt3to150` pominięty), `DRYRUN=1 convert_d23.sh`. Nie testowane: rucio, condor, rzeczywisty plik d23 (RSS przy 2× większej liczbie kandydatów — ustawiono STEP 50MB / 8 GB).
- Szacunek zasobów: ROOT ≈ 130–140 GB, h5 ≈ 90–100 GB; EOS 17.11: 280 GB + 43 GB z 1 TB.
- Po h5 potrzebne jeszcze (nieprzygotowane): obniżenie twardego pT ≥ 3 w `charmtag.samples.in_range`, klasyfikator/WP dla pT < 3 (H trenowany przy pT ≥ 3).

## 3. Ocena opłacalności (rząd wielkości)
Podstawa: MC `DLZ_dplus/SCAN_2026-11-17/cand/train/oof/base.npz` (D⁺), dane `LOWPT_calib/tables.md` T1, tabela kandydatów z dokumentacji produkcji v2.0.
- **Sygnał D⁺:** zrekonstruowany sygnał MC na GeV rośnie ku progowi coraz wolniej (3,9 → 3,1 GeV: 41,7k → 49,2k/GeV); ekstrapolacja płaska ≈ 50k/GeV daje
  S(2–3)/S(3–5) ≈ 0,6 (3–5: 79,3k). Czynnik niepewności ~1,5 (wydajność poniżej 3 GeV może się załamać).
- **Tło D⁺:** mb na GeV rośnie ×1,45 co 0,25 GeV; ekstrapolacja: B(2–3)/B(3–5) ≈ 3,6. Zgodne z liczbami z dokumentacji produkcji v2.0: mb_nd d23/d150 = 72,1/22,8 M = 3,2; dane 428770 206/109 M = 1,9 (wszystkie gatunki, d150 = cały zakres 3–150).
- **Wniosek D⁺:** S/B × ~0,17, istotność × ~0,3 względem 3–5 GeV przy tej samej sile odrzucania tła (której poniżej 3 GeV nie mamy).
  3–5 GeV, WP70: 427929 N/ε = 4331 ± 257 (Z ≈ 17, S/B 0,14), 428770 5917 ± 429 (Z ≈ 14, S/B 0,06) →
  2–3 GeV: N/ε ≈ 2,6k i 3,5k, **Z ≈ 4–5, S/B ≈ 0,01–0,025**. Pik D⁺ widoczny na granicy; pozycja z błędem kilku MeV.
- **Kalibracja ε_S skanem WP:** już w 3–5 GeV „nieużyteczna” (czułość 21–30 %); w 2–3 GeV ~3× gorzej → bez szans.
- **D\* (tag D⁰: odbicie, SF, weto):** π_s ma pT ≈ 0,07·pT(D*) = 0,14–0,21 GeV, czyli w większości pod progiem 0,2 GeV. W 3–5 GeV mamy 314 D* (427929, WP70), m(Kπ) z tagu Z 3,4–4,8.
  W 2–3 GeV spodziewam się kilkudziesięciu D* na run → **odbicie i wszystko oparte na tagu D\* niemierzalne** (szacunek z kinematyki, nie z MC).
- **Inkluzywne D⁰:** niestabilne już poniżej 8 GeV bez klasyfikatora; w 2–3 GeV nie liczyć.

**Rekomendacja:** pobrać mimo to — koszt mały (~1–2 h condor, ~230 GB EOS), repliki SCRATCHDISK i proxy wygasają, a MC fd_nd/mb_nd d23 jest jedynym sposobem, by zamienić powyższe ekstrapolacje w liczby
(wydajność, udział odbicia z prawdy, trening klasyfikatora < 3 GeV). Realny wynik na danych: pik D⁺ na ~4–5σ na run i udokumentowane „niemierzalne” dla odbicia/ε_S/tagu D*. Wariant oszczędny: najpierw tylko fd_nd + mb_nd (`SAMPLES`/skrócone txt; ROOT ~40 GB, h5 ~30 GB), dane po obejrzeniu MC — ale wtedy ryzyko wygaśnięcia replik danych.
