# PEAK_LXY010 — D0 z tagiem D* z h5 v2.0 (L_xy ≥ 0,10 mm) → parquet (lxplus/HTCondor), 18.11.2026

Projekt B. Przygotowane, **NIE wysłane**. Wzorzec: `Q_lxplus/hyp_extract` + etap `hyp_dstar` (to on zrobił wejścia PEAK_pos i SWAP_d0:
`HYP_d0/runs/lxplus_dstar`). Różnice względem wzorca:

| | hyp_dstar (wzorzec) | peak_lxy010 |
|---|---|---|
| wejście | h5 v1.2 `charm_h5_v12/parts` (L_xy ≥ 0,25 mm) | h5 **v2.0** `charm_h5_v20/parts` (D0 i D0-w-D*: L_xy ≥ 0,10 mm) |
| próbki | fd_nd, mb_nd, data_428770 | fd_nd, mb_nd, data_428770, **data_427929** |
| score D0_H | join z `scores/D0_H` (wiersze v1.1) | **modele foldów** `models/D0_H` (`SCORE_MODE=model`; wiersze v2.0 są inne) |
| wyjście EOS | `charm_q/hyp_d0_dstar` | `charm_q/peak_lxy010` |
| batch / logi | `charm_hyp`, `logs/hyp` | `charm_peaklxy`, `logs/peaklxy` |
| finish | metadane | + kontrola: `Lxy_min < 0,2` i udział `L_xy < 0,25` > 5 % w każdej próbce, brak NaN w score (inaczej STOP, bez DONE) |

Kod ekstrakcji bez zmian: `src/` = kopia `Q_lxplus/hyp_extract/src` (charmtag 2.0.3 + HYP.diff + `dstar_tag.py`), moduł
`charmtag.hyp_extract --dstar-only`. Wiersze: `decayType 10`, m ∈ [1,75; 2,00] GeV, pT 3–30 GeV, z dopasowanym D* (trackId córek,
Δm < 165 MeV); kolumny jak w `HYP_d0/runs/lxplus_dstar` (m, pt, eta, Lxy, score, dm, soft_pi_charge, hypo_from_tag, dstar_conflict,
truth, hypo_swap, cechy H, córki). Ekstrakcja NIE tnie po L_xy — cięcia 0,10 / 0,25 robi analiza. **Żadnej selekcji AFP** (dane inkluzywne).

## Użycie (lxplus, `~/private/charm_q/Q_lxplus/peak_lxy010`; normalnie przez ORCH `stages/peaklxy.sh`)
```
./submit_pl.sh check                 # wejścia v2.0 (liczby plików), models/D0_H, src, kolejka, DONE, AFS < 85 %
DRYRUN=1 ./submit_pl.sh submit       # tylko pl.dag + komendy węzłów
./submit_pl.sh submit                # condor_submit_dag -f -batch-name charm_peaklxy pl.dag (czysty env)
./submit_pl.sh status | fetch
```
Zmienne: `SAMPLES`, `PL_OUT`, `H5V20`, `NOMINAL` (= `$RUNS`, `charm_q/runs/lxplus_v11`), `PL_EXTRA=--dstar-only`, `DATA_FRAC=1`, `SCORE_MODE=model`.
Zasoby: 4 węzły (2 CPU; MC 6 GB, dane 12 GB; `workday`) + finish (`espresso`). Czas: lokalnie 226 plików 427929 (16,3 M wierszy) = 1,5 min;
na lxplus dominuje odczyt FUSE (428770: 97 plików, ~110 M wierszy D0+D*) — szacunek 0,5–2 h. Wyjście: 4 parquety, razem ~0,5–1 GB
(lokalnie 427929 z 226 plików = 133 MB → pełne ~270 MB; 428770 i fd_nd podobnego rzędu). Scratch `mktemp` + trap, wyniki tylko EOS, AFS tylko logi.

## Czego ta paczka NIE robi
- **Bin pT 2–3 GeV**: `charm_h5_v20` powstało wyłącznie z produkcji `dpt3to150` (`rucio_v20` filtruje `*dpt3to150*`), a `charmtag.samples.in_range`
  ma twarde pT ≥ 3. Potrzebne osobno: pobranie kontenerów `…Dpt2to3GeV.v2.0.Oct28_EXT0` (dane obu runów + mb/fd), konwersja h5 (`h5_v20` z innym
  wzorcem nazw i osobnym katalogiem — te same zdarzenia są w obu produkcjach, nie mieszać), obniżenie pT w `in_range`, a progi WP D0_H dla 2–3 GeV
  nie istnieją (model trenowany na pT ≥ 3). To osobny etap, nieprzygotowany.
- **D⁺**: v2.0 nie zmienia D⁺ (L_xy ≥ 0,5 mm w produkcji, wiersze identyczne z v1.2) — liczby PEAK_pos obowiązują bez zmian dla 428770;
  D⁺ w 427929 wymagałby odpowiednika `dp_extract` na v2.0.
- **D⁰ inkluzywne** (bez tagu): `PL_EXTRA=` (pusty) + `DATA_FRAC`/`MAX_GB` — duże pliki; PEAK_pos pokazał, że inkluzywne D⁰ nie nadaje się do tych systematyk.

## Po biegu (lokalnie)
ORCH pobiera do `PEAK_lxy010_2026-11-18/runs/lxplus/`. Potem `PEAK_lxy010_2026-11-18/RAPORT.md` §6 (podmiana wejść w `fit_lxy.py`, `refl_lxy.py`).
