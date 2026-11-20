# ISO-FULL — izolacja D⁰ z tagu D* na pełnej próbce (16.11, przygotowane, NIEWYSŁANE)

Cel: `charm_work/ISO_sf` (SF_iso = ε_dane/ε_MC dla ptcone40/100) miał tylko 2/39 plików danych z ptcone (S ≈ 176).
Ta paczka robi to samo złączenie co `ISO_sf/scripts/extract_data.py` + `evtools.py`, ale na wszystkich plikach h5 v1.2 na EOS.

- Wejścia (tylko odczyt): T4 `charm_q/hyp_d0_dstar/<próbka>.parquet + .meta.json` (D⁰ z tagiem D*, Δm, m, pT, score D0_H z join,
  `row` globalny = offset pliku + wiersz lokalny) i h5 v1.2 `charm_h5_v12/parts/<próbka>/<plik>_none.h5`.
- `iso_extract.py` (numpy/pandas/h5py/pyarrow, bez charmtag): na plik h5 jeden odczyt pól → ptcone40/60/80/100, actual_mu,
  truth_category z **wiersza D⁰**; n_cand i n_trkid (proxy krotności: różne trackId córek wszystkich kandydatów w bloku event_number;
  n_trk nie ma w h5). Kontrola: m bit w bit (|Δ| < 1e-4), decayType 10, liczba wierszy pliku = offsety meta.
  Kolumny: sample, file, row, evblk, run, event, m, pt, eta, phi, dm, score, score_src, fold, Lxy, pis_pt, soft_pi_charge, dstar_tag,
  dstar_conflict, dstar_n, dstar_row, truth, hypo_swap, hypo_from_tag, d_trk0/1, d_pt0/1, ptcone40/60/80/100, mu, truth_cat, n_cand, n_trkid.
- `iso_merge.py` (finish): scala części, sprawdza pokrycie plików i sumę wierszy = T4 → `<próbka>_tag.parquet` + `.meta.json`, `iso_summary.json`, `DONE`.
- Wyjście EOS `charm_q/iso_full/` (części w `parts/`, ORCH ich nie pobiera). Nazwy `<próbka>_tag.parquet` = jak w `ISO_sf/out` → skrypty ISO_sf czytają bez zmian.

## Uruchomienie (lxplus, `~/private/charm_q/Q_lxplus/iso_full`)
`./submit_iso.sh check|submit|status|fetch` — batch `charm_iso`, DAG: 3 węzły dane (po 13 plików) + fd_nd + mb_nd → finish.
Submit w czystym env (`env -u LD_LIBRARY_PATH -u PYTHONPATH`, `condor_submit_dag -f`), check robi `condor_q` i odmawia przy DONE.
Scratch `mktemp` + trap, wyniki tylko EOS, logi `Q_lxplus/logs/iso_full`.
Zmienne: `ISO_OUT`, `TAG_DIR`, `H5V12`, `ISO_SAMPLES` (domyślnie `data_428770 fd_nd mb_nd`), `ISO_NPARTS_DATA` (3), `ISO_DISK` (3GB).
ORCH: etap `iso_full` (`ORCH/stages/iso_full.sh`) → pobiera do `charm_work/ISO_sf/runs/lxplus/`.

## Test lokalny (16.11)
2 pliki danych + 8 fd (`FEAT_charmnet/data`, `--only-local`): 111 383 / 118 641 wierszy = `ISO_sf/out/*_tag.parquet`,
wszystkie wspólne kolumny identyczne (max |Δ| = 0); po scaleniu części też 0. ~10 s/plik danych, ~5 MB/plik (zstd).
Pełna analiza ISO_sf (`scripts/run_lxplus.sh 8` z RUN_DIR=test, ~30 min) odtwarza `out/sf_table.md` (ε, SF identyczne; różnice tylko w 3. cyfrze błędu toy wag — kolejność imap_unordered zmienia losowania).

## Szacunek
Dane: 2,21 mln wierszy → ~200 MB (zstd; gzip w LCG ~250 MB); fd_nd 0,31 mln → ~25 MB; mb_nd 0,39 mln → ~35 MB.
Czas: węzeł danych 13 × (~0,5–1 min przez FUSE) ≈ 10–20 min; MC ~10 min; finish ~2 min. Razem z kolejką ~0,5–1 h.

## Po biegu (lokalnie)
`charm_work/ISO_sf/scripts/run_lxplus.sh [NPROC]` → `ISO_sf/runs/lxplus/out/{sf_fits.json,sf_table.md,sf_summary.json}` (nie rusza `ISO_sf/out`).
