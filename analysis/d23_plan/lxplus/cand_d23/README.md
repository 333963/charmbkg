# CAND_D23 — kandydaci D z binu pT 2–3 GeV (h5 v2.0 `dpt2to3`) → parquet + histogramy, 18.11.2026 (projekt B). NIE wysłane.

Wejście: `/eos/user/k/klyczak/charm_h5_v20_d23/parts/{data_428770,data_427929,fd_nd}/*_none.h5` (etap `h5_d23`; MC minbias w tym binie nie ma).
Wyjście: `/eos/user/k/klyczak/charm_q/cand_d23/` — `dp_<próbka>.parquet`, `hist_<próbka>.npz`, `summary.json`, `DONE` (+ `parts/` części węzłów, nie pobierane).
Dane inkluzywne, bez selekcji AFP. Tylko pojedyncze D.

| co | treść |
|---|---|
| `dp_<s>.parquet` | D⁺→Kππ per kandydat: decayType 4, m ∈ [1,70; 2,20), **pT ∈ [2, 3)**, L_xy ≥ 0,25 (nieaktywne: produkcja ≥ 0,5); zawetowane D* zostają z flagą `dstar_veto`. Kolumny: `run, lb, event, entry, file_idx, row, evblk` (MC: zdarzenie = (`file_idx`, `entry`), nie (run, event)), `m, pt, eta, phi, charge, m_kpi1/2`, 18 cech H, `score` (model foldu), `score_avg`, `score_sd` (5 foldów), `fold`, `truth` (kategoria), `is_sig` (truth 0/1; dane −1), `D_tpdg`, `n_truth_ok`, `weight`, `nSelBL, mu, mu_act, npu, npv, ntrk, grl`, `trig_<HLT…>` (11), `tbp_<L1…>` (5) (−1 = brak kolumny), `ps_tier, ps_u, w_ps` |
| `hist_<s>.npz` | pełna statystyka: `dp_h[klasa, pT, score, m]` (po wecie D*; 1 MeV; score 140 binów), `d0_h[klasa, pT, m]` (D⁰, obie hipotezy, 1,70–2,05), `ds_h[RS/WS, klasa, pT, m(D0), Δm]`; klasy danych: wszystko / GRL∧mb_mbts / GRL∧L1J12 / GRL; MC: wszystko / truth 0,1 / truth ≥ 8 / inne; krawędzie w pliku |
| `summary.json` | wiersze, pliki, zakres pT (kontrola `pt_in_range`), liczniki D⁺/D⁰/D* (w tym kandydaci poza oknem pT w plikach), rozmiary, przerzedzenie, progi WP binu 3–5 (odniesienie) |

**Score poza zakresem treningu.** Modele `charm_q/dlz_data/train/models/H_s12345/fold{0..4}.json` (te same co dlz_data/SCAN) trenowane na pT ≥ 3 GeV;
`log_pt` < log1p(3) nie występował w treningu. Dane: fold = `fold_of(run, event)` jak `dlz_data/apply.py`. MC: `fold_of(1, file_idx·10⁷ + evblk)` — formuła wzorca,
ale pliki d23 ≠ d150, więc to nie fold treningowy (kandydaci pT < 3 i tak nie byli w treningu); `score_avg` = średnia foldów. WP nie są stosowane.

**Twardy próg pT ≥ 3** obejściem w KOPII `src/charmbkg` (migawka z `dlz_data/src`): `samples.pt_range()` czyta env `CHARMBKG_PT_RANGE="lo,hi"` (brak = 3,30);
używają jej `samples.in_range`, `hyp_extract.count_selected`, `bkgsrc`. `cd_extract.py` ustawia okno z `--pt-lo/--pt-hi` i przerywa, gdy wiersz jest poza oknem. Opis: `src/SNAPSHOT.txt`.

**Preskala wierszy (budżet < 2 GB).** Bin ma ~1,6·10⁸ kandydatów D⁺ → pełny zapis ≈ 25 GB. Zostają: sygnał MC (tier 3) i score ≥ `CD_SC_FULL`=0,93 (tier 2) w całości;
|m − m0| < 60 MeV z p = `CD_PS_WIN`=0,05 (tier 1); reszta p = `CD_PS_OUT`=0,01 (tier 0); wybór haszem zdarzenia `ps_u`, waga `w_ps` = 1/p. `finish` przy
sumie > `CD_MAX_GB`=1,8 przerzedza dalej po `ps_u` (tier 0/1, potem 2; `summary.json: thinning`), > 2,0 GB = STOP bez DONE. Uzysk i S/B w funkcji cięcia: z `hist_*.npz`.

## Użycie (lxplus `~/private/charm_q/Q_lxplus/cand_d23`; normalnie ORCH `stages/candd23.sh`)
```
./submit_cd.sh check            # pliki, DONE h5_d23, modele foldów, kolejka, AFS < 85 %, suchy przebieg na 1 pliku (dane + MC, 60 tys. wierszy)
DRYRUN=1 ./submit_cd.sh submit  # tylko cd.dag + komendy
./submit_cd.sh submit           # condor_submit_dag -f -batch-name charm_candd23 cd.dag (czyste env)
./submit_cd.sh status | fetch
```
DAG: `ext_<s>_<g>` (pliki dzielone na ≤ 12 / 10 / 3 części; 4 CPU, 8 GB, `workday`, RETRY 2) → `finish`. Logi `Q_lxplus/logs/cand_d23`. Scratch `mktemp` + `trap`.
Zmienne: `cd_common.sh` (`H5D23, CD_OUT, DP_MODELS, SAMPLES, CD_NGMAX, CD_PT_LO/HI, CD_SC_FULL, CD_PS_*, CD_MAX_GB`; `FORCE=1`, `CD_ALLOW_NODONE=1`, `CD_NODRY=1`).

## Nie zrobione
Tabele per kandydat D⁰ (score D0_H) i D⁰←D* — tylko liczniki/histogramy. W d23 są wyłącznie D* z pT(D*) ∈ [2, 3), więc tag D* dla D⁰ z tego binu byłby niepełny.

## Test lokalny (18.11)
`bash -n`, `py_compile`; pliki syntetyczne (h5 v2.0 427929 i v1.2 fd_nd z pT przeliczonym do 2–3 GeV, część wierszy zostawiona poza oknem): `make_dag_cd.sh` (6 + finish),
`cd_job.sh ext` × 6, `finish`, `submit_cd.sh check`, przerzedzanie budżetowe, brak części → STOP; etap ORCH na atrapach. Nie testowane: condor, EOS, prawdziwy plik d23, MC v2.0.
