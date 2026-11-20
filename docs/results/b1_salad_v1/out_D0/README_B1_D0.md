# B1-D0 — SALAD v1.0 dla D0: produkcja na danych + closure adwersarialny (2026-11-04)

Katalog `~/Documents/charm/charm_work/B1_salad_v1/out_D0/`. Dokument: `charm-B1-salad-v1.md` (pełna kopia dokumentu B1 z nową §6).
`AA_merge/src` NIE zmieniony (etap analityczny; pakiet staged jako `_in_src_d0.tgz`).

- `code/b1_prep_d0.py` — `S_bkggen/sub/{mb_p0..p3,fd_p0,fd_p1,data_p0,data_p1}.npz` → `sub/{mcbkg,sig,data}.npz` (cechy H, m, m_swap, pt, truth, parity, src;
  D0 decayType 10, 1,75–2,00 GeV, pT 3–30, weto refleksji |m_swap−m0|>25 MeV wszędzie; na maszynie roboczej, 230 MB, nie kopiowane).
- `code/b1_study_d0.py --mode smooth|cat_global|cat_local|cat_local_MM|cat_step|cuts_local|cat_local+smooth|swap_local|swap_plus|real_t27 --kind hgb|sk_mlp
  [--leak f] [--strength 0.8|1.6] [--swap_factor 1.5|0.5] [--mc_truth_min 8|2] [--seed n]`; `--stage report` → `b1_closure_D0.md`.
  Ta sama logika/metryki co `b1_study.py` (D⁺); okno wycieku 75–95 MeV; tryby swap_*/real_t27 i ε_B per kategoria — nowe (D0).
- `code/b1_data_d0.py [hgb|sk_mlp] [0|01]` — produkcja na danych → `salad_w_D0_data_p0_{hgb,sk_mlp}_v1.npz`, `salad_w_D0_data_p01_hgb_v1.npz`
  (wagi dla MC mb_p0+p1 truth≥8 w kolejności `b1_prep_d0` = mb_p0 ⧺ mb_p1 po zakresie/wecie, kolumna `row`).
- `code/b1_dataeff_d0.py [8|2]` — kontrola na danych (ε_B z fitu pol2 SB→okno po cięciu vs SALAD) → `data_eff_D0.json`, `data_eff_D0_mc2.json`.
- `closure_*.json` (27 biegów; `_seed1`, `_s1.6`, `_mc2`, `_x1.5/_x0.5` w nazwie), `wSR_*.npz`, `logs_b1_d0.tgz`, skrypty `run_*_d0.sh`.
Maszyna robocza, 2 CPU: prep 10 s, ~2,5 min/bieg closure, produkcja 1 min, data_eff 3 min.
