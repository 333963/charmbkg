# B1 — SALAD v1.0 + closure adwersarialny (2026-11-04)

Katalog `~/Documents/charm/charm_work/B1_salad_v1/`. Dokument etapu: `charm-B1-salad-v1.md` (kopia w tym katalogu).

- Addytywnie w `AA_merge/src`: `charmtag/bkggen/salad_v1.py`, `config/bkggen/salad_v1.yaml`, `tests/test_salad_v1.py` (nowe),
  `charmtag/bkggen/api.py` (+4 linie na końcu, `api_hook.diff`). Testy: `pytest tests/test_aa_merge.py tests/test_bkggen.py tests/test_salad_v1.py` → 32 passed (stan AA po wpięciu S8/DD2 z 13:02, + test_s8_salad).
- `b1_prep.py` — dp_*.npz (S6) → `sub/{mcbkg,sig,data}.npz` (cechy H, m, pt, truth, parity; na maszynie roboczej, nie kopiowane: 260 MB).
- `b1_study.py --mode smooth|cat_global|cat_local|cat_step|cuts_local|cat_local+smooth --kind hgb|sk_mlp [--leak 0.02]` — closure; `--stage report`.
- `b1_data.py [hgb|sk_mlp]` — bieg produkcyjny na danych → `out/salad_w_Dplus_data_p0_<kind>_v1.npz` (wagi dla MC mb_p0+p1 truth≥8
  w kolejności `b1_prep` = kolejność dp_mb_p0 ⧺ dp_mb_p1 po wecie D*/zakresie/truth≥2; kolumna `row` = indeks w tej ramce).
- `out/closure_*.json`, `out/b1_closure_Dplus.md`, `out/wSR_*.npz` (wagi SR do testu wycieku), `logs/`.
Uruchomienie na maszynie roboczej (2 CPU): closure ~10 min/tryb; `CHARMTAG_SRC` wskazuje `AA_merge/src`.
