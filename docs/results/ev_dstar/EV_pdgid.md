# EV — pełna lumi, prawda swapu closest_mass vs pdgid (v1.2), 2026-11-16

Wejście: `lxplus/ev_*.npz` (pełna lumi, v1.1, pobrane 14.11). `relabel_pdgid.py` dopina `Dmeson_hypo_swap` z wyciągu HYP v1.2
(`../HYP_d0/runs/lxplus/fd_nd.parquet`) po (plik, wiersz): 178 637 / 178 747 wierszy, masy zgodne 100 %. Wyniki:
`lxplus/ev_analysis_cm.{md,json}` (closest_mass) i `pdgid/ev_analysis.{md,json}` (pdgid). Progi WP bez zmian.

Relabel w parach D⁰←D* (fd): sygnał→swap 1 733, swap→sygnał 6 955; swapy MC z truth_ds 0/1 w oknie ±55 MeV 5–8 GeV: 1 598 → 83
(tag D* wybiera poprawną hipotezę; closest_mass oznaczał ~1,5 tys. prawdziwych D⁰ jako swap).

| bin | region | ε_dane WP95 | ε_MC WP95 (cm → pdgid) | SF95 (cm → pdgid) |
|---|---|---|---|---|
| 5–8 | ±36 MeV | 0,81 ± 0,08 | 0,95 → 0,95 | 0,85 ± 0,08 → 0,85 ± 0,08 |
| 5–8 | 36–55 MeV | 0,72 ± 0,31 | 0,93 → 0,94 | 0,77 ± 0,33 → 0,75 ± 0,23 |
| 5–8 | ±55 MeV | 0,79 ± 0,08 | 0,95 → 0,95 | 0,83 ± 0,08 → 0,83 ± 0,07 |
| 8–30 | ±55 MeV | 0,95 ± 0,04 | 0,95 → 0,95 | 1,00 ± 0,04 → 1,00 ± 0,03 |

Bin 3–5: fity pass/fail niestabilne (ε_dane > 1) — niemierzalny, jak w B2.

## Wniosek
- Prawda z v1.2 **nie tłumaczy** anomalii: ε_MC prawie bez zmian (+0,004), SF95 5–8 GeV = 0,83 ± 0,07 w obu wariantach.
- Na pełnej lumi anomalia jest mniejsza niż na podpróbce (0,74 → 0,79) i nie siedzi tylko w pierścieniu: rdzeń ±36 MeV też ma SF 0,85 ± 0,08.
  To ≈ 2,3σ od 1 w jednym binie pT (8–30: 1,00 ± 0,03). Hipoteza „szersze ogony m(D⁰) w danych” nie jest potwierdzona.
- Zalecenie: zamknąć jako systematykę ε_S w 5–8 GeV (|1 − SF| ≈ 0,17 przy WP95; przy WP70 SF 0,86–0,90 ± 0,09–0,12),
  zgodną z SF_rel B2 (0,97 ± 0,05 w 8–30). Bez dalszych badań, chyba że pojawi się więcej danych.
