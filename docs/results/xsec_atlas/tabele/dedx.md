# dE/dx: przesunięcie MC ∓0,023 (charmpy `training.yaml`: `dedx_mc_shift: -0.023`) → wydajność drzewa H przy stałym progu WP70

Sygnał prawdy (fd_nd + mb_nd, tabele INCAL `__sig_`), mianownik: preselekcja nominalna (L_xy, weto D*). Systematyka = max |ε(±)/ε − 1| (pełna zmiana; nasze MC nie ma przesunięcia).

| gatunek | bin | N sygnału MC | ε nominalna | max \|Δscore\| odtw. | ε(dE/dx − 0,023)/ε | ε(dE/dx + 0,023)/ε | systematyka ± [%] | ⟨dE/dx K⟩ |
|---|---|---|---|---|---|---|---|---|
| D0 | 3-5 | 129073 | 0.693 | 8.0e-04 | 0.9939 ± 0.0004 | 1.0051 ± 0.0004 | 0.6 | 0.932 |
| D0 | 5-8 | 68240 | 0.698 | 8.0e-04 | 0.9954 ± 0.0006 | 1.0047 ± 0.0006 | 0.5 | 0.914 |
| D0 | 8-30 | 29055 | 0.701 | 8.0e-04 | 0.9969 ± 0.0009 | 1.0016 ± 0.0009 | 0.3 | 0.910 |
| Dplus | 3-5 | 83913 | 0.695 | 3.9e-03 | 0.9952 ± 0.0005 | 1.0032 ± 0.0004 | 0.5 | 0.931 |
| Dplus | 5-8 | 50204 | 0.696 | 3.9e-03 | 0.9958 ± 0.0006 | 1.0009 ± 0.0006 | 0.4 | 0.896 |
| Dplus | 8-30 | 22783 | 0.696 | 3.9e-03 | 0.9986 ± 0.0010 | 1.0011 ± 0.0010 | 0.1 | 0.892 |
