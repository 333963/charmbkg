"""T0: odtworzenie liczb etapu 4 + podstawowe liczby (kandydaci/zdarzenie, blizniaki)."""
import numpy as np, json
from vlib import *
F = load("fd_d23")
print("fd_nd 2-3: sygnal poprawny", (F["cat"] == 0).sum(), "zamieniony", (F["cat"] == 1).sum(), "inne", (F["cat"] >= 2).sum(), "wagi", np.unique(F["w"])[:5])
for run in ("427929", "428770"):
    D = load(f"data_{run}_d23")
    print(run, "ntot", D["ntot"], "z blizniakiem w tabeli", D["ntw"], "wiersze", len(D["m"]))
    for j in (0, 1):
        for wp in (0.3, 0.2, 0.1):
            th = mc_wp(F, j, wp); q = templ(F, j, th)
            s = D["lg"][:, j] >= th
            f = fit(D["m"][s], q)
            ff = fit(D["m"][s], q, free=True)
            print(f"  {NM[j]} WP{int(wp*100)} thr {th:.4f} eps {q['eps']:.4f} frefl {q['frefl']:.3f} nsel {s.sum()} N {f['N']:.1f} ± {f['err']:.1f} Z {f['Z']:.2f} chi2 {f['chi2']:.1f}/{f['ndf']} y {f['N']/q['eps']:.0f}"
                  f" | wolny: N {ff['N']:.0f} ± {ff['err']:.0f} dmu {1e3*ff['dmu']:.2f} ± {1e3*ff['dmu_err']:.2f} MeV sc {ff['sc']:.3f} ± {ff['sc_err']:.3f}")
