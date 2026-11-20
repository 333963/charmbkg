"""PEAK_pos: kontrola D⁺ |η| — fit z k ustalonym na 1 (σ z MC) i wspólnym k z fitu D+D̄."""
import numpy as np
from common import load, M0
from fitlib import fit_mc, fit_data
D, S = load("Dplus", "data"), load("Dplus", "mc")
for nm, f in [("abseta<1", lambda L: np.abs(L.eta) < 1), ("abseta>1", lambda L: np.abs(L.eta) >= 1),
              ("eta<-1", lambda L: L.eta < -1), ("eta>1", lambda L: L.eta >= 1)]:
    sb, db = np.asarray(f(S)), np.asarray(f(D))
    shp = fit_mc(S.m.to_numpy()[sb], 1.74, 2.04)
    kall = fit_data(D.m.to_numpy()[db], 1.74, 2.04, shp, mu0=M0["Dplus"])["k"]
    out = []
    for kf, lab in ((1.0, "k=1"), (kall, f"k=k(D+D̄)={kall:.2f}")):
        r = {}
        for q in (1, -1):
            mc = fit_mc(S.m.to_numpy()[sb & (S.q == q).to_numpy()], 1.74, 2.04)
            x = fit_data(D.m.to_numpy()[db & (D.q == q).to_numpy()], 1.74, 2.04, shp, mu0=M0["Dplus"], k_fixed=kf)
            r[q] = ((x["mu"] - mc["mu"]) * 1e3, np.hypot(x["e_mu"], mc["e_mu"]) * 1e3, x["Ns"])
        a = r[1][0] - r[-1][0]; e = np.hypot(r[1][1], r[-1][1])
        out.append(f"{lab}: Δμ D {r[1][0]:+.2f}±{r[1][1]:.2f} (N={r[1][2]:.0f}), D̄ {r[-1][0]:+.2f}±{r[-1][1]:.2f} (N={r[-1][2]:.0f}), (D−D̄)−MC {a:+.2f}±{e:.2f}")
    print(nm, " | ".join(out))
