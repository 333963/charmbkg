import sys, numpy as np, time
from common import *
for s in SAMPLES:
    t0=time.time(); d = load(s); Dn, T, ca = d["Dn"], d["T"], d["cands"]
    print(s, len(Dn), len(T), len(ca), f"{time.time()-t0:.0f}s")
    print(" pt quantiles", np.quantile(ca.pt, [0,0.001,0.5,0.999,1]), "frac pt<3", (ca.pt<3).mean(), "pt>30", (ca.pt>30).mean())
    for code in (0,3,1):
        c = ca[ca.species==code]; print("  species", code, len(c), "m", c.m.min(), c.m.max(), "dm", c.dm.min(), c.dm.max(), "ptmin", c.pt.min())
    for sd in "AC":
        n = Dn[f"np_{sd}"].to_numpy(); x1 = Dn[f"xi1_{sd}"].to_numpy(); xc = Dn[f"xical_{sd}"].to_numpy()
        print(" side", sd, "np counts", np.bincount(np.minimum(n,4)))
        h = np.histogram(x1[n==1], [-2,0,0.01,0.02,0.025,0.03,0.04,0.05,0.06,0.08,0.10,0.12,0.15,0.2,1])[0]; print("  xi1 (np==1) hist", h)
        pool = xc>0.14
        h = np.histogram(x1[(n==1)&pool], [-2,0,0.01,0.02,0.025,0.03,0.04,0.05,0.06,0.08,0.10,0.12,0.15,0.2,1])[0]; print("  xi1 (np==1, pool) hist", h)
        print("  xical quantiles", np.quantile(xc,[0.001,0.01,0.1,0.5,0.9,0.99]), "frac pool", pool.mean(), "hasD pool", (pool&Dn.hasD.to_numpy()).sum())
