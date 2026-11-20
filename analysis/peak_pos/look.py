import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import load, PTB
fig, ax = plt.subplots(2, 3, figsize=(15, 8))
for i, sp in enumerate(["D0", "Dplus"]):
    d = load(sp, "data")
    for j, (lo, hi) in enumerate(PTB):
        a = (d.pt >= lo) & (d.pt < hi)
        for q, c in [(1, "C0"), (-1, "C3")]:
            ax[i, j].hist(d.m[a & (d.q == q)], bins=np.arange(1.7, 2.2, 0.004), histtype="step", color=c, label=f"q={q}")
        ax[i, j].set_title(f"{sp} {lo}-{hi} n={a.sum()}"); ax[i, j].legend()
plt.savefig("out/look_data.png", dpi=70)
