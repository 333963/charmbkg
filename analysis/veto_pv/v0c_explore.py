import numpy as np
R = "/home/k/Documents/charm/charm_work/B_trig/runs/lxplus/"
z = np.load(R + "trg_data_427929.npz")
npu, mu, grl, lb, ns, ev = z["npu"], z["actual_mu"], z["grl"], z["lb"], z["nSelBL"], z["event"]
g = (grl == 1) & (mu > 0)
print("GRL npu counts", np.bincount(npu[g])[:25].tolist())
for a, b in [(0, 1), (1, 2), (2, 3), (3, 6), (6, 12), (12, 40)]:
    m = g & (npu >= a) & (npu < b)
    print(f"npu {a}-{b}: n={m.sum()} nSelBL med={np.median(ns[m]):.0f} mean={ns[m].mean():.1f} mu med={np.median(mu[m]):.4f} mean={mu[m].mean():.4f}",
          "sptrk4", z["tp_HLT_mb_sptrk_pt4_L1MBTS_2"][m].mean().round(3), "J12", z["tp_HLT_noalg_L1J12"][m].mean().round(4), "mbts1", z["tp_HLT_mb_mbts_L1MBTS_1"][m].mean().round(3), "j20", z["tp_HLT_j20_L1MBTS_2"][m].mean().round(4))
# npu>=1 fraction in fine mu bins, all GRL events and nSelBL<10
ed = np.quantile(mu[g], np.linspace(0, 1, 13))
for a, b in zip(ed[:-1], ed[1:]):
    m = g & (mu >= a) & (mu < b); l = m & (ns < 10); h = m & (ns >= 20)
    print(f"mu {a:.4f}-{b:.4f} n={m.sum()} P(npu=1)={(npu[m]==1).mean()*100:.3f}% P(npu=2)={(npu[m]==2).mean()*100:.4f}% P(npu>=3)={(npu[m]>=3).mean()*100:.3f}% | nSel<10: P(npu>=1)={(npu[l]>=1).mean()*100:.3f}% | nSel>=20: P(npu=1,2)={((npu[h]>=1)&(npu[h]<3)).mean()*100:.3f}%")
# time dependence of the tail
for a, b in [(119, 170), (170, 220), (220, 270), (270, 325)]:
    m = g & (lb >= a) & (lb < b); print("LB", a, b, "P(npu>=3) %", (npu[m] >= 3).mean() * 100, "P(npu=1) %", (npu[m] == 1).mean() * 100, "mu", mu[m].mean())
# event number adjacency of tail events (same crossing population?)
t = g & (npu >= 6); print("tail mu uniq", len(np.unique(mu[t])), "of", len(np.unique(mu[g])), "; top mu values", np.unique(mu[t], return_counts=True)[1].max())
z2 = np.load(R + "trg_data_428770.npz"); npu2, mu2, ns2 = z2["npu"], z2["actual_mu"], z2["nSelBL"]
print("428770 npu counts", np.bincount(npu2)[:14].tolist())
for a, b in [(0, 1), (1, 2), (2, 3), (3, 20)]:
    m = (npu2 >= a) & (npu2 < b); print(f"428770 npu {a}-{b}: n={m.sum()} nSelBL mean={ns2[m].mean():.1f} mu mean={mu2[m].mean():.4f}")
f = np.load(R + "bt_fd_nd.npz")
print("fd all_event range", f["all_event"].min(), f["all_event"].max(), "ev_event", f["ev_event"].min(), f["ev_event"].max(), "pairs event", f["event"].min(), f["event"].max(),
      "all_file_idx uniq", np.unique(f["all_file_idx"])[:5], "ev_file_idx", np.unique(f["ev_file_idx"])[:5], "file_idx", np.unique(f["file_idx"])[:5])
print("all_w uniq", np.unique(f["all_w"])[:5], "weight uniq", np.unique(f["weight"])[:5])
