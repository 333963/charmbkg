import sys, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/k/Documents/charm/charm_work/B_trig/ISO_origin")
import iso_common as I
print("SIG", I.SIG, "import ok")
R = "/home/k/Documents/charm/charm_work/B_trig/runs/lxplus/"
z = np.load(R + "bt_data_427929.npz")
print([k for k in z.files if k.startswith("all_")])
mu, npu, grl, lb = z["ev_actual_mu"], z["ev_npu"], z["ev_grl"], z["ev_lb"]
print("427929 mu==0:", (mu == 0).sum(), "GRL mu pct", np.percentile(mu[grl == 1], [0, 1, 50, 99, 99.9, 100]), "nonGRL mu pct", np.percentile(mu[grl == 0], [0, 1, 25, 50, 75, 99, 100]))
h = (grl == 0) & (mu > 5)
print("nonGRL high-mu: n", h.sum(), "<npu>", npu[h].mean(), "<mu>", mu[h].mean(), "ratio", npu[h].mean() / mu[h].mean(), "LB", np.unique(lb[h])[:5], np.unique(lb[h])[-5:])
l = (grl == 0) & (mu <= 5); print("nonGRL low-mu n", l.sum(), "P(npu>0)", (npu[l] > 0).mean(), "LB", np.unique(lb[l])[:10], "mu med", np.median(mu[l]))
g = grl == 1; print("GRL LB", lb[g].min(), lb[g].max(), "P(npu>0)", (npu[g] > 0).mean(), "npu counts", np.bincount(npu[g])[:8])
z = np.load(R + "bt_data_428770.npz")
mu, npu, lb = z["ev_actual_mu"], z["ev_npu"], z["ev_lb"]
print("428770 mu==0:", (mu == 0).sum(), "LB uniq", np.unique(lb)[:30], np.unique(lb)[-5:])
for a, b in [(0, 100), (100, 250), (250, 290), (290, 400)]:
    s = (lb >= a) & (lb < b); print(" LB", a, b, s.sum(), "mu", mu[s].mean() if s.any() else None, "P(npu>0)", (npu[s] > 0).mean() if s.any() else None)
# ev vs all
ae = np.unique(z["all_event"]); print("ev", len(z["ev_event"]), "uniq all events", len(ae), "pairs events", len(np.unique(z["event"])))
