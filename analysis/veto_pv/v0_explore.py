#!/usr/bin/env python3
"""Rozpoznanie: rozkład nPileupVertices, actual_mu, LB, GRL w plikach B_trig (tylko odczyt)."""
import numpy as np, json
R = "/home/k/Documents/charm/charm_work/B_trig/runs/lxplus/"
for s in ["data_428770", "data_427929", "fd_nd", "mb_nd"]:
    z = np.load(R + f"bt_{s}.npz", allow_pickle=True)
    npu, mu, lb, grl = z["ev_npu"], z["ev_actual_mu"], z["ev_lb"], z["ev_grl"]
    print("==", s, "ev", len(npu), "grl", grl.mean().round(4))
    print("  npu counts", dict(zip(*[x.tolist() for x in np.unique(npu, return_counts=True)])))
    print("  mu pct [0,1,5,25,50,75,95,99,100]", np.percentile(mu, [0, 1, 5, 25, 50, 75, 95, 99, 100]).round(4), "mean", mu.mean().round(4), "nuniq", len(np.unique(mu)))
    print("  lb range", lb.min(), lb.max(), "nLB", len(np.unique(lb)))
    print("  pairs", len(z["m"]), "pt", z["pt"].min(), z["pt"].max(), "Lxy min", z["Lxy"].min(), "dm", z["dm"].min(), z["dm"].max(), "m", z["m"].min(), z["m"].max(), "score<0.05", (z["score"] < .05).mean().round(3))
    print("  all", len(z["all_m"]), "pt", z["all_pt"].min(), z["all_pt"].max(), "m", z["all_m"].min(), z["all_m"].max())
    print("  truth", dict(zip(*[x.tolist() for x in np.unique(z["truth"], return_counts=True)])), "truth_ds", dict(zip(*[x.tolist() for x in np.unique(z["truth_ds"], return_counts=True)])))
    m = json.loads(str(z["meta"])) ; print("  meta keys", list(m.keys()), str(m)[:600])
for s in ["data_428770", "data_427929"]:
    z = np.load(R + f"trg_{s}.npz", allow_pickle=True)
    npu, mu, grl = z["npu"], z["actual_mu"], z["grl"]
    print("== trg", s, len(npu), "grl", grl.mean().round(4), "P(npu=0)", (npu == 0).mean().round(5), "| GRL:", (npu[grl == 1] == 0).mean().round(5))
    print("  mu pct", np.percentile(mu, [0, 1, 5, 25, 50, 75, 95, 99, 100]).round(4))
    print("  keys nds", [k for k in z.files if k.startswith("nds")], "sum sr", {k: int(z[k].sum()) for k in z.files if k.startswith("nds")})
    print(str(json.loads(str(z["meta"])))[-800:])
