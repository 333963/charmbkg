#!/usr/bin/env python3
"""s0: cache + kontrole definicji tagu i zakresów."""
import numpy as np, common as C
for r in C.RUNS:
    d = C.load(r, force=True)
    n = len(d["lb"]); print(r, "zdarzeń", n, "GRL", int(d["grl"].sum()), "pt_min per gatunek", d["pt_min"], "D0 Lxy min", d["d0_lxy_min"])
    for s in C.SIDES:
        ref = (d[f"np_{s}"] == 1) & (d[f"nw_{s}"] == 1)
        inw = (d[f"xi1_{s}"] >= 0.03) & (d[f"xi1_{s}"] <= 0.10)
        print("  strona", s, "tag", int(d[f"tag_{s}"].sum()), "np==1&nw==1", int(ref.sum()), "zgodne", bool((ref == d[f"tag_{s}"]).all()),
              "tag bez ξ1 w oknie", int((d[f"tag_{s}"] & ~inw).sum()))
    T0, T1, T2 = C.regions(d)
    print("  T0,T1,T2:", int(T0.sum()), int(T1.sum()), int(T2.sum()), " bez GRL tag&npu0:", int(((d['tag_A'] | d['tag_C']) & (d['npu'] == 0)).sum()),
          " dwustronne w T2:", int((T2 & d['tag_A'] & d['tag_C']).sum()))
    print("  mu: GRL śr", float(d["mu"][d["grl"]].mean()), "poza GRL śr", float(d["mu"][~d["grl"]].mean()) if (~d["grl"]).any() else None)
