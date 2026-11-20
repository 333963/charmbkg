#!/usr/bin/env python3
"""B_trig merge: parts/bt_<próbka>_<g>.npz → bt_<próbka>.npz (konkatenacja; file_idx/all_file_idx/ev_file_idx + g·1000) + bt_summary.json.
Kontrola: każda próbka z --need obecna; tabela ev_* bez duplikatów (run, event) w danych (bramka G1)."""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument("--parts", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--need", default="")
a = ap.parse_args()
by = {}
for f in sorted(glob.glob(os.path.join(a.parts, "bt_*_*.npz"))):
    m = re.match(r"bt_(\w+?)_(\d+)\.npz$", os.path.basename(f))
    if m:
        by.setdefault(m.group(1), []).append((int(m.group(2)), f))
miss = [s for s in a.need.split() if s not in by]
if miss:
    sys.exit(f"brak części: {miss}")
summ = {}
for s, lst in by.items():
    cols, metas = {}, []
    for g, f in sorted(lst):
        z = np.load(f, allow_pickle=False)
        metas.append({"group": g, **json.loads(str(z["meta"]))})
        for k in z.files:
            if k == "meta":
                continue
            v = z[k]
            if k in ("file_idx", "all_file_idx", "ev_file_idx"):
                v = v + g * 1000
            cols.setdefault(k, []).append(v)
    keys = set.intersection(*[set(np.load(f).files) for _, f in lst]) - {"meta"}
    out = {k: np.concatenate(v) for k, v in cols.items() if k in keys}
    dup = 0
    if s.startswith("data_") and "ev_run" in out:
        k = (out["ev_run"].astype(np.int64) << 43) | out["ev_event"].astype(np.int64)
        dup = int(len(k) - len(np.unique(k)))
    m0 = metas[0]
    out["meta"] = np.array(json.dumps({"sample": s, "groups": metas, "features": m0["features"], "thresholds": m0["thresholds"], "m0": m0["m0"],
                                       "event_columns": m0.get("event_columns"), "n_pairs": int(len(out.get("m", []))),
                                       "n_all": int(len(out.get("all_m", []))), "ev_dup_run_event": dup}))
    np.savez_compressed(os.path.join(a.out, f"bt_{s}.npz"), **out)
    summ[s] = {"n_groups": len(lst), "n_pairs": int(len(out.get("m", []))), "n_all": int(len(out.get("all_m", []))),
               "n_events": int(len(out.get("ev_run", []))), "ev_dup_run_event": dup, "event_columns": m0.get("event_columns")}
json.dump(summ, open(os.path.join(a.out, "bt_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
