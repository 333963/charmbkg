#!/usr/bin/env python3
"""DLZ_DATA: sklejenie części apply (dlz_<próbka>_<g>.npz) → dlz_<próbka>.npz (układ slim/data_<run>.npz lokalnej analizy) + dlz_summary.json."""
import argparse
import glob
import json
import os

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--ngroups", required=True)
a = ap.parse_args()
ng = dict(x.split(":") for x in a.ngroups.split())
summ = {}
for s in a.samples.split():
    fs = sorted(glob.glob(os.path.join(a.parts, f"dlz_{s}_*.npz")), key=lambda f: int(f.rsplit("_", 1)[1][:-4]))
    if len(fs) != int(ng[s]):
        raise SystemExit(f"{s}: części {len(fs)} ≠ {ng[s]}")
    acc, off, fi, nf, nb, nv = {}, 0, [], 0, 0, 0
    hall = {"hall_any": 0, "hall_mb": 0}
    names = None
    for f in fs:
        z = np.load(f)
        for k in z.files:
            if k in ("fi", "feat_names", "n_files", "n_before_veto", "n_after_veto", "hall_any", "hall_mb"):
                continue
            acc.setdefault(k, []).append(z[k])
        fi.append(z["fi"] + off)
        off += len(z["m"])
        nf += int(z["n_files"]); nb += int(z["n_before_veto"]); nv += int(z["n_after_veto"])
        for k in hall:
            hall[k] = hall[k] + z[k]
        names = z["feat_names"]
    out = {k: np.concatenate(v) for k, v in acc.items()}
    out.update(fi=np.concatenate(fi), feat_names=names, n_files=nf, n_before_veto=nb, n_after_veto=nv, **hall)
    assert len(out["feat"]) == len(out["fi"])
    np.savez_compressed(os.path.join(a.out, f"dlz_{s}.npz"), **out)
    summ[s] = dict(n_files=nf, n_before_veto=nb, n_after_veto=nv, n_loose=int(len(out["m"])), n_feat=int(len(out["fi"])),
                   n_all_hist=int(np.sum(hall["hall_any"])), n_mb_hist=int(np.sum(hall["hall_mb"])))
    print(s, summ[s], flush=True)
json.dump(summ, open(os.path.join(a.out, "dlz_summary.json"), "w"), indent=1)
