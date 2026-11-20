#!/usr/bin/env python3
"""MCTMPL merge: sumuje histogramy <parts>/mctmpl_<próbka>_<g>.npz → <out>/mctmpl_<próbka>.npz + mctmpl_summary.json.
Kontrola: komplet grup, każdy plik h5 dokładnie raz, liczba plików = --nfiles, zgodne progi/klucze/kształty. Użycie:
  mctmpl_merge.py --parts DIR --out DIR --samples "a b" --ngroups "a:10 b:3" [--nfiles "a:97 b:24"]"""
import argparse
import json
import os
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--ngroups", required=True)
ap.add_argument("--nfiles", default="")
a = ap.parse_args()


def kv(s):
    return dict((x.split(":")[0], int(x.split(":")[1])) for x in s.replace(",", " ").split())


ng, nf = kv(a.ngroups), kv(a.nfiles)
S, err = {}, []
for s in a.samples.replace(",", " ").split():
    ps = [os.path.join(a.parts, f"mctmpl_{s}_{g}.npz") for g in range(ng[s])]
    miss = [os.path.basename(p) for p in ps if not os.path.exists(p)]
    if miss:
        err.append(f"{s}: brak części {miss}")
        continue
    H, files, m0, cnt = {}, [], None, {}
    for p in ps:
        z = np.load(p, allow_pickle=False)
        m = json.loads(str(z["meta"]))
        if m0 is None:
            m0 = m
        elif any(m[k] != m0[k] for k in ("keys", "thresholds", "spec", "axes", "lxy", "mctmpl")):
            err.append(f"{s}: niezgodne meta w {os.path.basename(p)}")
        for k in z.files:
            if k != "meta":
                H[k] = H[k] + z[k] if k in H else z[k].astype(np.int64)
        files += m["files"]
        for k, v in m["counts"].items():
            cnt[k] = cnt.get(k, 0) + v
    if len(set(files)) != len(files):
        err.append(f"{s}: pliki powtórzone")
    if s in nf and nf[s] != len(files):
        err.append(f"{s}: plików {len(files)} ≠ oczekiwane {nf[s]}")
    m0.update(files=files, counts=cnt, n_groups=len(ps), sample=s)
    o = os.path.join(a.out, f"mctmpl_{s}.npz")
    np.savez_compressed(o, meta=np.array(json.dumps(m0)), **H)
    S[s] = {"n_files": len(files), "bytes": os.path.getsize(o), "counts": cnt, "entries": {k: int(v.sum()) for k, v in H.items()}}
    print(s, cnt, f"{S[s]['bytes'] / 1e6:.2f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump({"samples": S, "bytes_total": sum(v["bytes"] for v in S.values())}, open(os.path.join(a.out, "mctmpl_summary.json"), "w"), indent=1)
