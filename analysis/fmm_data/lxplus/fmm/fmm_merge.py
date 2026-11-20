#!/usr/bin/env python3
"""FMM merge: skleja <parts>/fmm_<próbka>_<g>.npz → <out>/fmm_<próbka>.npz (file_idx przenumerowany globalnie) + fmm_summary.json.
Kontrola: komplet grup (--ngroups "a:10 b:3"), każdy plik h5 dokładnie raz, liczba plików = --nfiles. Użycie:
  fmm_merge.py --parts DIR --out DIR --samples "a b" --ngroups "a:10 b:3" [--nfiles "a:97 b:24"]"""
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
    ps = [os.path.join(a.parts, f"fmm_{s}_{g}.npz") for g in range(ng[s])]
    miss = [os.path.basename(p) for p in ps if not os.path.exists(p)]
    if miss:
        err.append(f"{s}: brak części {miss}")
        continue
    cols, files, metas, rpf = {}, [], [], {}
    for p in ps:
        z = np.load(p, allow_pickle=False)
        m = json.loads(str(z["meta"]))
        metas.append(m)
        for k in z.files:
            if k == "meta":
                continue
            v = z[k]
            if k == "file_idx":
                v = (v.astype(np.int32) + len(files)).astype(np.int16)
            cols.setdefault(k, []).append(v)
        files += m["files"]
        rpf.update(m["rows_per_file"])
    if len(set(files)) != len(files):
        err.append(f"{s}: pliki powtórzone")
    if s in nf and nf[s] != len(files):
        err.append(f"{s}: plików {len(files)} ≠ oczekiwane {nf[s]}")
    bad = sorted(k for k, v in cols.items() if len(v) != len(ps))
    if bad:
        err.append(f"{s}: kolumny nie we wszystkich częściach: {bad}")
        continue
    out = {k: np.concatenate(v) for k, v in cols.items()}
    n = len(out["m"])
    m0 = dict(metas[0])
    m0.update(files=files, rows_per_file=rpf, n=n, n_groups=len(ps))
    out["meta"] = np.array(json.dumps(m0))
    o = os.path.join(a.out, f"fmm_{s}.npz")
    np.savez_compressed(o, **out)
    S[s] = {"n": n, "n_files": len(files), "bytes": os.path.getsize(o), "columns": sorted(cols), "frac": m0["frac"],
            "n_twin": int((out["n_grp"] == 2).sum()), "n_sr": int((np.abs(out["m"] - m0["m0"]) < 0.055).sum())}
    print(s, n, "wierszy", f"{S[s]['bytes'] / 1e6:.1f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump({"samples": S, "bytes_total": sum(v["bytes"] for v in S.values())}, open(os.path.join(a.out, "fmm_summary.json"), "w"), indent=1)
