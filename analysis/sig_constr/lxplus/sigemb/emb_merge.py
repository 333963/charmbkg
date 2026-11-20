#!/usr/bin/env python3
"""SIGEMB merge: parts/emb_<próbka>_p<k>.npz → emb_<próbka>.npz + emb_summary.json. Kontrola: komplet części, każdy plik ROOT dokładnie raz.
Użycie: emb_merge.py --parts DIR --out WD --samples "a b" --nparts "a:12 b:5" [--partial]"""
import argparse, json, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True); ap.add_argument("--out", required=True); ap.add_argument("--samples", required=True); ap.add_argument("--nparts", required=True)
ap.add_argument("--partial", action="store_true", help="TEST: tylko istniejące części")
a = ap.parse_args()
NP = dict((kv.split(":")[0], int(kv.split(":")[1])) for kv in a.nparts.replace(",", " ").split())
summ, err = {}, []
for s in a.samples.replace(",", " ").split():
    fs = [os.path.join(a.parts, f"emb_{s}_p{k}.npz") for k in range(NP[s])]
    if a.partial:
        fs = [f for f in fs if os.path.exists(f)]
    miss = [os.path.basename(f) for f in fs if not os.path.exists(f)]
    if miss or not fs:
        err.append(f"{s}: brak części {miss}"); continue
    cols, metas, keys = {}, [], None
    for f in fs:
        z = np.load(f, allow_pickle=False)
        metas.append(json.loads(str(z["meta"])))
        k = set(z.files) - {"meta"}
        keys = k if keys is None else keys & k
        for c in k:
            cols.setdefault(c, []).append(z[c])
    files = [x for m in metas for x in m["files"]]
    if not a.partial and (len(set(files)) != len(files) or any(m["n_files_total"] != len(files) for m in metas)):
        err.append(f"{s}: pliki {len(files)} (unikalnych {len(set(files))}), oczekiwane {sorted(set(m['n_files_total'] for m in metas))}")
    out = {c: np.concatenate(v) for c, v in cols.items() if c in keys}
    n = len(out["var"])
    info = {"sample": s, "n_rows": n, "n_carriers": sum(m["n_carriers"] for m in metas), "n_entries": sum(m["n_entries"] for m in metas), "n_files": len(files),
            "prescale": metas[0]["prescale"], "frac_ok": {str(v): float(np.mean(out["ok"][out["var"] == v])) for v in range(3)},
            "eps_cone_lt2": {str(v): float(np.mean(out["v_s40_p2"][(out["var"] == v) & (out["ok"] == 1)] < 2.0)) for v in range(3)},
            "cpu_seconds": sum(m["seconds"] for m in metas)}
    out["meta"] = np.array(json.dumps(dict(info, files=files, names=metas[0]["names"], variants=metas[0]["variants"],
                                           note="SIGEMB: embedding (pT, η) D0 z MC w nośniki; wiersz = (nośnik, wariant kierunku); kolumny jądra TRKJET2")))
    p = os.path.join(a.out, f"emb_{s}.npz")
    np.savez_compressed(p, **out)
    info["bytes"] = os.path.getsize(p); summ[s] = info
    print(s, n, "wierszy", f"{info['bytes'] / 1e6:.1f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump(summ, open(os.path.join(a.out, "emb_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
