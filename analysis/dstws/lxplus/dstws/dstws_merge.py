#!/usr/bin/env python3
"""DSTWS merge: parts/dw_<próbka>_p<k>.npz → dw_<próbka>.npz (konkatenacja; c_d0 przesunięte o liczbę D0 wcześniejszych części; histogramy
sumowane) + dw_summary.json. Kontrola: komplet części; każdy plik h5 dokładnie raz i liczba plików = n_files_total; zero niezgodności
eventNumber; c_d0 w zakresie.
Użycie: dstws_merge.py --parts DIR --out WD --samples "a b" --nparts "a:12 b:3" [--partial]"""
import argparse
import json
import os
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True); ap.add_argument("--out", required=True); ap.add_argument("--samples", required=True)
ap.add_argument("--nparts", required=True)
ap.add_argument("--partial", action="store_true", help="TEST: tylko istniejące części, bez kontroli kompletu (nie używać w produkcji)")
a = ap.parse_args()
NP = dict((kv.split(":")[0], int(kv.split(":")[1])) for kv in a.nparts.replace(",", " ").split())
summ, err = {}, []
for s in a.samples.replace(",", " ").split():
    fs = [os.path.join(a.parts, f"dw_{s}_p{k}.npz") for k in range(NP[s])]
    if a.partial:
        fs = [f for f in fs if os.path.exists(f)]
        if not fs:
            continue
    miss = [os.path.basename(f) for f in fs if not os.path.exists(f)]
    if miss:
        err.append(f"{s}: brak części {miss}"); continue
    cols, metas, H = {}, [], {}
    nd = 0
    for f in fs:
        z = np.load(f, allow_pickle=False)
        m = json.loads(str(z["meta"])); metas.append(m)
        if m["nparts"] != NP[s] and not a.partial:
            err.append(f"{s}: {os.path.basename(f)} liczone z nparts={m['nparts']} ≠ {NP[s]}")
        for c in z.files:
            if c == "meta":
                continue
            v = z[c]
            if c.startswith("H_"):
                H[c] = v if c not in H else H[c] + v
                continue
            if c == "c_d0":
                if len(v) and (v.min() < 0 or v.max() >= m["n_d0"]):
                    err.append(f"{s}: {os.path.basename(f)} c_d0 poza zakresem")
                v = v + nd
            cols.setdefault(c, []).append(v)
        nd += m["n_d0"]
    files = [x for m in metas for x in m["files"]]
    if a.partial:
        print(f"UWAGA: {s} --partial: {len(fs)} części, {len(files)} plików — wynik TESTOWY")
    elif len(set(files)) != len(files) or any(m["n_files_total"] != len(files) for m in metas):
        err.append(f"{s}: pliki {len(files)} (unikalnych {len(set(files))}), części widziały {sorted(set(m['n_files_total'] for m in metas))}")
    out = {c: np.concatenate(v) for c, v in cols.items()}
    out.update(H)
    nmis = sum(m["n_event_mismatch"] for m in metas)
    if nmis:
        err.append(f"{s}: niezgodności eventNumber {nmis}")
    pf = {k: v for m in metas for k, v in m["per_file"].items()}
    sec = sorted(v.get("sec", 0.0) for v in pf.values())
    typ = out["c_typ"]
    info = {"sample": s, "n_files": len(files), "nparts": NP[s], "partial": bool(a.partial), "n_d0_rows": int(len(out["d_pt"])),
            "n_d0_window": int(sum(v["n_d0_win"] for v in pf.values())), "n_comb": int(len(typ)),
            "n_by_type": {t: int((typ == i).sum()) for i, t in enumerate(metas[0]["types"])},
            "n_nominal_pairs": int(len(out["n_dm"])), "nominal_found_rs": float(np.mean(out["n_found"] == 1)) if len(out["n_dm"]) else None,
            "frac_daughters_matched": float(np.mean(out["d_nmatch"] == 2)) if len(out["d_pt"]) else None,
            "sec_per_file_median": sec[len(sec) // 2] if sec else None, "sec_per_file_max": sec[-1] if sec else None,
            "cpu_seconds": sum(m["seconds"] for m in metas)}
    keep = ("m0", "m_win", "nom_win", "lxy_min", "dm_max", "row_dm", "dm_max_nom", "pt_bins", "sc_n", "dm_edges", "nrot", "nmix", "match_dr",
            "match_rel", "types", "thresholds", "trig_bits", "flag_bits")
    out["meta"] = np.array(json.dumps(dict(info, files=files, **{k: metas[0][k] for k in keep},
                                           note="DSTWS: D0 (selekcja par B_trig, score H) × tor z Tracks: RS/WS/ROT/MIX; c_d0 = indeks w d_*")))
    p = os.path.join(a.out, f"dw_{s}.npz")
    np.savez_compressed(p, **out)
    info["bytes"] = os.path.getsize(p)
    summ[s] = info
    print(s, info["n_comb"], "kombinacji", f"{info['bytes'] / 1e6:.1f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump(summ, open(os.path.join(a.out, "dw_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
