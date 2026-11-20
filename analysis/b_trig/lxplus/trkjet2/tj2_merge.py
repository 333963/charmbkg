#!/usr/bin/env python3
"""TRKJET2 merge: parts/tj2_<próbka>_p<k>.npz → tj2_<próbka>.npz (konkatenacja, posortowane po bt_idx) + tj2_summary.json.
Kontrola: komplet części; każdy plik h5 dokładnie raz i liczba plików = n_files_total; każda para bt dokładnie raz (bt_idx = 0..N−1, N z bt npz);
zero niezgodności eventNumber. Kolumny x_tp_/x_tbp_ (bity dżetowe z ROOT): część wspólna części (MC może ich nie mieć).
Użycie: tj2_merge.py --parts DIR --out WD --samples "a b" --nparts "a:12 b:3" [--bt BT_DIR]"""
import argparse
import json
import os
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True); ap.add_argument("--out", required=True); ap.add_argument("--samples", required=True)
ap.add_argument("--nparts", required=True); ap.add_argument("--bt", default="")
ap.add_argument("--partial", action="store_true", help="TEST: tylko istniejące części, bez kontroli kompletu (nie używać w produkcji)")
a = ap.parse_args()
NP = dict((kv.split(":")[0], int(kv.split(":")[1])) for kv in a.nparts.replace(",", " ").split())
summ, err = {}, []
for s in a.samples.replace(",", " ").split():
    fs = [os.path.join(a.parts, f"tj2_{s}_p{k}.npz") for k in range(NP[s])]
    if a.partial:
        fs = [f for f in fs if os.path.exists(f)]
    miss = [os.path.basename(f) for f in fs if not os.path.exists(f)]
    if miss:
        err.append(f"{s}: brak części {miss}"); continue
    cols, metas, keys = {}, [], None
    for f in fs:
        z = np.load(f, allow_pickle=False)
        m = json.loads(str(z["meta"])); metas.append(m)
        if m["nparts"] != NP[s] and not a.partial:
            err.append(f"{s}: {os.path.basename(f)} liczone z nparts={m['nparts']} ≠ {NP[s]}")
        k = set(z.files) - {"meta"}
        keys = k if keys is None else keys & k
        for c in k:
            cols.setdefault(c, []).append(z[c])
    files = [x for m in metas for x in m["files"]]
    if a.partial:
        print(f"UWAGA: {s} --partial: {len(fs)} części, {len(files)} plików — wynik TESTOWY")
    elif len(set(files)) != len(files) or any(m["n_files_total"] != len(files) for m in metas):
        err.append(f"{s}: pliki {len(files)} (unikalnych {len(set(files))}), części widziały {sorted(set(m['n_files_total'] for m in metas))}")
    out = {c: np.concatenate(v) for c, v in cols.items() if c in keys}
    o = np.argsort(out["bt_idx"], kind="stable")
    out = {c: v[o] for c, v in out.items()}
    n = len(out["bt_idx"])
    nbt = None
    if a.bt and os.path.exists(os.path.join(a.bt, f"bt_{s}.npz")):
        nbt = int(json.loads(str(np.load(os.path.join(a.bt, f"bt_{s}.npz"))["meta"]))["n_pairs"])
    if not a.partial and (not np.array_equal(out["bt_idx"], np.arange(n)) or (nbt is not None and nbt != n)):
        err.append(f"{s}: pary {n}, bt {nbt}, bt_idx niekompletny/powtórzony")
    nmis = sum(m["n_event_mismatch"] for m in metas)
    if nmis:
        err.append(f"{s}: niezgodności eventNumber {nmis}")
    sec = sorted(v for m in metas for v in m["sec_per_file"].values())
    ok = (out["nmatch"] + out["nd_lo"] == 2) & (out["same_entry"] == 1)
    d = out["v_s40_p2"][ok] - out["ptcone40"][ok]
    info = {"sample": s, "n_pairs": n, "n_files": len(files), "nparts": NP[s], "columns": sorted(out), "xtrig": sorted(c for c in out if c.startswith("x_t")),
            "frac_same_entry": float(np.mean(out["same_entry"] == 1)) if n else None,
            "frac_daughters_ok": float(np.mean(out["nmatch"] + out["nd_lo"] == 2)) if n else None,
            "frac_pis_match": float(np.mean(out["pis_match"] == 1)) if n else None,
            "frac_rnd_ok": float(np.mean(out["r_ok"] == 1)) if n else None, "xi": sorted(c for c in out if c.startswith("xi_cal")),
            "frac_xi_nan": float(np.mean(~np.isfinite(out["xi_cal_A"]))) if n and "xi_cal_A" in out else None,
            "ptcone40_own_vs_branch": {"frac_abs_lt_1MeV": float(np.mean(np.abs(d) < 1e-3)), "frac_abs_lt_100MeV": float(np.mean(np.abs(d) < 0.1)),
                                       "median": float(np.median(d)), "rms": float(np.sqrt(np.mean(d ** 2)))} if ok.any() else None,
            "n_pis_mismatch": sum(m["n_pis_mismatch"] for m in metas), "n_d0pt_mismatch": sum(m["n_d0pt_mismatch"] for m in metas),
            "sec_per_file_median": sec[len(sec) // 2] if sec else None, "cpu_seconds": sum(m["seconds"] for m in metas)}
    out["meta"] = np.array(json.dumps(dict(info, files=files, names=metas[0]["names"], jet=metas[0]["jet"], thresholds=metas[0]["thresholds"],
                                           m0=metas[0]["m0"], note="TRKJET2: pary z bt_<próbka>.npz (b_trig) + dżety torowe anty-kT R=0,4 + stożki prostopadłe/losowy, obszary Δφ, ξ_trk, ξ_cal, bity TBP "
                                           "(własna implementacja), stożki, półkule; bt_idx = indeks pary w bt npz")))
    p = os.path.join(a.out, f"tj2_{s}.npz")
    np.savez_compressed(p, **out)
    info["bytes"] = os.path.getsize(p); info.pop("columns")
    summ[s] = info
    print(s, n, "par", f"{info['bytes'] / 1e6:.1f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump(summ, open(os.path.join(a.out, "tj2_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
