#!/usr/bin/env python3
"""CAND_D23 finish: części (parts/dp_<s>_<g>.parquet, hist_<s>_<g>.npz, part_<s>_<g>.json) → dp_<s>.parquet, hist_<s>.npz, summary.json.
Kontrole (błąd = brak DONE): komplet części (--ngroups), pT D⁺ ∈ [lo, hi) w każdym wierszu, brak NaN w score, > 0 wierszy na próbkę,
suma plików wynikowych ≤ --hard-gb. Budżet --max-gb: gdy suma części parquet go przekracza, dodatkowe przerzedzenie po `ps_u`
(spójne z preskalą ekstrakcji; w_ps przeliczane): najpierw ps_tier 0/1 czynnikiem g01 ≥ --g-min, potem ps_tier 2 czynnikiem g2;
ps_tier 3 (sygnał MC) nietknięty. Histogramy (pełna statystyka) nie są przerzedzane."""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--ngroups", required=True, help="próbka:liczba_części ...")
ap.add_argument("--max-gb", type=float, default=1.8)
ap.add_argument("--hard-gb", type=float, default=2.0)
ap.add_argument("--g-min", type=float, default=0.05)
ap.add_argument("--thr", default="", help="thr.json treningu dlz_data (progi WP binu 3–5 GeV — tylko odniesienie)")
a = ap.parse_args()
ng = dict(x.split(":") for x in a.ngroups.replace(",", " ").split())
samples = a.samples.replace(",", " ").split()
gi = lambda f: int(f.rsplit("_", 1)[1].split(".")[0])  # noqa: E731
os.makedirs(a.out, exist_ok=True)
P = {}
for s in samples:
    js = sorted(glob.glob(os.path.join(a.parts, f"part_{s}_*.json")), key=gi)
    if len(js) != int(ng[s]) or [gi(f) for f in js] != list(range(int(ng[s]))):
        sys.exit(f"STOP {s}: części {len(js)} ≠ {ng[s]}")
    P[s] = [json.load(open(f)) for f in js]
# ---- budżet
rows = {t: 0 for t in range(4)}
byts = 0
for s in samples:
    for j in P[s]:
        byts += j["dp_bytes"]
        for t in range(4):
            rows[t] += j["counts"][f"dp_tier{t}"]
ntot = sum(rows.values())
bpr = byts / max(ntot, 1)
budget = a.max_gb * 1e9
g01 = g2 = 1.0
if byts > budget:
    b3, b2, b01 = bpr * rows[3], bpr * rows[2], bpr * (rows[0] + rows[1])
    g01 = float(np.clip((budget - b3 - b2) / max(b01, 1), a.g_min, 1.0))
    if b3 + b2 + g01 * b01 > budget:
        g2 = float(np.clip((budget - b3 - g01 * b01) / max(b2, 1), 0.01, 1.0))
    print(f"budżet: części {byts/1e9:.2f} GB > {a.max_gb} GB → g01={g01:.4f} g2={g2:.4f} (wiersze tier0..3: {rows})", flush=True)
G = np.array([g01, g01, g2, 1.0])
S, tot = {}, 0
for s in samples:
    lo, hi = P[s][0]["pt_range"]
    cnt, hist, W, n, nnan, files, fidx, miss = {}, {}, None, 0, 0, [], [], set()
    ptmin, ptmax = np.inf, -np.inf
    tiers = np.zeros(4, np.int64)
    out = os.path.join(a.out, f"dp_{s}.parquet")
    for g, j in enumerate(P[s]):
        for k, v in j["counts"].items():
            cnt[k] = cnt.get(k, 0) + int(v)
        files += j["files"]; fidx += j["file_idx"]; miss |= set(j["missing_columns"])
        z = np.load(os.path.join(a.parts, f"hist_{s}_{g}.npz"))
        for k in z.files:
            hist[k] = hist[k] + z[k] if k in hist and k.endswith("_h") else z[k]
        f = os.path.join(a.parts, f"dp_{s}_{g}.parquet")
        if not os.path.exists(f):
            if j["counts"]["dp_kept"]:
                sys.exit(f"STOP: brak {f} (dp_kept={j['counts']['dp_kept']})")
            continue
        t = pq.read_table(f)
        if len(t) != j["counts"]["dp_kept"]:
            sys.exit(f"STOP: {f}: {len(t)} wierszy ≠ dp_kept {j['counts']['dp_kept']}")
        if g01 < 1 or g2 < 1:
            tier = t["ps_tier"].to_numpy()
            w = t["w_ps"].to_numpy()
            gg = G[tier]
            keep = t["ps_u"].to_numpy() < np.round(gg / w * 1_000_000)
            t = t.set_column(t.schema.get_field_index("w_ps"), "w_ps", pa.array((w / gg).astype(np.float32))).filter(pa.array(keep))
        if not len(t):
            continue
        ptmin, ptmax = min(ptmin, pc.min(t["pt"]).as_py()), max(ptmax, pc.max(t["pt"]).as_py())
        nnan += int(pc.sum(pc.cast(pc.is_nan(t["score"]), "int64")).as_py() or 0)
        tiers += np.bincount(t["ps_tier"].to_numpy(), minlength=4)[:4]
        if W is None:
            codec = next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")
            W = pq.ParquetWriter(out, t.schema, compression=codec)
        elif t.schema != W.schema:
            t = t.cast(W.schema)
        W.write_table(t)
        n += len(t)
    if W is None:
        sys.exit(f"STOP {s}: 0 wierszy D+")
    W.close()
    np.savez_compressed(os.path.join(a.out, f"hist_{s}.npz"), **hist)
    b = os.path.getsize(out) + os.path.getsize(os.path.join(a.out, f"hist_{s}.npz"))
    tot += b
    S[s] = {"n_files": len(files), "file_idx": fidx, "n_parts": len(P[s]), "dp_rows": n, "dp_rows_by_tier": tiers.tolist(),
            "dp_pt_min": ptmin, "dp_pt_max": ptmax, "pt_in_range": bool(ptmin >= lo and ptmax < hi), "dp_score_nan": nnan,
            "dp_hist_all_after_veto": int(hist["dp_h"][0].sum()), "d0_hist_all": int(hist["d0_h"][0].sum()),
            "ds_hist_all_rs": int(hist["ds_h"][0, 0].sum()), "ds_hist_all_ws": int(hist["ds_h"][1, 0].sum()),
            "counts": cnt, "bytes": b, "missing_columns": sorted(miss), "columns": pq.read_schema(out).names}
    print(s, {k: v for k, v in S[s].items() if k not in ("file_idx", "columns", "counts")}, flush=True)
thr = None
if a.thr and os.path.exists(a.thr):
    thr = (json.load(open(a.thr)).get("H_s12345") or {}).get("3-5")
summ = {"samples": S, "total_bytes": tot, "max_gb": a.max_gb, "parts_parquet_bytes": byts, "thinning": {"g_tier01": g01, "g_tier2": g2},
        "pt_range": P[samples[0]][0]["pt_range"], "settings": P[samples[0]][0]["settings"], "model": P[samples[0]][0]["model"],
        "ref_thresholds_H_s12345_pt3-5": thr,
        "note": "score H poza zakresem treningu (pT >= 3 GeV); progi WP binu 3-5 tylko jako odniesienie. Wiersze dp_* sa preskalowane "
                "(w_ps, ps_tier, ps_u) — pelna statystyka w hist_*.npz. D0/D*: tylko liczniki i histogramy; d23 zawiera D* tylko z pT(D*) w oknie.",
        "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
json.dump(summ, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
bad = [s for s in samples if not S[s]["pt_in_range"]]
if bad:
    sys.exit(f"STOP: pT poza oknem: {[(s, S[s]['dp_pt_min'], S[s]['dp_pt_max']) for s in bad]}")
if any(S[s]["dp_score_nan"] for s in samples):
    sys.exit("STOP: score NaN")
if tot > a.hard_gb * 1e9:
    sys.exit(f"STOP: wyjście {tot/1e9:.2f} GB > {a.hard_gb} GB (zmniejsz CD_PS_WIN / CD_PS_OUT albo podnieś CD_SC_FULL)")
print(f"OK razem {tot/1e9:.3f} GB", flush=True)
