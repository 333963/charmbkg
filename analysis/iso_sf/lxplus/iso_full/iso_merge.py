#!/usr/bin/env python3
"""ISO-FULL finish: scala <out>/parts/<próbka>_tag.p<k>.parquet → <wd>/<próbka>_tag.parquet (+ .meta.json) i iso_summary.json.
Kontrola kompletności: wszystkie części obecne, pliki z meta T4 pokryte dokładnie raz, suma wierszy == wiersze T4 (n_rows_out).
Użycie: iso_merge.py --out ISO_OUT --tag-dir TAG_DIR --wd WD --samples "data_428770 fd_nd mb_nd" --nparts-data 3
"""
import argparse, json, os, sys
import pyarrow as pa
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--tag-dir", required=True); ap.add_argument("--wd", required=True)
ap.add_argument("--samples", required=True); ap.add_argument("--nparts-data", type=int, required=True)
a = ap.parse_args()
codec = next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")
S, err = {}, []
for s in a.samples.replace(",", " ").split():
    n = a.nparts_data if s.startswith("data_") else 1
    T = json.load(open(os.path.join(a.tag_dir, f"{s}.meta.json")))
    parts = [os.path.join(a.out, "parts", f"{s}_tag.p{k}") for k in range(n)]
    miss = [p for p in parts if not os.path.exists(p + ".parquet")]
    if miss:
        err.append(f"{s}: brak części {miss}"); continue
    metas = [json.load(open(p + ".meta.json")) for p in parts]
    files = sorted(f for m in metas for f in m["files"])
    rows = sum(m["n_rows_out"] for m in metas)
    if files != sorted(T["files"]):
        err.append(f"{s}: pliki części ≠ pliki T4 ({len(files)} vs {len(T['files'])})")
    if rows != T["n_rows_out"]:
        err.append(f"{s}: wierszy {rows} ≠ T4 {T['n_rows_out']}")
    tabs = [pq.read_table(p + ".parquet") for p in parts]
    t = pa.concat_tables([x.cast(tabs[0].schema) for x in tabs])
    o = os.path.join(a.wd, f"{s}_tag.parquet")
    pq.write_table(t, o, compression=codec, row_group_size=500_000)
    rpf = {}
    for m in metas:
        rpf.update(m["rows_per_file"])
    S[s] = {"n_rows_out": t.num_rows, "n_rows_t4": T["n_rows_out"], "n_files": len(files), "nparts": n, "bytes": os.path.getsize(o),
            "codec": codec, "columns": t.schema.names}
    json.dump(dict(S[s], sample=s, files=T["files"], file_offsets=T["file_offsets"], rows_per_file=rpf,
                   note="ISO-FULL: wiersze T4 hyp_d0_dstar + ptcone40/60/80/100, mu, truth_cat (wiersz D0), n_cand, n_trkid (proxy krotności)"),
              open(os.path.join(a.wd, f"{s}_tag.meta.json"), "w"), indent=1)
    print(s, S[s]["n_rows_out"], f"{S[s]['bytes'] / 1e6:.1f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump({"samples": S}, open(os.path.join(a.wd, "iso_summary.json"), "w"), indent=1)
