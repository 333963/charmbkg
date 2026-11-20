#!/usr/bin/env python3
"""AFP_EXT finish: scala <out>/parts/<próbka>.p<k>.{events,cands}.parquet → <wd>/<próbka>_{events,cands}.parquet + <próbka>.meta.json
i afp_summary.json. Kontrola: wszystkie części, każdy plik h5 dokładnie raz, liczba plików = oczekiwana, sumy wierszy = meta części.
Strumieniowo (ParquetWriter, grupy wierszy) — mało RAM. Użycie: afp_merge.py --out AFP_OUT --wd WD --samples "a b" --nparts 8 --nfiles "a:97 b:458"
"""
import argparse, json, os, sys
import pyarrow as pa
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--wd", required=True); ap.add_argument("--samples", required=True)
ap.add_argument("--nparts", type=int, required=True); ap.add_argument("--nfiles", default="")
a = ap.parse_args()
nexp = dict((kv.split(":")[0], int(kv.split(":")[1])) for kv in a.nfiles.replace(",", " ").split())
codec = next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")
S, err = {}, []
for s in a.samples.replace(",", " ").split():
    parts = [os.path.join(a.out, "parts", f"{s}.p{k}") for k in range(a.nparts)]
    miss = [p for p in parts if not (os.path.exists(p + ".meta.json") and os.path.exists(p + ".events.parquet") and os.path.exists(p + ".cands.parquet"))]
    if miss:
        err.append(f"{s}: brak części {[os.path.basename(p) for p in miss]}"); continue
    metas = [json.load(open(p + ".meta.json")) for p in parts]
    files = sorted(f for m in metas for f in m["files"])
    if len(set(files)) != len(files):
        err.append(f"{s}: pliki powtórzone w częściach")
    if any(m["n_files_total"] != len(files) for m in metas) or (s in nexp and nexp[s] != len(files)):
        err.append(f"{s}: plików {len(files)}, części widziały {sorted(set(m['n_files_total'] for m in metas))}, oczekiwane {nexp.get(s)}")
    if len(set(m["run_number"] for m in metas)) != 1:
        err.append(f"{s}: różne run_number w częściach")
    info = {"sample": s, "run_number": metas[0]["run_number"], "n_files": len(files), "nparts": a.nparts, "codec": codec}
    for kind, key in (("events", "n_events"), ("cands", "n_cands")):
        o = os.path.join(a.wd, f"{s}_{kind}.parquet")
        w, n = None, 0
        for p in parts:
            pf = pq.ParquetFile(f"{p}.{kind}.parquet")
            for b in pf.iter_batches(batch_size=2_000_000):
                t = pa.Table.from_batches([b])
                if w is None:
                    w = pq.ParquetWriter(o, t.schema, compression=codec)
                w.write_table(t.cast(w.schema)); n += t.num_rows
        w.close()
        exp = sum(m[key] for m in metas)
        if n != exp:
            err.append(f"{s}: {kind} wierszy {n} ≠ meta części {exp}")
        info[key] = n; info[f"bytes_{kind}"] = os.path.getsize(o)
        info[f"columns_{kind}"] = pq.ParquetFile(o).schema_arrow.names
    PS = {}
    for m in metas:
        for r in m["trigps"]:
            PS[tuple(r[:-1])] = PS.get(tuple(r[:-1]), 0) + r[-1]
    nby = {k: [sum(m["n_cand_by_decayType"][k][i] for m in metas) for i in (0, 1)] for k in metas[0]["n_cand_by_decayType"]}
    rpf = {}
    for m in metas:
        rpf.update(m["rows_per_file"])
    S[s] = info
    json.dump(dict(info, files=files, rows_per_file=rpf, trig=metas[0]["trig"], l1=metas[0]["l1"], sel=metas[0]["sel"], quant=metas[0]["quant"],
                   n_cand_by_decayType=nby, trigps_columns=["lumiBlock"] + metas[0]["trig"] + ["n_events"],
                   trigps=[list(k) + [n] for k, n in sorted(PS.items())],
                   note="AFP_EXT: slim evtable + kandydaci D+/D0/D* z h5 v2.0 (bez score); run_number stały = meta; n_cand_by_decayType = [w h5, po slim]"),
              open(os.path.join(a.wd, f"{s}.meta.json"), "w"), indent=1)
    print(s, info["n_events"], "zdarzeń", info["n_cands"], "kandydatów", f"{info['bytes_events'] / 1e6:.1f} + {info['bytes_cands'] / 1e6:.1f} MB", flush=True)
if err:
    sys.exit("BŁĄD: " + "; ".join(err))
json.dump({"samples": S, "bytes_total": sum(v["bytes_events"] + v["bytes_cands"] for v in S.values())}, open(os.path.join(a.wd, "afp_summary.json"), "w"), indent=1)
