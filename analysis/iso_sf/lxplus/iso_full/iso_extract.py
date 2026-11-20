#!/usr/bin/env python3
"""ISO-FULL (16.11): slim tabela D0 (wiersz Kπ, decayType 10) z tagiem D* + izolacja i krotność z h5 v1.2.

Ta sama logika co charm_work/ISO_sf/scripts/extract_data.py + evtools.py (bez charmtag; numpy/pandas/h5py/pyarrow):
  - wiersze: parquet T4 hyp_d0_dstar/<próbka>.parquet (tylko D0 z tagiem D*, score D0_H z join), kolumna `row` = wiersz
    GLOBALNY (offset pliku z <próbka>.meta.json "file_offsets" + wiersz lokalny),
  - z h5 v1.2 (wiersz D0): Dmeson_ptcone40/60/80/100, actual_mu, Dmeson_truth_category; kontrola m (|Δ| < 1e-4) i decayType == 10,
  - krotność (brak n_trk w h5): n_cand = liczba wierszy w bloku kolejnych wierszy o tym samym event_number (cały plik h5, wszystkie
    decayType), n_trkid = liczba różnych trackId córek (0–2) wszystkich kandydatów w tym bloku (proxy n_trk).
Wyjście: parquet (jedna część = podzbiór plików: indeks pliku % nparts == part) + <out>.meta.json.

Użycie:
  iso_extract.py --sample data_428770 --tag DIR/data_428770.parquet --h5dir H5V12/data_428770 --out X.parquet [--part k --nparts n]
                 [--only-local]   # test: tylko pliki z meta obecne w --h5dir (lokalnie 2/39); bez tego brak pliku = błąd
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
import h5py
import pyarrow as pa
import pyarrow.parquet as pq

TAG_COLS = ["sample", "file", "row", "evblk", "run", "event", "m", "pt", "eta", "phi", "dm", "score", "score_src", "fold", "Lxy",
            "pis_pt", "soft_pi_charge", "dstar_tag", "dstar_conflict", "dstar_n", "dstar_row", "truth", "hypo_swap", "hypo_from_tag",
            "d_trk0", "d_trk1", "d_pt0", "d_pt1"]
H5_ROW = {"Dmeson_ptcone40": "ptcone40", "Dmeson_ptcone60": "ptcone60", "Dmeson_ptcone80": "ptcone80",
          "Dmeson_ptcone100": "ptcone100", "actual_mu": "mu", "Dmeson_truth_category": "truth_cat"}
H5_CHK = ["Dmeson_mass", "Dmeson_decayType"]
H5_EVT = ["event_number"] + [f"Dmeson_daughter{k}_trackId" for k in range(3)]


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def event_mult(ev, tid):
    """= ISO_sf/scripts/evtools.event_mult, ale na gotowych tablicach (ev: event_number, tid: (n,3) trackId)."""
    n = len(ev)
    brk = np.r_[0, np.flatnonzero(ev[1:] != ev[:-1]) + 1, n]
    blk = np.repeat(np.arange(len(brk) - 1), np.diff(brk))
    ncand = np.diff(brk)[blk]
    b = np.repeat(blk, 3); t = tid.astype(np.int64).ravel(); ok = t >= 0
    key = b[ok].astype(np.int64) * (1 << 24) + t[ok]
    u = np.unique(key)
    ntr = np.bincount((u >> 24).astype(np.int64), minlength=len(brk) - 1)
    return ncand.astype(np.int32), ntr[blk].astype(np.int32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tag", required=True, help="parquet T4 (hyp_d0_dstar/<próbka>.parquet); meta obok: <próbka>.meta.json")
    ap.add_argument("--h5dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--only-local", action="store_true")
    a = ap.parse_args()
    meta = json.load(open(a.tag[:-len(".parquet")] + ".meta.json"))
    files, offs = meta["files"], [int(x) for x in meta["file_offsets"]]
    fidx = [i for i in range(len(files)) if i % a.nparts == a.part]
    if a.only_local:
        have = set(os.listdir(a.h5dir)); fidx = [i for i in fidx if files[i] in have]
    else:
        miss = [files[i] for i in fidx if not os.path.exists(os.path.join(a.h5dir, files[i]))]
        if miss:
            sys.exit(f"brak plików h5 v1.2: {miss}")
    if not fidx:
        sys.exit("0 plików do przetworzenia")
    pf = pq.ParquetFile(a.tag)
    cols = [c for c in TAG_COLS if c in pf.schema_arrow.names]
    log(f"{a.sample}: część {a.part}/{a.nparts}, {len(fidx)} plików, kolumny tagu {len(cols)}")
    codec = next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")
    tmp = a.out + ".tmp"
    W, nout, per_file, nmis = None, 0, {}, 0
    for fi in fidx:
        t0 = time.time()
        s = pq.read_table(a.tag, columns=cols, filters=[("file", "=", fi)]).to_pandas()
        if s.empty:
            per_file[files[fi]] = 0; continue
        lr = s.row.to_numpy(np.int64) - offs[fi]
        with h5py.File(os.path.join(a.h5dir, files[fi]), "r") as f:
            ds = f["events"]
            nrow = ds.shape[0]
            if fi + 1 < len(offs) and offs[fi + 1] - offs[fi] != nrow:
                sys.exit(f"{files[fi]}: {nrow} wierszy ≠ offsety meta {offs[fi + 1] - offs[fi]}")
            if lr.min() < 0 or lr.max() >= nrow:
                sys.exit(f"{files[fi]}: row poza plikiem (min {lr.min()}, max {lr.max()}, n {nrow})")
            A = ds.fields(H5_EVT + H5_CHK + list(H5_ROW))[:]          # jeden przebieg po pliku
        ncand, ntrk = event_mult(A["event_number"], np.stack([A[c] for c in H5_EVT[1:]], 1))
        mm = np.abs(A["Dmeson_mass"][lr].astype(np.float64) - s.m.to_numpy(np.float64)) > 1e-4
        if mm.any() or (A["Dmeson_decayType"][lr] != 10).any():
            sys.exit(f"{files[fi]}: niezgodne m ({int(mm.sum())}) albo decayType ≠ 10 — to nie te same wiersze")
        for k, v in H5_ROW.items():
            s[v] = A[k][lr]
        s["n_cand"], s["n_trkid"] = ncand[lr], ntrk[lr]
        for c in s.columns:
            if s[c].dtype == np.float64 and c != "m":
                s[c] = s[c].astype(np.float32)
        tb = pa.Table.from_pandas(s, preserve_index=False)
        if W is None:
            W = pq.ParquetWriter(tmp, tb.schema, compression=codec)
        elif tb.schema != W.schema:
            tb = tb.cast(W.schema)
        W.write_table(tb)
        nout += len(s); per_file[files[fi]] = len(s)
        log(f"  {files[fi]}: {len(s)} wierszy, n_h5 {nrow}, {time.time() - t0:.0f} s (razem {nout})")
        del A
    if W is None:
        sys.exit("0 wierszy")
    W.close(); os.replace(tmp, a.out)
    m = {"sample": a.sample, "tag": a.tag, "h5dir": a.h5dir, "part": a.part, "nparts": a.nparts, "only_local": a.only_local,
         "files": [files[i] for i in fidx], "file_index": fidx, "rows_per_file": per_file, "n_rows_out": nout,
         "codec": codec, "bytes": os.path.getsize(a.out), "time": time.strftime("%FT%T"),
         "note": "row = wiersz globalny (offset meta T4); ptcone z wiersza D0; n_trkid = proxy krotności (ISO_sf/evtools)"}
    json.dump(m, open(a.out[:-len(".parquet")] + ".meta.json", "w"), indent=1)
    log(f"OK {a.out}: {nout} wierszy, {m['bytes'] / 1e6:.1f} MB ({codec})")


if __name__ == "__main__":
    main()
