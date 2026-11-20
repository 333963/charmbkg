#!/usr/bin/env python3
"""AFP_EXT (17.11, projekt B): slim tabela zdarzeń (evtable) + slim kandydaci D⁺/D⁰/D* z h5 v2.0 → parquet (wejście AFP-MIX,
charmbkg.afp_mix; analiza: run_data.py obok). Tylko numpy/h5py/pyarrow.

Zdarzenia (1 wiersz / wpis drzewa, także bez kandydata; nazwy jak w h5 v2.0 = afp_mix.COLMAP):
  event_number u8, lumiBlock u4, actual_mu f4, passGRL_MinBiasAFP i1, nPileupVertices i2, nSelBL i2, xi_cal_A/C f4,
  afp_nProtons_A/C u1 (obcięte do 255), afp_nInWindow_A/C u1, afp_xi1/xi2_A/C f4 (−1 brak), afp_tag_A/C u1,
  trigPassed_<11 łańcuchów> i1, trigTBP_<5 L1> i1, n_cand u2 (wszyscy kandydaci zdarzenia w h5), n_cand_sel u2 (po slim).
  run_number i trigPS (stałe w LB) → meta json (trigps: wiersz na (lumiBlock, zestaw preskal)).
Kandydaci (bez score — brak w h5; luźne okna masy), kolumny: event_number u8, Dmeson_species u1 (0 D⁺, 1 D*, 3 D⁰),
  Dmeson_mass f4 [GeV], Dmeson_delta_mass f4 (D*; inaczej 0), Dmeson_pt f4, Dmeson_eta f4, Dmeson_Lxy f4 [mm].
  Selekcja: D⁺ (decayType 4) m ∈ [1,72; 2,02]; D⁰ (10) m ∈ [1,75; 2,00] i L_xy ≥ 0,25 mm; D* (1) Δm ∈ [0,139; 0,160] (L_xy jak w v2.0, ≥ 0,10).
  Kwantyzacja (lepsza kompresja): m, Δm do 1e-4 GeV; pT 0,01 GeV; η 0,001; L_xy 0,001 mm.
Użycie: afp_extract.py --sample S --h5dir DIR --part K --nparts N --out-prefix WD/S.pK [--max-files M] [--wd WD]
"""
import argparse, glob, json, os, subprocess, sys, time
import numpy as np
import h5py
import pyarrow as pa
import pyarrow.parquet as pq

TRIG = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
        "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_2"]
L1 = ["L1_MBTS_1", "L1_MBTS_2", "L1_J12", "L1_AFP_A_OR_C", "L1_RD0_FILLED"]
EV = {"event_number": "u8", "lumiBlock": "u4", "actual_mu": "f4", "passGRL_MinBiasAFP": "i1", "nPileupVertices": "i2", "nSelBL": "i2",
      "xi_cal_A": "f4", "xi_cal_C": "f4", "afp_nProtons_A": "u1", "afp_nProtons_C": "u1", "afp_nInWindow_A": "u1", "afp_nInWindow_C": "u1",
      "afp_xi1_A": "f4", "afp_xi1_C": "f4", "afp_xi2_A": "f4", "afp_xi2_C": "f4", "afp_tag_A": "u1", "afp_tag_C": "u1"}
CLIP = {"u1": (0, 255), "i2": (-32768, 32767)}
CAND_IN = ["event_number", "Dmeson_decayType", "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_pt", "Dmeson_eta", "Dmeson_Lxy"]
SPEC = {4: 0, 1: 1, 10: 3}
SEL = {"Dplus_m": (1.72, 2.02), "D0_m": (1.75, 2.00), "D0_Lxy_min": 0.25, "Dstar_dm": (0.139, 0.160)}
QUANT = {"Dmeson_mass": 1e-4, "Dmeson_delta_mass": 1e-4, "Dmeson_pt": 0.01, "Dmeson_eta": 0.001, "Dmeson_Lxy": 0.001}
XRD = os.environ.get("XRD", "root://eosuser.cern.ch")


def codec():
    return next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")


def local_copy(path, wd):
    """xrdcp do scratch (szybciej i pewniej niż FUSE); gdy się nie da → ścieżka oryginalna."""
    if wd and path.startswith("/eos/"):
        dst = os.path.join(wd, os.path.basename(path))
        try:
            if subprocess.run(["xrdcp", "-s", "-f", f"{XRD}/{path}", dst], timeout=900).returncode == 0:
                return dst, True
        except Exception:
            pass
        if os.path.exists(dst):
            os.remove(dst)
    return path, False


def read_fields(d, names):
    """Odczyt podzbioru pól compound porcjami (bez kolumn object)."""
    dt = np.dtype([(n, d.dtype[n]) for n in names])
    out = np.empty(d.shape[0], dt)
    step = 500_000
    for i in range(0, d.shape[0], step):
        out[i:i + step] = d.fields(names)[i:i + step]
    return out


def one_file(fn):
    with h5py.File(fn, "r") as f:
        E = f["evtable"]
        miss = [n for n in list(EV) + ["run_number"] + [f"trigPassed_{t}" for t in TRIG] + [f"trigTBP_{t}" for t in L1]
                + [f"trigPS_{t}" for t in TRIG] if n not in E.dtype.names]
        if miss:
            sys.exit(f"{fn}: brak kolumn evtable {miss}")
        e = E[:]
        c = read_fields(f["events"], CAND_IN)
        n_read = int(f.attrs.get("n_evtable", len(e)))
    if n_read != len(e):
        sys.exit(f"{fn}: n_evtable {n_read} ≠ wierszy {len(e)}")
    runs = np.unique(e["run_number"])
    if len(runs) != 1:
        sys.exit(f"{fn}: run_number niejednoznaczny {runs[:5]}")
    if len(np.unique(e["event_number"])) != len(e):
        sys.exit(f"{fn}: event_number nieunikalny w evtable")
    ev = {}
    for n, t in EV.items():
        v = e[n]
        if t in CLIP:
            v = np.clip(v, *CLIP[t])
        ev[n] = v.astype(t)
    for t in TRIG:
        ev[f"trigPassed_{t}"] = e[f"trigPassed_{t}"].astype("i1")
    for t in L1:
        ev[f"trigTBP_{t}"] = e[f"trigTBP_{t}"].astype("i1")
    # preskale: stałe w LB → meta
    ps = np.stack([e["lumiBlock"].astype(np.float64)] + [e[f"trigPS_{t}"].astype(np.float64) for t in TRIG], 1)
    psu, psn = np.unique(ps, axis=0, return_counts=True)
    # kandydaci
    dt, m, dm, lxy = c["Dmeson_decayType"], c["Dmeson_mass"], c["Dmeson_delta_mass"], c["Dmeson_Lxy"]
    sel = ((dt == 4) & (m >= SEL["Dplus_m"][0]) & (m <= SEL["Dplus_m"][1])) \
        | ((dt == 10) & (m >= SEL["D0_m"][0]) & (m <= SEL["D0_m"][1]) & (lxy >= SEL["D0_Lxy_min"])) \
        | ((dt == 1) & (dm >= SEL["Dstar_dm"][0]) & (dm <= SEL["Dstar_dm"][1]))
    nby = {str(k): [int((dt == k).sum()), int((sel & (dt == k)).sum())] for k in SPEC}
    idx = np.searchsorted(np.sort(e["event_number"]), c["event_number"])
    es = np.sort(e["event_number"])
    if len(c) and not np.array_equal(es[np.clip(idx, 0, len(es) - 1)], c["event_number"]):
        sys.exit(f"{fn}: kandydat spoza evtable")
    order = np.argsort(e["event_number"], kind="stable")
    ncand = np.zeros(len(e), np.int64); ncs = np.zeros(len(e), np.int64)
    np.add.at(ncand, order[idx], 1); np.add.at(ncs, order[idx[sel]], 1)
    ev["n_cand"] = np.clip(ncand, 0, 65535).astype("u2"); ev["n_cand_sel"] = np.clip(ncs, 0, 65535).astype("u2")
    cs = c[sel]
    spec = np.zeros(len(cs), "u1")
    for k, v in SPEC.items():
        spec[cs["Dmeson_decayType"] == k] = v
    ca = {"event_number": cs["event_number"].astype("u8"), "Dmeson_species": spec}
    for n, q in QUANT.items():
        v = cs[n].astype(np.float64)
        if n == "Dmeson_delta_mass":
            v = np.where(cs["Dmeson_decayType"] == 1, v, 0.0)
        ca[n] = (np.round(v / q) * q).astype("f4")
    return int(runs[0]), ev, ca, psu, psn, nby


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--h5dir", required=True)
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--out-prefix", required=True); ap.add_argument("--max-files", type=int, default=0)
    ap.add_argument("--wd", default="")
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.h5dir, f"{a.sample}_*_none.h5")))
    if not files:
        sys.exit(f"brak plików {a.h5dir}/{a.sample}_*_none.h5")
    mine = files[a.part::a.nparts]
    if a.max_files:
        mine = mine[:a.max_files]
    cd = codec()
    we = wc = None
    meta = {"sample": a.sample, "part": a.part, "nparts": a.nparts, "n_files_total": len(files), "files": [], "n_events": 0, "n_cands": 0,
            "rows_per_file": {}, "n_cand_by_decayType": {str(k): [0, 0] for k in SPEC}, "trig": TRIG, "l1": L1, "sel": SEL, "quant": QUANT,
            "codec": cd, "run_number": None}
    PS = {}
    t0 = time.time()
    for fn in mine:
        loc, tmp = local_copy(fn, a.wd)
        try:
            run, ev, ca, psu, psn, nby = one_file(loc)
        finally:
            if tmp:
                os.remove(loc)
        if meta["run_number"] not in (None, run):
            sys.exit(f"{fn}: run {run} ≠ {meta['run_number']}")
        meta["run_number"] = run
        te, tc = pa.table(ev), pa.table(ca)
        if we is None:
            we = pq.ParquetWriter(a.out_prefix + ".events.parquet", te.schema, compression=cd)
            wc = pq.ParquetWriter(a.out_prefix + ".cands.parquet", tc.schema, compression=cd)
        we.write_table(te); wc.write_table(tc)
        b = os.path.basename(fn)
        meta["files"].append(b); meta["rows_per_file"][b] = [te.num_rows, tc.num_rows]
        meta["n_events"] += te.num_rows; meta["n_cands"] += tc.num_rows
        for k, (n0, n1) in nby.items():
            meta["n_cand_by_decayType"][k][0] += n0; meta["n_cand_by_decayType"][k][1] += n1
        for r, n in zip(psu, psn):
            PS[tuple(r.tolist())] = PS.get(tuple(r.tolist()), 0) + int(n)
        print(f"{b}: zdarzeń {te.num_rows} kandydatów {tc.num_rows} ({time.time() - t0:.0f} s)", flush=True)
    we.close(); wc.close()
    meta["trigps"] = [list(k) + [n] for k, n in sorted(PS.items())]     # [lumiBlock, PS × len(TRIG), n_zdarzeń]
    meta["bytes"] = [os.path.getsize(a.out_prefix + ".events.parquet"), os.path.getsize(a.out_prefix + ".cands.parquet")]
    meta["seconds"] = round(time.time() - t0, 1)
    json.dump(meta, open(a.out_prefix + ".meta.json", "w"), indent=1)
    print(f"OK {a.sample} część {a.part}/{a.nparts}: {len(mine)} plików, {meta['n_events']} zdarzeń, {meta['n_cands']} kandydatów, "
          f"{meta['bytes'][0] / 1e6:.1f} + {meta['bytes'][1] / 1e6:.1f} MB, {meta['seconds']} s")


if __name__ == "__main__":
    main()
