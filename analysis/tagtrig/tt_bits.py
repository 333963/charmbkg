#!/usr/bin/env python3
"""TAGTRIG / węzeł bits (18.11.2026, projekt B): bity WSZYSTKICH łańcuchów trygera z ROOT v2.0 (drzewo CharmAnalysis) → parquet.
Tylko numpy/awkward/uproot/pyarrow. Tylko zdarzenia z passGRL_MinBiasAFP = 1.

Wyjścia (części, scalane przez tt_merge.py):
  <o>.events.parquet — wiersz na zdarzenie z tagiem AFP (T0 = GRL ∧ (tag_A ∨ tag_C); tag = dokładnie 1 proton po stronie i 0,030 ≤ ξ ≤ 0,100):
      run_number u4, event_number u8, lumiBlock u4, npu i2 (nPileupVertices), actual_mu f4, average_mu f4, nSelBL i2,
      tag_A/C u1, afp_nProtons_A/C u1, xi_A/C f4 (ξ protonu tagu; −1 gdy strona bez tagu), xi_cal_A/C f4,
      n_D0, n_Dp, n_Dst u2 — LICZNIKI kandydatów po selekcji slim AFP_mix (D⁺ m ∈ [1,72; 2,02]; D⁰ m ∈ [1,75; 2,00] i L_xy ≥ 0,25 mm;
      D* Δm ∈ [0,139; 0,160]; wszędzie L_xy ≥ 0) — żadnych mas, żadnych histogramów,
      trig_w0.. u8 (bit j łańcucha chains[j] → słowo j // 64, bit j % 64; trigPassed po preskali), tbp_w0.. u8 (L1 „before prescale”, l1tbp[j]).
  <o>.ps.parquet — zdarzenia GRL, per (kind 0 = HLT trigPassed/trigPS, 1 = L1 trigTBP; idx; lumiBlock; wartość preskali): n_ev, n_pass
      (ps = −1: brak gałęzi trigPS / pozycja L1).
  <o>.lb.parquet — per lumiBlock: n_all (wpisy drzewa), n_grl, sum_mu_grl.
  <o>.meta.json — pliki, liczby, lista łańcuchów (chains) i pozycji L1 (l1tbp).
Użycie: tt_bits.py --sample data_428770 --rootdir DIR --part K --nparts N --out-prefix WD/S.pK [--wd WD] [--no-copy] [--max-files M] [--max-events N]
"""
import argparse, glob, json, os, subprocess, sys, time
import numpy as np
import awkward as ak
import uproot
import pyarrow as pa
import pyarrow.parquet as pq

GeV = 0.001
XI_MIN, XI_MAX = 0.030, 0.100                      # okno tagu (charmpp Structs.h tag_reco; jak charmpp_h5_v20.py)
SEL = {"Dplus_m": (1.72, 2.02), "D0_m": (1.75, 2.00), "D0_Lxy_min": 0.25, "Dstar_dm": (0.139, 0.160)}   # = afp_ext/afp_extract.py SEL
PRE = {"P": "EventInfo_trigPassed_", "S": "EventInfo_trigPS_", "T": "EventInfo_trigTBP_"}
EVB = ["EventInfo_runNumber", "EventInfo_eventNumber", "EventInfo_lumiBlock", "EventInfo_actualInteractionsPerCrossing",
       "EventInfo_averageInteractionsPerCrossing", "EventInfo_passGRL_MinBiasAFP", "EventInfo_nPileupVertices", "EventInfo_nSelBL",
       "EventInfo_xi_cal_A", "EventInfo_xi_cal_C", "AfpProtons_xi", "AfpProtons_side",
       "DMesons_decayType", "DMesons_m", "DMesons_DeltaMass", "DMesons_fitOutput__Lxy"]
XRD = os.environ.get("XRD", "root://eosuser.cern.ch")


def strip(n):
    return n[:-6] if n.endswith("_NOSYS") else n


def resolve(keys, name):
    if name in keys:
        return name
    if name + "_NOSYS" in keys:
        return name + "_NOSYS"
    raise KeyError(name)


def discover(keys):
    """Gałęzie trygera w drzewie: {nazwa łańcucha: gałąź} dla trigPassed (P), trigPS (S), trigTBP (T)."""
    return {k: {strip(b[len(p):]): b for b in keys if b.startswith(p)} for k, p in PRE.items()}


def pack_bits(M):
    """bool (n, k) → uint64 (n, max(1, ceil(k / 64))); bit j → słowo j // 64, pozycja j % 64."""
    M = np.asarray(M, bool)
    n, k = M.shape
    W = np.zeros((n, max(1, (k + 63) // 64)), np.uint64)
    for j in range(k):
        W[:, j // 64] |= M[:, j].astype(np.uint64) << np.uint64(j % 64)
    return W


def unpack_bits(W, k):
    W = np.asarray(W, np.uint64)
    M = np.zeros((W.shape[0], k), bool)
    for j in range(k):
        M[:, j] = ((W[:, j // 64] >> np.uint64(j % 64)) & np.uint64(1)).astype(bool)
    return M


def afp_sides(xi, side):
    """Per strona (0 = A, 1 = C): n protonów, n w oknie, najmniejsze ξ (−1 brak), tag = (n == 1 ∧ n_w_oknie == 1). Jak charmpp_h5_v20.afp_summary."""
    n = np.minimum(ak.to_numpy(ak.num(xi)), ak.to_numpy(ak.num(side)))
    xi = xi[ak.local_index(xi) < n]
    side = side[ak.local_index(side) < n]
    out = {}
    for s, S in ((0, "A"), (1, "C")):
        v = ak.sort(xi[side == s], axis=1)
        ns = ak.to_numpy(ak.num(v)).astype(np.int64)
        nw = ak.to_numpy(ak.sum((v >= XI_MIN) & (v <= XI_MAX), axis=1)).astype(np.int64)
        x1 = ak.to_numpy(ak.fill_none(ak.firsts(v, axis=1), -1.0)).astype(np.float32)
        out[S] = (ns, nw, x1, (ns == 1) & (nw == 1))
    return out


def cand_counts(dt, m, dm, lxy):
    """Liczniki kandydatów na zdarzenie (D⁰, D⁺, D*) po selekcji slim AFP_mix; m, Δm w MeV (ROOT)."""
    m = m * GeV
    dm = dm * GeV
    base = lxy >= 0.0
    d0 = base & (dt == 10) & (m >= SEL["D0_m"][0]) & (m <= SEL["D0_m"][1]) & (lxy >= SEL["D0_Lxy_min"])
    dp = base & (dt == 4) & (m >= SEL["Dplus_m"][0]) & (m <= SEL["Dplus_m"][1])
    ds = base & (dt == 1) & (dm >= SEL["Dstar_dm"][0]) & (dm <= SEL["Dstar_dm"][1])
    return [np.clip(ak.to_numpy(ak.sum(x, axis=1)), 0, 65535).astype(np.uint16) for x in (d0, dp, ds)]


def acc_ps(acc, kind, lb, passed, ps):
    """acc[(kind, idx, lumiBlock, ps)] += [n_ev, n_pass]; ps None → −1."""
    if len(lb) == 0:
        return
    ulb, ilb = np.unique(lb, return_inverse=True)
    for j in range(passed.shape[1]):
        if ps is None:
            up, ip = np.array([-1.0]), np.zeros(len(lb), np.int64)
        else:
            up, ip = np.unique(ps[:, j], return_inverse=True)
        code = ilb * len(up) + ip
        n = np.bincount(code, minlength=len(ulb) * len(up))
        k = np.bincount(code, weights=passed[:, j], minlength=len(ulb) * len(up))
        for c in np.nonzero(n)[0]:
            key = (kind, j, int(ulb[c // len(up)]), float(up[c % len(up)]))
            v = acc.setdefault(key, [0, 0])
            v[0] += int(n[c]); v[1] += int(round(k[c]))


def process_chunk(a, chains, l1, have_ps, acc, lbacc):
    """a: {nazwa logiczna: tablica}; zwraca kolumny wierszy zdarzeń tagowanych; uzupełnia acc (preskale per LB) i lbacc."""
    npf = lambda k, t: ak.to_numpy(a[k]).astype(t)
    lb = npf("EventInfo_lumiBlock", np.int64)
    grl = npf("EventInfo_passGRL_MinBiasAFP", np.int64) > 0
    mu = npf("EventInfo_actualInteractionsPerCrossing", np.float64)
    n = len(lb)
    P = np.zeros((n, len(chains)), bool)
    S = np.full((n, len(chains)), -1.0)
    for j, c in enumerate(chains):
        P[:, j] = ak.to_numpy(a["P:" + c]) > 0
        if c in have_ps:
            S[:, j] = np.nan_to_num(ak.to_numpy(a["S:" + c]).astype(np.float64), nan=-1.0)
    T = np.zeros((n, len(l1)), bool)
    for j, c in enumerate(l1):
        T[:, j] = ak.to_numpy(a["T:" + c]) > 0
    acc_ps(acc, 0, lb[grl], P[grl], S[grl])
    acc_ps(acc, 1, lb[grl], T[grl], None)
    for u in np.unique(lb):
        m = lb == u
        v = lbacc.setdefault(int(u), [0, 0, 0.0])
        v[0] += int(m.sum()); v[1] += int((m & grl).sum()); v[2] += float(mu[m & grl].sum())
    afp = afp_sides(a["AfpProtons_xi"], a["AfpProtons_side"])
    tagA, tagC = afp["A"][3], afp["C"][3]
    sel = grl & (tagA | tagC)
    nD0, nDp, nDs = cand_counts(a["DMesons_decayType"], a["DMesons_m"], a["DMesons_DeltaMass"], a["DMesons_fitOutput__Lxy"])
    ev = {"run_number": npf("EventInfo_runNumber", np.uint32)[sel], "event_number": npf("EventInfo_eventNumber", np.uint64)[sel],
          "lumiBlock": lb[sel].astype(np.uint32), "npu": np.clip(npf("EventInfo_nPileupVertices", np.int64), -32768, 32767)[sel].astype(np.int16),
          "actual_mu": mu[sel].astype(np.float32), "average_mu": npf("EventInfo_averageInteractionsPerCrossing", np.float32)[sel],
          "nSelBL": np.clip(npf("EventInfo_nSelBL", np.int64), -32768, 32767)[sel].astype(np.int16)}
    for S_ in ("A", "C"):
        ns, _, x1, tag = afp[S_]
        ev["tag_" + S_] = tag[sel].astype(np.uint8)
        ev["afp_nProtons_" + S_] = np.clip(ns, 0, 255)[sel].astype(np.uint8)
        ev["xi_" + S_] = np.where(tag, x1, np.float32(-1.0))[sel].astype(np.float32)
        ev["xi_cal_" + S_] = npf("EventInfo_xi_cal_" + S_, np.float32)[sel]
    ev["n_D0"], ev["n_Dp"], ev["n_Dst"] = nD0[sel], nDp[sel], nDs[sel]
    W = pack_bits(P[sel])
    for i in range(W.shape[1]):
        ev[f"trig_w{i}"] = W[:, i]
    W = pack_bits(T[sel])
    for i in range(W.shape[1]):
        ev[f"tbp_w{i}"] = W[:, i]
    return ev, n, int(grl.sum())


def codec():
    return next((c for c in ("zstd", "gzip", "snappy") if pa.Codec.is_available(c)), "none")


def local_copy(path, wd):
    """xrdcp do scratch (pewniej niż FUSE); gdy się nie da → ścieżka oryginalna."""
    if wd and path.startswith("/eos/"):
        dst = os.path.join(wd, os.path.basename(path))
        try:
            if subprocess.run(["xrdcp", "-s", "-f", f"{XRD}/{path}", dst], timeout=1800).returncode == 0:
                return dst, True
        except Exception:
            pass
        if os.path.exists(dst):
            os.remove(dst)
    return path, False


def list_files(rootdir, ptprod):
    return sorted(glob.glob(os.path.join(rootdir, "**", f"*{ptprod}*.root"), recursive=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--rootdir", required=True)
    ap.add_argument("--ptprod", default="dpt3to150")
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--out-prefix", required=True); ap.add_argument("--wd", default="")
    ap.add_argument("--no-copy", action="store_true"); ap.add_argument("--step", default="100 MB")
    ap.add_argument("--max-files", type=int, default=0); ap.add_argument("--max-events", type=int, default=0)
    a = ap.parse_args()
    files = list_files(a.rootdir, a.ptprod)
    if not files:
        sys.exit(f"brak plików {a.rootdir}/**/*{a.ptprod}*.root")
    mine = files[a.part::a.nparts]
    if a.max_files:
        mine = mine[:a.max_files]
    cd = codec()
    step = a.step if not str(a.step).isdigit() else int(a.step)
    if isinstance(step, str) and step.upper().endswith("MB") and " " not in step:
        step = step[:-2] + " MB"
    meta = {"sample": a.sample, "part": a.part, "nparts": a.nparts, "n_files_total": len(files), "files": [], "rows_per_file": {},
            "n_events": 0, "n_grl": 0, "n_tag": 0, "run_number": None, "chains": None, "l1tbp": None, "ps_missing": [], "codec": cd,
            "sel": SEL, "xi_window": [XI_MIN, XI_MAX], "max_events": a.max_events}
    acc, lbacc, wr, t0 = {}, {}, None, time.time()
    for fn in mine:
        if a.max_events and meta["n_events"] >= a.max_events:
            break
        loc, tmp = (fn, False) if a.no_copy else local_copy(fn, a.wd)
        try:
            with uproot.open(loc) as f:
                t = f["CharmAnalysis"]
                keys = set(t.keys())
                D = discover(keys)
                chains, l1 = sorted(D["P"]), sorted(D["T"])
                if not chains:
                    sys.exit(f"{fn}: brak gałęzi EventInfo_trigPassed_* (to nie dane v2.0?)")
                if meta["chains"] is None:
                    meta["chains"], meta["l1tbp"] = chains, l1
                    meta["ps_missing"] = [c for c in chains if c not in D["S"]]
                elif chains != meta["chains"] or l1 != meta["l1tbp"]:
                    sys.exit(f"{fn}: inna lista łańcuchów niż w pierwszym pliku części ({len(chains)}/{len(l1)} vs {len(meta['chains'])}/{len(meta['l1tbp'])})")
                have_ps = set(D["S"])
                want = {b: resolve(keys, b) for b in EVB}
                want.update({"P:" + c: D["P"][c] for c in chains})
                want.update({"S:" + c: D["S"][c] for c in chains if c in have_ps})
                want.update({"T:" + c: D["T"][c] for c in l1})
                stop = None if not a.max_events else max(0, min(t.num_entries, a.max_events - meta["n_events"]))
                nf = [0, 0]
                if t.num_entries and stop != 0:
                    for arr in t.iterate(sorted(set(want.values())), step_size=step, entry_stop=stop, library="ak"):
                        ev, n, ng = process_chunk({k: arr[b] for k, b in want.items()}, chains, l1, have_ps, acc, lbacc)
                        runs = np.unique(ak.to_numpy(arr[want["EventInfo_runNumber"]]))
                        if len(runs) != 1 or meta["run_number"] not in (None, int(runs[0])):
                            sys.exit(f"{fn}: run_number niejednoznaczny {runs[:5]} (dotąd {meta['run_number']})")
                        meta["run_number"] = int(runs[0])
                        tb = pa.table(ev)
                        if wr is None:
                            wr = pq.ParquetWriter(a.out_prefix + ".events.parquet", tb.schema, compression=cd)
                        wr.write_table(tb)
                        nf[0] += n; nf[1] += tb.num_rows
                        meta["n_events"] += n; meta["n_grl"] += ng; meta["n_tag"] += tb.num_rows
        finally:
            if tmp:
                os.remove(loc)
        b = os.path.basename(fn)
        meta["files"].append(b); meta["rows_per_file"][b] = nf
        print(f"{b}: zdarzeń {nf[0]}, tagowanych (GRL) {nf[1]} ({time.time() - t0:.0f} s)", flush=True)
    if wr is None:
        sys.exit(f"{a.sample} część {a.part}: żadnego zdarzenia w {len(mine)} plikach")
    wr.close()
    K = sorted(acc)
    pq.write_table(pa.table({"kind": np.array([k[0] for k in K], np.uint8), "idx": np.array([k[1] for k in K], np.uint16),
                             "lumiBlock": np.array([k[2] for k in K], np.uint32), "ps": np.array([k[3] for k in K], np.float64),
                             "n_ev": np.array([acc[k][0] for k in K], np.int64), "n_pass": np.array([acc[k][1] for k in K], np.int64)}),
                   a.out_prefix + ".ps.parquet", compression=cd)
    L = sorted(lbacc)
    pq.write_table(pa.table({"lumiBlock": np.array(L, np.uint32), "n_all": np.array([lbacc[k][0] for k in L], np.int64),
                             "n_grl": np.array([lbacc[k][1] for k in L], np.int64), "sum_mu_grl": np.array([lbacc[k][2] for k in L], np.float64)}),
                   a.out_prefix + ".lb.parquet", compression=cd)
    meta["bytes"] = [os.path.getsize(a.out_prefix + x) for x in (".events.parquet", ".ps.parquet", ".lb.parquet")]
    meta["seconds"] = round(time.time() - t0, 1)
    json.dump(meta, open(a.out_prefix + ".meta.json", "w"), indent=1)
    print(f"OK {a.sample} część {a.part}/{a.nparts}: {len(meta['files'])} plików, {meta['n_events']} zdarzeń, GRL {meta['n_grl']}, tag {meta['n_tag']}, "
          f"łańcuchów {len(meta['chains'])} (+{len(meta['l1tbp'])} L1 TBP), {meta['bytes'][0] / 1e6:.1f} MB, {meta['seconds']} s")


if __name__ == "__main__":
    main()
