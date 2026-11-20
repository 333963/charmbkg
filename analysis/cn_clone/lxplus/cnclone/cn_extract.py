#!/usr/bin/env python3
"""CN_clone — ekstrakcja z h5 v2.0 (compound `events`) dla klonu CharmNet D0. Bez charmtag (tylko numpy + h5py).
Na wiersz D0 (decayType 10) po preselekcji treningu CharmNet (pT ≥ 3 GeV, |η| < 2,2, min pT córki > 0,8 GeV) i L_xy ≥ --lxy-min:
  a_*  kandydaci inkluzywni, m ∈ [--all-lo, --all-hi]: X (14 wejść), m, pt, eta, truth, hypo_swap, uid, bity triggera;
       podpróbka zdarzeń --all-frac (hash event_uid); fd: tylko sygnał z prawdą (truth 0/1), chyba że --fd-bkg;
  p_*  pary D0←D* (jak B_trig/btrig_extract: |m − m0| < 75 MeV, 0 < Δm < 165 MeV, |m_D0 − (m_D* − Δm)| < 1 MeV,
       0,8·pT(D*) < pT(D0) ≤ 1,001·pT(D*), |Δη| < 0,1), 100 % zdarzeń: X, m, pt, eta, dm, truth, truth_ds, pt_ds, pis_pt, …
Klucz zdarzenia: event_uid (MC: (run, event) nie jest kluczem). Dane: truth = −1 (h5 v2.0 ma tam 0).
Użycie: cn_extract.py --tag data|fd|mb --out OUT.npz [--all-frac 0.05] [--lxy-min 0.10] [--max-rows N] FILE.h5 [...]"""
import argparse
import json
import os
import sys
import time

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cn_feat as CF  # noqa: E402

DM_MAX, M_WIN, TOL_M = 0.165, 0.075, 0.001


def read_rows(path, chunk, max_rows):
    """Wiersze D0 i D* (tylko potrzebne kolumny) + numer wiersza w pliku."""
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        ev = {k: c for k, c in CF.EVENT_COLS.items() if c in names}
        tp = {c: "trigPassed_" + c for c in CF.TRIG_CHAINS if "trigPassed_" + c in names}
        miss = [c for c in CF.RAW_COLS if c not in names]
        if miss:
            raise SystemExit(f"{path}: brak kolumn {miss}")
        cols = CF.RAW_COLS + list(ev.values()) + list(tp.values())
        n = ds.shape[0] if not max_rows else min(ds.shape[0], max_rows)
        out, rows = [], []
        for a in range(0, n, chunk):
            b = min(n, a + chunk)
            x = ds.fields(cols)[a:b]
            k = (x["Dmeson_decayType"] == 10) | (x["Dmeson_decayType"] == 1)
            out.append(x[k])
            rows.append(np.arange(a, b, dtype=np.int64)[k])
    return np.concatenate(out), np.concatenate(rows), ev, tp, n


def table(R, rows, X, ev, tp, tag, fi, prefix):
    o = {"uid": R["event_uid"].astype(np.uint64), "run": R["run_number"].astype(np.int64), "event": R["event_number"].astype(np.int64),
         "row": rows, "file_idx": np.full(len(R), fi, np.int32), "m": R["Dmeson_mass"].astype(np.float32),
         "pt": R["Dmeson_pt"].astype(np.float32), "eta": R["Dmeson_eta"].astype(np.float32),
         "truth": (R["Dmeson_truth_category"].astype(np.int16) if tag != "data" else np.full(len(R), -1, np.int16)),
         "swap": R["Dmeson_hypo_swap"].astype(np.int8), "X": X}
    for k, c in ev.items():
        o["e_" + k] = R[c].astype(np.float32 if k == "lead_pt" else np.int32)
    for c, col in tp.items():
        o["e_tp_" + c] = R[col].astype(np.int8)
    return {prefix + k: v for k, v in o.items()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, choices=["data", "fd", "mb"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=400_000)
    ap.add_argument("--all-lo", type=float, default=1.75)
    ap.add_argument("--all-hi", type=float, default=2.00)
    ap.add_argument("--all-frac", type=float, default=1.0)
    ap.add_argument("--lxy-min", type=float, default=0.10, help="L_xy(D0), L_xy(D*) ≥ [mm]; 0,10 = podłoga plików (jak u CharmNet)")
    ap.add_argument("--fd-bkg", action="store_true", help="fd: zachowaj też nie-sygnał w a_* (domyślnie tylko truth 0/1)")
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    P = CF.PRESEL
    parts, t0, stat = {}, time.time(), []
    evk, tpk = {}, {}
    for fi, path in enumerate(a.files):
        R, rows, ev, tp, nread = read_rows(path, a.chunk, a.max_rows)
        evk.update(ev); tpk.update(tp)
        is0 = R["Dmeson_decayType"] == 10
        R0, rows0 = R[is0], rows[is0]
        X, dmin = CF.features(R0)
        m = R0["Dmeson_mass"].astype(np.float64)
        ok = (R0["Dmeson_pt"] >= P["pt_min"]) & (np.abs(R0["Dmeson_eta"]) < P["abs_eta_max"]) & (dmin > P["dau_min_pt"]) & \
             (R0["Dmeson_Lxy"] >= a.lxy_min)
        R0, rows0, X, m = R0[ok], rows0[ok], X[ok], m[ok]
        # --- inkluzywne
        aw = (m >= a.all_lo) & (m <= a.all_hi)
        if a.all_frac < 1.0:
            aw &= CF.unit_hash(R0["event_uid"], salt=3) < a.all_frac
        if a.tag == "fd" and not a.fd_bkg:
            aw &= np.isin(R0["Dmeson_truth_category"], CF.SIGNAL_CATS)
        T = table(R0[aw], rows0[aw], X[aw], ev, tp, a.tag, fi, "a_")
        # --- pary D0←D*
        Rs, rows_s = R[~is0], rows[~is0]
        dm = Rs["Dmeson_delta_mass"].astype(np.float64)
        ks = (dm > 0) & (dm < DM_MAX) & (Rs["Dmeson_Lxy"] >= a.lxy_min)
        Rs, rows_s, dm = Rs[ks], rows_s[ks], dm[ks]
        w0 = np.abs(m - CF.M0_D0) < M_WIN
        i0 = np.nonzero(w0)[0]
        u0 = R0["event_uid"][i0]
        o = np.argsort(u0, kind="mergesort")
        i0, u0 = i0[o], u0[o]
        lo, hi = np.searchsorted(u0, Rs["event_uid"], "left"), np.searchsorted(u0, Rs["event_uid"], "right")
        cnt = hi - lo
        js = np.repeat(np.arange(len(Rs)), cnt)
        j0 = i0[np.repeat(lo, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))] if cnt.sum() else np.zeros(0, int)
        m0f = Rs["Dmeson_mass"].astype(np.float64)[js] - dm[js]
        pt0, pts = R0["Dmeson_pt"].astype(np.float64)[j0], Rs["Dmeson_pt"].astype(np.float64)[js]
        good = (np.abs(m[j0] - m0f) < TOL_M) & (pt0 <= pts * 1.001) & (pt0 > 0.8 * pts) & \
               (np.abs(R0["Dmeson_eta"][j0] - Rs["Dmeson_eta"][js]) < 0.1)
        j0, js = j0[good], js[good]
        Pp = table(R0[j0], rows0[j0], X[j0], ev, tp, a.tag, fi, "p_")
        pis = CF.daughter_pt_gev(Rs, 2)
        Pp.update({"p_dm": dm[js].astype(np.float32), "p_row_ds": rows_s[js], "p_pt_ds": pts[good].astype(np.float32),
                   "p_truth_ds": (Rs["Dmeson_truth_category"][js].astype(np.int16) if a.tag != "data" else np.full(len(js), -1, np.int16)),
                   "p_pis_pt": pis[js].astype(np.float32)})
        for k, v in {**T, **Pp}.items():
            parts.setdefault(k, []).append(v)
        stat.append({"file": os.path.basename(path), "rows_read": int(nread), "d0_presel": int(len(R0)), "a": int(aw.sum()), "dstar": int(len(Rs)),
                     "pairs": int(len(j0))})
        print(f"[{a.tag}] {os.path.basename(path)}: wiersze {nread}  D0(presel) {len(R0)}  a {int(aw.sum())}  D* {len(Rs)}  pary {len(j0)}  "
              f"{time.time() - t0:.0f}s", flush=True)
    keys = [k for k in parts if len(parts[k]) == len(a.files)]   # kolumny obecne we wszystkich plikach
    out = {k: np.concatenate(parts[k]) for k in keys}
    out["meta"] = np.array(json.dumps({"tag": a.tag, "files": stat, "feats": CF.FEATS, "presel": P, "lxy_min": a.lxy_min,
                                       "all_window": [a.all_lo, a.all_hi], "all_frac": a.all_frac, "fd_bkg": bool(a.fd_bkg),
                                       "dm_max": DM_MAX, "m_win": M_WIN, "event_cols": sorted(evk), "trig": sorted(tpk),
                                       "n_a": int(len(out.get("a_m", []))), "n_p": int(len(out.get("p_m", []))), "cnclone": "1.0"}))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.tag}] zapisano {a.out}: a {len(out.get('a_m', []))}, pary {len(out.get('p_m', []))}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
