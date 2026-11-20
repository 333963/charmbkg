#!/usr/bin/env python3
"""DNET — ekstrakcja z h5 v2.0 → npz części (klucze <gatunek>__<a_|p_><kolumna>). Tylko pojedyncze D (D0→Kπ, D+→Kππ); bez par D–D.
  a_*  kandydaci inkluzywni (pT 3–30 GeV, zakres mas gatunku, L_xy ≥ --lxy-min): X (wszystkie cechy dn_feat.ALLF), blok córek Dd,
       m, pt, eta, phi, truth, swap, m_swap, ok (D+: weto D*), uid, entry, ftag, row, bity trygera. Podpróbka ZDARZEŃ (hash event_uid):
       --frac-d0 / --frac-dp; fd: tylko sygnał z prawdą (truth 0/1).
  p_*  (tylko D0) pary D0←D* jak cnclone/B_trig (|m − m0| < 75 MeV, 0 < Δm < 165 MeV, |m_D0 − (m_D* − Δm)| < 1 MeV,
       0,8·pT(D*) < pT(D0) ≤ 1,001·pT(D*), |Δη| < 0,1, L_xy(D*) ≥ --lxy-min), 100 % zdarzeń → pomiar ε_S z tagu D*.
Użycie: dn_extract.py --tag data|fd|mb --out OUT.npz [--frac-d0 f] [--frac-dp f] [--lxy-min 0.25] FILE.h5 [...]"""
import argparse
import json
import os
import sys
import time

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dn_feat as DF  # noqa: E402

DM_MAX, M_WIN, TOL_M = 0.165, 0.075, 0.001
DS_COLS = ["Dmeson_decayType", "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_pt", "Dmeson_eta", "Dmeson_Lxy", "Dmeson_truth_category",
           "Dmeson_daughter2_pt", "event_uid"]


def dstar_rows(path):
    with h5py.File(path, "r") as f:
        ds = f["events"]
        dt = ds.fields(["Dmeson_decayType"])[:]["Dmeson_decayType"]
        idx = np.flatnonzero(dt == 1)
        x = ds.fields(DS_COLS)[:][idx] if len(idx) else ds.fields(DS_COLS)[:0]
    return x, idx


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, choices=["data", "fd", "mb"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--frac-d0", type=float, default=1.0); ap.add_argument("--frac-dp", type=float, default=1.0)
    ap.add_argument("--lxy-min", type=float, default=0.25)
    ap.add_argument("--species", default="D0,Dplus")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    parts, stat, t0 = {}, [], time.time()
    for path in a.files:
        ftag = DF.ftag_of(path)
        with h5py.File(path, "r") as f:
            names = set(f["events"].dtype.names)
        ev = {k: c for k, c in DF.EVENT_COLS.items() if c in names}
        tp = {c: "trigPassed_" + c for c in DF.TRIG_CHAINS if "trigPassed_" + c in names}
        xc = list(ev.values()) + list(tp.values())
        st = {"file": os.path.basename(path)}
        for spn in a.species.split(","):
            L, rows = DF.read_species(path, spn, a.lxy_min, xc)
            m = L["m"].to_numpy(float)
            frac = a.frac_d0 if spn == "D0" else a.frac_dp
            aw = np.ones(len(L), bool)
            if frac < 1.0:
                aw &= DF.unit_hash(L["uid"].to_numpy(np.uint64), salt=3) < frac
            if a.tag == "fd":
                aw &= np.isin(L["truth"].to_numpy(float), DF.SIGNAL_CATS)
            T = DF.table(L[aw], rows[aw], spn, a.tag, ftag, "a_", ev, tp)
            st[spn] = {"presel": int(len(L)), "a": int(aw.sum())}
            if spn == "D0":
                Rs, rows_s = dstar_rows(path)
                dm = Rs["Dmeson_delta_mass"].astype(np.float64)
                ks = (dm > 0) & (dm < DM_MAX) & (Rs["Dmeson_Lxy"] >= a.lxy_min)
                Rs, rows_s, dm = Rs[ks], rows_s[ks], dm[ks]
                i0 = np.nonzero(np.abs(m - DF.cfg["species"]["D0"]["m0"]) < M_WIN)[0]
                u0 = L["uid"].to_numpy(np.uint64)[i0]
                o = np.argsort(u0, kind="mergesort")
                i0, u0 = i0[o], u0[o]
                us = Rs["event_uid"].astype(np.uint64)
                lo, hi = np.searchsorted(u0, us, "left"), np.searchsorted(u0, us, "right")
                cnt = hi - lo
                js = np.repeat(np.arange(len(Rs)), cnt)
                j0 = i0[np.repeat(lo, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))] if cnt.sum() else np.zeros(0, int)
                pt0, eta0 = L["pt"].to_numpy(float), L["eta"].to_numpy(float)
                pts = Rs["Dmeson_pt"].astype(np.float64)[js]
                good = (np.abs(m[j0] - (Rs["Dmeson_mass"].astype(np.float64)[js] - dm[js])) < TOL_M) & (pt0[j0] <= pts * 1.001) & \
                       (pt0[j0] > 0.8 * pts) & (np.abs(eta0[j0] - Rs["Dmeson_eta"][js]) < 0.1)
                j0, js = j0[good], js[good]
                P = DF.table(L.iloc[j0], rows[j0], spn, a.tag, ftag, "p_", ev, tp)
                pis = Rs["Dmeson_daughter2_pt"].astype(np.float64)[js]
                if len(pis) and np.nanmedian(pis) > 100:
                    pis = pis / 1e3
                P.update({"p_dm": dm[js].astype(np.float32), "p_row_ds": rows_s[js], "p_pt_ds": pts[good].astype(np.float32),
                          "p_truth_ds": (Rs["Dmeson_truth_category"][js].astype(np.int16) if a.tag != "data" else np.full(len(js), -1, np.int16)),
                          "p_pis_pt": pis.astype(np.float32)})
                T.update(P)
                st[spn]["pairs"] = int(len(j0))
            for k, v in T.items():
                parts.setdefault(f"{spn}__{k}", []).append(v)
        stat.append(st)
        print(f"[{a.tag}] {st}  {time.time() - t0:.0f}s", flush=True)
    out = {k: np.concatenate(v) for k, v in parts.items() if len(v) == len(a.files)}
    out["meta"] = np.array(json.dumps({"tag": a.tag, "files": stat, "feats": DF.ALLF, "dau": DF.DAU, "lxy_min": a.lxy_min,
                                       "frac": {"D0": a.frac_d0, "Dplus": a.frac_dp}, "dnet": "1.0"}))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.tag}] zapisano {a.out}: " + ", ".join(f"{k}={len(v)}" for k, v in out.items() if k.endswith("_m")) + f"  {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
