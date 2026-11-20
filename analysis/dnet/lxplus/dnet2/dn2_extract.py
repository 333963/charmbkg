#!/usr/bin/env python3
"""DNET etap 2 — tory z otoczenia kandydata D (drzewo `Tracks` ROOT v2.0) dla wierszy tabel etapu 1 (dn_<próbka>.npz; bez czytania h5).
Dla kandydata: tory zdarzenia (wpis `entry` pliku ROOT o numerze `ftag`) w ΔR < 1,0 od kierunku D, BEZ córek (ΔR(tor, córka) < 0,01 i
|pT_tor/pT_córki − 1| < 0,1), NT o największym pT → E [n, NT, 6] = (pT [GeV], Δη, Δφ, ładunek, dE/dx, npts_dEdx), float16; nE = liczba torów
w ΔR < 0,4 / 1,0; pc40 = Σ pT (pT > 0,5 GeV, ΔR < 0,4) — kontrola względem gałęzi Dmeson_ptcone40. Drzewo Tracks NIE ma d0/z0.
Wybór wierszy: a_ — podpróbka zdarzeń (hash event_uid) --frac-a; p_ (pary D0←D*) — okno |m − m0| < 55 MeV, 100 %.
Użycie: dn2_extract.py --sample S --data DIR --rawdir DIR --part k --nparts N --out OUT.npz [--species D0] [--frac-a f] [--max-files n]"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

NT, RMAX = 20, 1.0
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("DN1") or os.path.join(os.path.dirname(HERE), "dnet1"))
BR = ["pt", "eta", "phi", "charge", "dEdx", "npts_dEdx"]


def neighbours(tr, off, pev, dax, dau):
    """tr: dict płaskich tablic torów porcji; off: offsety zdarzeń; pev: zdarzenie kandydata; dax [n, 2] (η, φ) D; dau [n, nd, 3] (pT, η, φ)."""
    n = len(pev)
    E = np.zeros((n, NT, 6), np.float16); nE = np.zeros((n, 2), np.int16); pc = np.zeros(n, np.float32)
    cnt = (off[pev + 1] - off[pev]).astype(np.int64)
    if cnt.sum() == 0:
        return E, nE, pc
    ci = np.repeat(np.arange(n), cnt)
    ti = np.repeat(off[pev], cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    pt, eta, phi = tr["pt"][ti], tr["eta"][ti], tr["phi"][ti]
    deta = eta - dax[ci, 0]
    dphi = (phi - dax[ci, 1] + np.pi) % (2 * np.pi) - np.pi
    keep = np.hypot(deta, dphi) < RMAX
    for d in range(dau.shape[1]):
        dd = np.hypot(eta - dau[ci, d, 1], (phi - dau[ci, d, 2] + np.pi) % (2 * np.pi) - np.pi)
        keep &= ~((dd < 0.01) & (np.abs(pt / np.maximum(dau[ci, d, 0], 1e-6) - 1) < 0.1))
    ci, ti, pt, deta, dphi = ci[keep], ti[keep], pt[keep], deta[keep], dphi[keep]
    dr = np.hypot(deta, dphi)
    nE[:, 1] = np.bincount(ci, minlength=n)
    nE[:, 0] = np.bincount(ci[dr < 0.4], minlength=n)
    pc[:] = np.bincount(ci[(dr < 0.4) & (pt > 0.5)], weights=pt[(dr < 0.4) & (pt > 0.5)], minlength=n)
    o = np.lexsort((-pt, ci))
    ci, ti, pt, deta, dphi = ci[o], ti[o], pt[o], deta[o], dphi[o]
    first = np.searchsorted(ci, ci, "left")
    rank = np.arange(len(ci)) - first
    k = rank < NT
    E[ci[k], rank[k]] = np.stack([pt[k], deta[k], dphi[k], tr["charge"][ti[k]], tr["dEdx"][ti[k]], tr["npts_dEdx"][ti[k]]], 1).astype(np.float16)
    return E, nE, pc


def uid_to_ftag(d):
    """event_uid → numer pliku (z nazwy <próbka>_NNNNNN_none.h5); zdarzenie w > 1 pliku = błąd."""
    import re
    import h5py
    U, T = [], []
    fs = sorted(glob.glob(os.path.join(d, "*_none.h5")))
    if not fs:
        raise RuntimeError(f"brak h5 w {d}")
    for p in fs:
        with h5py.File(p, "r") as f:
            u = np.unique(f["events"].fields(["event_uid"])[:]["event_uid"].astype(np.uint64))
        U.append(u); T.append(np.full(len(u), int(re.search(r"_(\d{6})_none", os.path.basename(p)).group(1)), np.int32))
    U, T = np.concatenate(U), np.concatenate(T)
    o = np.argsort(U, kind="stable"); U, T = U[o], T[o]
    dup = int((U[1:] == U[:-1]).sum())
    if dup:
        raise RuntimeError(f"{d}: {dup} event_uid w więcej niż jednym pliku")
    print(f"uid→ftag: {len(fs)} plików h5, {len(U)} zdarzeń", flush=True)
    return U, T


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--data", required=True); ap.add_argument("--rawdir", required=True)
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
    ap.add_argument("--out", required=True); ap.add_argument("--species", default="D0")
    ap.add_argument("--frac-a", type=float, default=1.0); ap.add_argument("--max-files", type=int, default=0)
    ap.add_argument("--step", type=int, default=20000)
    ap.add_argument("--h5dir", default=os.environ.get("H5V20", "/eos/user/k/klyczak/charm_h5_v20/parts"))
    a = ap.parse_args(argv)
    import awkward as ak
    import uproot
    import dn_feat as DF
    t0 = time.time()
    z = np.load(os.path.join(a.data, f"dn_{a.sample}.npz"), allow_pickle=False)
    UT = [None]
    out, stat = {}, {"sample": a.sample, "part": a.part, "nparts": a.nparts, "files": {}, "NT": NT, "RMAX": RMAX, "cols": BR}
    for sp in a.species.split(","):
        m0, win = DF.cfg["species"][sp]["m0"], DF.cfg["species"][sp]["signal_window"]
        for pre in ("a_", "p_"):
            k = f"{sp}__{pre}"
            if k + "m" not in z.files:
                continue
            ftag, entry, uid, m = z[k + "ftag"], z[k + "entry"], z[k + "uid"], z[k + "m"]
            if a.sample.startswith("data"):
                # etap 1 zapisał dla danych ftag = numer runu (regex trafił w „data_<run>_”); numer pliku odtwarzany z event_uid w h5
                if UT[0] is None:
                    UT[0] = uid_to_ftag(os.path.join(a.h5dir, a.sample))
                U, Tg = UT[0]
                j = np.clip(np.searchsorted(U, uid), 0, len(U) - 1)
                bad = int((U[j] != uid).sum())
                if bad:
                    raise RuntimeError(f"{a.sample} {k}: {bad} wierszy bez pliku h5 (event_uid)")
                ftag = Tg[j]
            sel = (DF.unit_hash(uid, salt=31) < a.frac_a) if pre == "a_" else (np.abs(m - m0) < win)
            tags = np.unique(ftag)[a.part::a.nparts]
            if a.max_files:
                tags = tags[:a.max_files]
            idx = np.flatnonzero(sel & np.isin(ftag, tags))
            idx = idx[np.lexsort((entry[idx], ftag[idx]))]
            n = len(idx)
            E = np.zeros((n, NT, 6), np.float16); nE = np.zeros((n, 2), np.int16); pc = np.zeros(n, np.float32)
            dax = np.stack([z[k + "eta"][idx], z[k + "phi"][idx]], 1).astype(np.float64)
            dau = z[k + "Dd"][idx][:, :, :3].astype(np.float64)
            ft, en = ftag[idx], entry[idx]
            for tag in tags:
                tf = time.time()
                lo, hi = np.searchsorted(ft, [tag, tag + 1])
                if hi == lo:
                    continue
                rf = [f for f in glob.glob(os.path.join(a.rawdir, f"*._{tag:06d}.*.root"))]
                if len(rf) > 1:
                    rf = [f for f in rf if "dpt3to150" in f]
                if len(rf) != 1:
                    raise RuntimeError(f"{a.sample} {tag:06d}: plików ROOT {len(rf)} (oczekiwany 1)")
                ents = en[lo:hi]
                ue = np.unique(ents)
                with uproot.open(rf[0]) as R:
                    T = R["Tracks"]
                    if ue[-1] >= T.num_entries or ue[0] < 0:
                        raise RuntimeError(f"{rf[0]}: wpis {ue[-1]} poza {T.num_entries}")
                    for s in range(0, T.num_entries, a.step):
                        e = min(T.num_entries, s + a.step)
                        j0, j1 = np.searchsorted(ue, [s, e])
                        if j1 == j0:
                            continue
                        arr = T.arrays([f"Tracks_{b}_NOSYS" for b in BR], entry_start=s, entry_stop=e)[ue[j0:j1] - s]
                        cnt = ak.to_numpy(ak.num(arr["Tracks_pt_NOSYS"]))
                        off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
                        tr = {b: ak.to_numpy(ak.flatten(arr[f"Tracks_{b}_NOSYS"])).astype(np.float64) for b in BR}
                        tr["pt"] *= 1e-3
                        p0, p1 = np.searchsorted(ents, [s, e])
                        pev = np.searchsorted(ue[j0:j1], ents[p0:p1])
                        sl = slice(lo + p0, lo + p1)
                        E[sl], nE[sl], pc[sl] = neighbours(tr, off, pev, dax[sl], dau[sl])
                stat["files"][f"{sp}{pre}{tag:06d}"] = [int(hi - lo), round(time.time() - tf, 1)]
            out[k + "idx"], out[k + "E"], out[k + "nE"], out[k + "pc40"] = idx, E, nE, pc
            if n and "iso" in "".join(DF.ALLF[sp]):
                ref = z[k + "X"][idx][:, DF.ALLF[sp].index("ptcone40_iso")].astype(np.float64)
                d = pc - ref
                stat[f"{sp}{pre}ptcone40_check"] = {"n": int(n), "median_diff": float(np.median(d)), "frac_abs_lt_50MeV": float(np.mean(np.abs(d) < 0.05)),
                                                   "corr": float(np.corrcoef(pc, ref)[0, 1]) if n > 10 else None}
            print(f"[{a.sample}] {sp} {pre}: wierszy {n} z {len(m)}, plików {len(tags)}, {time.time() - t0:.0f}s  {stat.get(f'{sp}{pre}ptcone40_check')}", flush=True)
    out["meta"] = np.array(json.dumps(stat))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample}] zapisano {a.out}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
