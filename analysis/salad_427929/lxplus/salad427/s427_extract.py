#!/usr/bin/env python3
"""SALAD427 — ekstrakcja kandydatów D0 i D+ z h5 v2.0 (compound `events`) z 18 cechami H (charmbkg.features.compute).

Selekcja (jak B1 / nominał): decayType gatunku, pT 3-30 GeV, okno masy gatunku (D0 1,75-2,00; D+ 1,70-2,20), L_xy >= --lxy-min
(0,25 mm), D+: weto D* (charmbkg.samples.dstar_veto_mask). D0: weto refleksji NIE jest tu stosowane, zapisywane jest m_swap
(|m_swap - m0| > 25 MeV tnie s427_study.py). MC: truth z Dmeson_truth_category; D0 z hypo_swap == 1 i truth 0/1 -> truth += 100.
  --tag data : wszystkie zdarzenia; dla 428770 klasa J12only (J12 i brak MBTS) podpróbkowana po zdarzeniach (--frac-j12, hash
               event_uid), reszta 100 %; `sw` = 1/frac. Dla 427929 tylko zdarzenia z GRL.
  --tag mb   : wszystko (sygnał truth 0/1 zostaje, odcina go study);   --tag fd : tylko truth 0/1 (+ swapy 100/101).
Klucz zdarzenia: event_uid (dane i MC). Wyjście npz: <sp>_X [n,18] f4, _m f8, _pt, _truth, _uid, _mswap (D0), _tbits (bity
łańcuchów wg s427_lib.CHAINS), _grl, _nsel, _sw; meta (json).
Użycie: s427_extract.py --tag data|mb|fd --sample NAZWA --out OUT.npz [--frac-j12 0.25] [--max-rows N] plik.h5 ..."""
import argparse
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s427_lib as SL  # noqa: E402

SL.add_src()
from charmbkg import config as C, features as F  # noqa: E402
from charmbkg.bkggen import salad_v1 as SV  # noqa: E402
from charmbkg.kinematics import swapped_mass_sp  # noqa: E402
from charmbkg.samples import colmap, dstar_veto_mask  # noqa: E402

PFX = {"D0": "D0", "Dplus": "Dp"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, choices=["data", "mb", "fd"])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lxy-min", type=float, default=SL.LXY_MIN)
    ap.add_argument("--frac-j12", type=float, default=0.25)
    ap.add_argument("--chunk", type=int, default=400_000)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    cfg = C.load()
    scfg = SV.load_cfg()
    feats = list(scfg["features"])
    units = cfg.get("units")
    acc = {s: {} for s in SL.SPECIES}
    stat, t0 = [], time.time()
    for path in a.files:
        cm = colmap(cfg, path, 3)                       # mapa dla 3 córek zawiera też 2
        with h5py.File(path, "r") as f:
            ds = f["events"]
            names = set(ds.dtype.names)
            tcols = {c: "trigPassed_" + c for c in SL.CHAINS if "trigPassed_" + c in names}
            extra = [c for c in (SL.UID_COL, SL.GRL_COL, "nSelBL") if c in names]
            if SL.UID_COL not in names:
                raise SystemExit(f"{path}: brak {SL.UID_COL}")
            if a.tag == "data" and len(tcols) < len(SL.CHAINS):
                raise SystemExit(f"{path}: brak kolumn trigPassed_ dla {[c for c in SL.CHAINS if c not in tcols]}")
            cols = [c for c in dict.fromkeys(F.raw_columns(cm, 3) + extra + list(tcols.values())) if c in names]
            n = ds.shape[0] if not a.max_rows else min(ds.shape[0], a.max_rows)
            cnt = {s: 0 for s in SL.SPECIES}
            for lo in range(0, n, a.chunk):
                x = ds.fields(cols)[lo:min(n, lo + a.chunk)]
                dt = np.asarray(x[cm["decay_type"]])
                for s in SL.SPECIES:
                    sp = cfg["species"][s]
                    nd = int(sp["n_daughters"])
                    k = np.flatnonzero(dt == int(sp["decay_type"]))
                    if not len(k):
                        continue
                    raw = pd.DataFrame({c: np.asarray(x[c])[k] for c in cols})
                    L = F.to_logical(raw, cm, nd)
                    L["weight"] = L["weight"].fillna(1.0)
                    L = F.compute(L, {**sp, "n_daughters": nd}, units)
                    m, pt = L["m"].to_numpy(float), L["pt"].to_numpy(float)
                    keep = (pt >= 3) & (pt < 30) & (m >= sp["mass_range"][0]) & (m < sp["mass_range"][1]) \
                        & (L["Lxy"].to_numpy(float) >= a.lxy_min)
                    if s == "Dplus":
                        keep &= dstar_veto_mask(L, sp)
                    truth = L["truth"].to_numpy(float)
                    truth = np.where(np.isfinite(truth), truth, -1).astype(np.int16)
                    if a.tag == "data":
                        truth[:] = -1
                    elif s == "D0" and "hypo_swap" in L:
                        sw_ = np.isin(truth, [0, 1]) & (L["hypo_swap"].to_numpy(float) == 1)
                        truth[sw_] += 100
                    if a.tag == "fd":
                        keep &= np.isin(truth, [0, 1, 100, 101])
                    uid = raw[SL.UID_COL].to_numpy().astype(np.uint64)
                    tb = np.zeros(len(raw), np.int32)
                    for i, c in enumerate(SL.CHAINS):
                        if c in tcols:
                            tb |= (raw[tcols[c]].to_numpy() == 1).astype(np.int32) << i     # MC / brak gałęzi = -1 -> 0
                    grl = raw[SL.GRL_COL].to_numpy().astype(np.int8) if SL.GRL_COL in raw else np.full(len(raw), -1, np.int8)
                    sw = np.ones(len(raw), np.float32)
                    if a.tag == "data" and "427929" in a.sample:
                        keep &= grl == 1
                    if a.tag == "data" and "428770" in a.sample and a.frac_j12 < 1:
                        cm_ = SL.class_masks(tb, grl, "data_428770")
                        j = cm_["J12only"]
                        keep &= ~j | (SL.unit_hash(uid, salt=3) < a.frac_j12)
                        sw[j] = 1.0 / a.frac_j12
                    L, raw_k = L[keep], keep
                    o = {"X": L[feats].to_numpy(np.float32), "m": L["m"].to_numpy(np.float64), "pt": L["pt"].to_numpy(np.float32),
                         "truth": truth[raw_k], "uid": uid[raw_k], "tbits": tb[raw_k], "grl": grl[raw_k], "sw": sw[raw_k],
                         "nsel": (raw["nSelBL"].to_numpy()[raw_k].astype(np.int32) if "nSelBL" in raw else np.full(int(raw_k.sum()), -1, np.int32))}
                    if s == "D0":
                        o["mswap"] = np.asarray(swapped_mass_sp(L, sp), np.float32)
                    for kk, v in o.items():
                        acc[s].setdefault(kk, []).append(v)
                    cnt[s] += int(raw_k.sum())
        stat.append({"file": os.path.basename(path), "rows": int(n), **cnt})
        print(f"[{a.sample}] {os.path.basename(path)}: wiersze {n}  D0 {cnt['D0']}  D+ {cnt['Dplus']}  {time.time() - t0:.0f}s", flush=True)
    out = {}
    for s in SL.SPECIES:
        for kk, v in acc[s].items():
            out[f"{PFX[s]}_{kk}"] = np.concatenate(v)
    out["meta"] = np.array(json.dumps({"sample": a.sample, "tag": a.tag, "files": stat, "feats": feats, "lxy_min": a.lxy_min,
                                       "frac_j12": a.frac_j12, "chains": SL.CHAINS, "version": SL.VERSION,
                                       "n": {s: int(len(out.get(f"{PFX[s]}_m", []))) for s in SL.SPECIES}}))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"[{a.sample}] zapisano {a.out}: " + ", ".join(f"{s} {len(out.get(PFX[s] + '_m', []))}" for s in SL.SPECIES)
          + f"  {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
