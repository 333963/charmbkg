#!/usr/bin/env python3
"""DLZ_DATA: ramki treningowe D+ z h5 v2.0 → npz (cechy H + 8 DLZ, m, pT, |η|, waga, truth, klucze).
  --tag fd : sygnał (truth 0/1) po weto D*            --tag mb : wszystko po weto D* (rzeźbienie, closure)
  --tag sb : dane, SB |m−m0| > 70 MeV po weto D*, podpróbka zdarzeń --frac (hasz (run,event))
Użycie: dlz_frames.py --tag fd|mb|sb [--frac 0.06] --out X.npz plik.h5 ..."""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dlzlib as D  # noqa: E402

FEATS = D.H + D.DLZ8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True, choices=["fd", "mb", "sb"])
    ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--lxy-min", type=float, default=0.25)
    ap.add_argument("--out", required=True)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    acc = {}
    n_in = 0
    for path in a.files:
        L = D.read_h5(path, a.lxy_min)
        n_in += len(L)
        keep = L["dstar_ok"].to_numpy(bool).copy()
        if a.tag == "fd":
            keep &= np.isin(L["truth"].to_numpy(), [0, 1])
        if a.tag == "sb":
            keep &= np.abs(L["m"].to_numpy() - D.M0) > D.SP["sideband_min"]
            if a.frac < 1:
                h = D.fold_of(L["run"].to_numpy(np.int64) + 7, L["event"].to_numpy().astype(np.int64), 10000)  # niezależny od foldów (inne ziarno)
                keep &= h.astype(np.int64) % 10000 < int(round(a.frac * 10000))
        L = L[keep].reset_index(drop=True)
        o = {"X": L[FEATS].to_numpy(np.float32), "m": L["m"].to_numpy(np.float64), "pt": L["pt"].to_numpy(np.float32),
             "abs_eta": L["abs_eta"].to_numpy(np.float32), "weight": L["weight"].to_numpy(np.float32),
             "truth": L["truth"].to_numpy(np.int32) if "truth" in L else np.zeros(len(L), np.int32),
             "run": L["run"].to_numpy(np.int64), "event": L["event"].to_numpy().astype(np.int64), "evkey": L["evkey"].to_numpy(np.int64)}
        for k, v in o.items():
            acc.setdefault(k, []).append(v)
        print(f"{os.path.basename(path)}: {len(L)} / {int(keep.size)}", flush=True)
    out = {k: np.concatenate(v) for k, v in acc.items()}
    out["names"] = np.array(FEATS)
    out["n_files"] = len(a.files)
    out["n_in"] = n_in
    np.savez_compressed(a.out, **out)
    print("zapisano", a.out, len(out["m"]), flush=True)


if __name__ == "__main__":
    main()
