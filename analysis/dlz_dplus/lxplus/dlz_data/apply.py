"""DLZ_dplus/DATA: ekstrakcja D+ z h5 (v2.0/v1.2) + score wszystkich modeli (fold po (run,event) → model foldu, który
tego zdarzenia nie widział) → slim npz. Ten sam skrypt lokalnie i w paczce lxplus (Q_lxplus/dlz_data).

Użycie: python apply.py --models MDIR --thr THR.json --out OUT.npz [--tag data|mc] [--nthread 3] plik1.h5 plik2.h5 ...
Wyjście: wiersze „luźne” (po weto D*, L_xy ≥ 0,25, score ≥ próg WP95 własnego binu pT dla co najmniej jednego modelu seed 12345):
  m, pt, run, event, fold, grl, mbts, j12, nselbl, s_<wariant>_s<seed>; dla wierszy z s_H ≥ WP95(H): indeks `fi` + 26 cech `feat`;
  hall_any / hall_mb: histogramy masy WSZYSTKICH kandydatów po weto (3 biny pT × 500 binów 1 MeV; mb = GRL ∧ HLT_mb_mbts_L1MBTS_1).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dlzlib as D  # noqa: E402

SETS = [("H", 12345), ("HDLZ", 12345), ("Hdal", 12345), ("H", 23456), ("HDLZ", 23456)]
FEATS = D.H + D.DLZ8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True)
    ap.add_argument("--thr", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nthread", type=int, default=3)
    ap.add_argument("--lxy-min", type=float, default=0.25)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    T = json.load(open(a.thr))
    sets = [(v, s) for v, s in SETS if os.path.isdir(os.path.join(a.models, f"{v}_s{s}"))]
    t95 = {f"{v}_s{s}": np.array([T[f"{v}_s{s}"][f"{lo}-{hi}"]["0.95"] for lo, hi in D.BINS]) for v, s in sets if s == 12345}
    models = {(v, s): D.load_models(a.models, v, s) for v, s in sets}
    acc = {}
    n_in = n_veto = 0
    E1 = np.linspace(D.SP["mass_range"][0], D.SP["mass_range"][1], 501)  # 1 MeV
    hall = {"any": np.zeros((3, 500), np.int64), "mb": np.zeros((3, 500), np.int64)}
    for path in a.files:
        L = D.read_h5(path, a.lxy_min)
        n_in += len(L)
        L = L[L["dstar_ok"].to_numpy(bool)].reset_index(drop=True)
        n_veto += len(L)
        if not len(L):
            continue
        fold = D.data_fold(L)
        o = {"m": L["m"].to_numpy(np.float64), "pt": L["pt"].to_numpy(np.float32), "run": L["run"].to_numpy(np.int64),
             "event": L["event"].to_numpy().astype(np.int64), "fold": fold}
        for nm, c in (("grl", "passGRL_MinBiasAFP"), ("mbts", "trigPassed_HLT_mb_mbts_L1MBTS_1"), ("j12", "trigPassed_HLT_noalg_L1J12"),
                      ("nselbl", "nSelBL")):
            o[nm] = L[c].to_numpy().astype(np.int16) if c in L else np.full(len(L), -1, np.int16)
        ib = D.ptbin(o["pt"])
        mbc = (o["grl"] == 1) & (o["mbts"] == 1)
        for bi in range(3):
            hall["any"][bi] += np.histogram(o["m"][ib == bi], E1)[0]
            hall["mb"][bi] += np.histogram(o["m"][(ib == bi) & mbc], E1)[0]
        loose = np.zeros(len(L), bool)
        for (v, s), mm in models.items():
            o[f"s_{v}_s{s}"] = D.predict(mm, L, D.VARIANTS[v], fold, a.nthread)
            if s == 12345:
                loose |= o[f"s_{v}_s{s}"] >= t95[f"{v}_s{s}"][ib]
        o = {k: v[loose] for k, v in o.items()}
        keep = o["s_H_s12345"] >= t95["H_s12345"][ib[loose]]
        o["fi_local"] = keep
        o["feat"] = L.loc[np.flatnonzero(loose)[keep], FEATS].to_numpy(np.float32)
        for k, v in o.items():
            acc.setdefault(k, []).append(v)
        print(f"{os.path.basename(path)}: {len(L)} kand. po weto, luźnych {int(loose.sum())}, z cechami {int(keep.sum())}", flush=True)
    out = {k: np.concatenate(v) for k, v in acc.items()}
    out["fi"] = np.flatnonzero(out.pop("fi_local"))
    out["feat_names"] = np.array(FEATS)
    out["n_files"] = len(a.files)
    out["n_before_veto"] = n_in
    out["n_after_veto"] = n_veto
    out["hall_any"], out["hall_mb"] = hall["any"], hall["mb"]
    np.savez_compressed(a.out, **out)
    print("zapisano", a.out, len(out["m"]), flush=True)


if __name__ == "__main__":
    main()
