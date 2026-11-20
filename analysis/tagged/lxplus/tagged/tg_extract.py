#!/usr/bin/env python3
"""TAGGED (19.11.2026, projekt B) — masa × score × tag AFP per kandydat (kryteria: TAGGED_2026-11-19/KRYTERIA.md, zapisane przed ekstrakcją).
Wejście (tylko odczyt): pliki app etapu 4 (app_<próbka>_*.npz: uid, m, pt, ok, cls, lg), evtable h5 v2.0 (mapa event_uid → zdarzenie: numer, GRL, npu,
tag AFP, ξ, ξ_cal, bity 11 łańcuchów), progi WP z wps_hist.npz (drzewo H). Wyjście per run:
  tg_cand_<run>.npz  kandydaci D0 / D+ w zdarzeniach GRL ∧ (tag A ∨ tag C) [T0] — kolumny per wiersz (bez selekcji npu / ξ_cal / score);
  tg_incl_<run>.npz  histogramy inkluzywne masa (1 MeV) × przedział score, klucze <gatunek>|<bin>|<klasa> [nS, nm] uint32
                     (klasy 427929: MB, MB_npu0, MBAFP_npu0; 428770: all, all_npu0, J12_npu0) + meta.
Program nie drukuje żadnych liczb o masie w zdarzeniach z tagiem (tylko liczby wierszy).
Użycie: tg_extract.py --app DIR --h5 DIR --wps WPS.npz --out DIR [--runs 427929,428770] [--max-files N]"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

BINS = [(3, 5), (5, 8), (8, 30)]
SP = {"D0": {"mr": (1.70, 2.05), "tag": "D0.xgb.H.sb.n.s1"}, "Dplus": {"mr": (1.70, 2.20), "tag": "Dplus.xgb.H.sb.n.s1"}}
MBW = 0.001
TRIG = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1MBTS_2",
        "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2"]
EVC = ["event_uid", "event_number", "lumiBlock", "passGRL_MinBiasAFP", "nPileupVertices", "nSelBL", "xi_cal_A", "xi_cal_C", "afp_nProtons_A", "afp_nProtons_C",
       "afp_xi1_A", "afp_xi1_C", "afp_tag_A", "afp_tag_C"] + ["trigPassed_" + t for t in TRIG]


def fileno(uid):
    return int((int(uid) >> 22) & 0x3FFF)


class EvMap:
    def __init__(self, h5dir, sample):
        import h5py
        self.h5py = h5py
        self.paths = {}
        for f in sorted(glob.glob(os.path.join(h5dir, sample, f"{sample}_*_none.h5"))):
            with h5py.File(f, "r") as h:
                u = h["evtable"].fields("event_uid")[0:1]
            self.paths[(int(u[0]) >> 22)] = f                      # klucz: (task, nr pliku)
        self.cache = {}

    def get(self, key):
        if key not in self.cache:
            if len(self.cache) > 40:
                self.cache.clear()
            with self.h5py.File(self.paths[key], "r") as h:
                e = h["evtable"].fields(EVC)[:]
            e = e[np.argsort(e["event_uid"])]
            tw = np.zeros(len(e), np.uint16)
            for i, t in enumerate(TRIG):
                tw |= ((e["trigPassed_" + t] > 0).astype(np.uint16) << i)
            self.cache[key] = (e, tw)
        return self.cache[key]

    def lookup(self, uid):
        """→ dict kolumn zdarzenia per wiersz."""
        n = len(uid)
        out = {"event_number": np.zeros(n, np.uint64), "lb": np.zeros(n, np.uint16), "grl": np.zeros(n, bool), "npu": np.zeros(n, np.int16), "nsel": np.zeros(n, np.int16),
               "xical_A": np.zeros(n, np.float32), "xical_C": np.zeros(n, np.float32), "np_A": np.zeros(n, np.int16), "np_C": np.zeros(n, np.int16),
               "xi_A": np.zeros(n, np.float32), "xi_C": np.zeros(n, np.float32), "tag_A": np.zeros(n, bool), "tag_C": np.zeros(n, bool), "trig": np.zeros(n, np.uint16)}
        k = (uid >> np.uint64(22)).astype(np.int64)
        for key in np.unique(k):
            m = k == key
            e, tw = self.get(int(key))
            ix = np.searchsorted(e["event_uid"], uid[m])
            if ix.max(initial=0) >= len(e) or not np.array_equal(e["event_uid"][ix], uid[m]):
                raise RuntimeError(f"uid spoza evtable (klucz {key})")
            r = e[ix]
            out["event_number"][m] = r["event_number"]; out["lb"][m] = r["lumiBlock"]; out["grl"][m] = r["passGRL_MinBiasAFP"] > 0
            out["npu"][m] = r["nPileupVertices"]; out["nsel"][m] = np.clip(r["nSelBL"], -1, 32000)
            out["xical_A"][m] = r["xi_cal_A"]; out["xical_C"][m] = r["xi_cal_C"]; out["np_A"][m] = r["afp_nProtons_A"]; out["np_C"][m] = r["afp_nProtons_C"]
            out["xi_A"][m] = r["afp_xi1_A"]; out["xi_C"][m] = r["afp_xi1_C"]; out["tag_A"][m] = r["afp_tag_A"] > 0; out["tag_C"][m] = r["afp_tag_C"] > 0
            out["trig"][m] = tw[ix]
        return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--app", required=True); ap.add_argument("--h5", required=True); ap.add_argument("--wps", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--runs", default="427929,428770"); ap.add_argument("--max-files", type=int, default=0)
    a = ap.parse_args(argv)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    W = np.load(a.wps, allow_pickle=False)
    winfo = json.loads(str(W["meta"]))
    eps = winfo["eps_desc"]
    nS = len(eps) + 1
    thr = {(sp, b): W[f"{sp}|{c['tag']}|{b[0]}-{b[1]}|thr"] for sp, c in SP.items() for b in BINS}
    os.makedirs(a.out, exist_ok=True)
    iMB, iAFP, iJ12 = TRIG.index("HLT_mb_mbts_L1MBTS_1"), TRIG.index("HLT_noalg_L1AFP_A_OR_C"), TRIG.index("HLT_noalg_L1J12")
    for run in a.runs.split(","):
        sample = f"data_{run}"
        EM = EvMap(a.h5, sample)
        log(f"{sample}: {len(EM.paths)} plików h5")
        fs = sorted(glob.glob(os.path.join(a.app, f"app_{sample}_*.npz")))
        if a.max_files:
            fs = fs[:a.max_files]
        H, C, info = {}, {}, {"run": run, "eps_desc": eps, "bins": BINS, "sp": SP, "trig": TRIG, "mbw": MBW, "files": len(fs), "rows": {}, "rows_T0": {}, "cls_mismatch": 0}
        for f in fs:
            z = np.load(f, allow_pickle=False)
            tags = json.loads(str(z["meta"]))["tags"]
            for isp, (sp, c) in enumerate(SP.items()):
                if f"{sp}__m" not in z.files:
                    continue
                jt = tags[sp].index(c["tag"])
                uid, m, pt, ok, cls = z[f"{sp}__uid"], z[f"{sp}__m"], z[f"{sp}__pt"], z[f"{sp}__ok"], z[f"{sp}__cls"]
                lg = z[f"{sp}__lg"][:, jt]
                ev = EM.lookup(uid)
                mb1 = (ev["trig"] >> iMB) & 1 > 0
                info["cls_mismatch"] += int((((cls & 1) > 0) != ev["grl"]).sum() + (((cls & 2) > 0) != mb1).sum())
                lo, hi = c["mr"]
                nm = int(round((hi - lo) / MBW))
                ib = np.floor((m.astype(np.float64) - lo) / MBW + 1e-9).astype(np.int64)
                okm = (ib >= 0) & (ib < nm)
                s = np.zeros(len(m), np.int8); ipt = np.full(len(m), -1, np.int8)
                for k, b in enumerate(BINS):
                    inb = (pt >= b[0]) & (pt < b[1])
                    s[inb] = np.searchsorted(thr[(sp, b)], lg[inb].astype(np.float64), side="right"); ipt[inb] = k
                npu0 = ev["npu"] == 0
                if run == "427929":
                    cl = {"MB": ev["grl"] & mb1, "MB_npu0": ev["grl"] & mb1 & npu0, "MBAFP_npu0": ev["grl"] & mb1 & npu0 & ((ev["trig"] >> iAFP) & 1 > 0)}
                else:
                    cl = {"all": np.ones(len(m), bool), "all_npu0": ev["grl"] & npu0, "J12_npu0": ev["grl"] & npu0 & ((ev["trig"] >> iJ12) & 1 > 0)}
                base = ok & okm & (ipt >= 0)
                for k, b in enumerate(BINS):
                    for cn, cm in cl.items():
                        mm = base & (ipt == k) & cm
                        key = f"{sp}|{b[0]}-{b[1]}|{cn}"
                        h = H.setdefault(key, np.zeros(nS * nm, np.uint32))
                        h += np.bincount(s[mm].astype(np.int64) * nm + ib[mm], minlength=nS * nm).astype(np.uint32)
                t0m = ev["grl"] & (ev["tag_A"] | ev["tag_C"])
                row = {"sp": np.full(int(t0m.sum()), isp, np.int8), "m": m[t0m], "pt": pt[t0m], "ok": ok[t0m], "s": s[t0m], "ipt": ipt[t0m], "lg": lg[t0m], "uid": uid[t0m]}
                row.update({k: v[t0m] for k, v in ev.items()})
                for k, v in row.items():
                    C.setdefault(k, []).append(v)
                info["rows"][sp] = info["rows"].get(sp, 0) + int(len(m)); info["rows_T0"][sp] = info["rows_T0"].get(sp, 0) + int(t0m.sum())
        log(f"{sample}: {len(fs)} plików app, wiersze {info['rows']}, T0 {info['rows_T0']}, niezgodności cls {info['cls_mismatch']}")
        info["sec"] = round(time.time() - t0, 1)
        out = {k: v.reshape(nS, -1) for k, v in H.items()}
        out["meta"] = np.array(json.dumps(info))
        np.savez_compressed(os.path.join(a.out, f"tg_incl_{run}.npz"), **out)
        cand = {k: np.concatenate(v) for k, v in C.items()}
        cand["meta"] = np.array(json.dumps(info))
        np.savez_compressed(os.path.join(a.out, f"tg_cand_{run}.npz"), **cand)
        log(f"zapisano tg_incl_{run}.npz, tg_cand_{run}.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
