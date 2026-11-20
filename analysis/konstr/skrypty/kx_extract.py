#!/usr/bin/env python3
"""KONSTR — ekstrakcja (tylko odczyt wejść). Per próbka:
  G|<próbka>|<bin>|<kat>  [nS, 70, 110] uint32: D0, m (5 MeV, 1,70–2,05) × m_swap (5 MeV, 1,60–2,15) × przedział score drzewa H (własny wiersza);
      kat: dane 'all' (427929: MB∧GRL); MC: t0 sygnał poprawny, t1 sygnał zamieniony, t2 truth 2..7, t3 truth ≥ 8.
  --rows (tylko dane z h5): wiersze D0/D+ klasy MB∧GRL∧npu=0 ze score ≥ WP90: sp, m, ms, pt, s, ipt, nsel, uid, t0 (tag A∨C).
Użycie: kx_extract.py --app DIR --wps WPS --sample S --out F.npz [--h5 DIR --tg DIR_z_tg_extract --rows] [--nproc 2]"""
import argparse, glob, json, os, sys, time
import numpy as np
from multiprocessing import Pool
BINS = [(3, 5), (5, 8), (8, 30)]
TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
MR = {"D0": (1.70, 2.05), "Dplus": (1.70, 2.20)}
ME = np.arange(1.70, 2.05 + 1e-9, 0.005); SE = np.arange(1.60, 2.15 + 1e-9, 0.005)
A = THR = EM = None

def init(a):
    global A, THR, EM
    A = a
    W = np.load(a.wps, allow_pickle=False)
    THR = {(sp, b): W[f"{sp}|{TAG[sp]}|{b[0]}-{b[1]}|thr"] for sp in TAG for b in BINS}
    if a.rows:
        sys.path.insert(0, a.tg)
        import tg_extract as T
        EM = T.EvMap(a.h5, a.sample)

def one(f):
    z = np.load(f, allow_pickle=False)
    tags = json.loads(str(z["meta"]))["tags"]
    isdata = A.sample.startswith("data")
    nS = len(THR[("D0", BINS[0])]) + 1
    G, R = {}, {}
    for isp, sp in enumerate(TAG):
        if f"{sp}__m" not in z.files:
            continue
        jt = tags[sp].index(TAG[sp])
        uid, m, pt, ok, cls = (z[f"{sp}__{c}"] for c in ("uid", "m", "pt", "ok", "cls"))
        lg = z[f"{sp}__lg"][:, jt].astype(np.float64)
        ms = z["D0__m_swap"] if sp == "D0" else np.zeros(len(m), np.float32)
        s = np.zeros(len(m), np.int8); ipt = np.full(len(m), -1, np.int8)
        for k, b in enumerate(BINS):
            inb = (pt >= b[0]) & (pt < b[1])
            s[inb] = np.searchsorted(THR[(sp, b)], lg[inb], side="right"); ipt[inb] = k
        base = ok & (ipt >= 0) & (m >= MR[sp][0]) & (m < MR[sp][1])
        if A.sample == "data_427929":
            base &= ((cls & 1) > 0) & ((cls & 2) > 0)
        if sp == "D0":
            if isdata:
                cats = {"all": np.ones(len(m), bool)}
            else:
                truth, swap = z["D0__truth"], z["D0__swap"]
                sg = np.isin(truth, (0, 1))
                cats = {"t0": sg & (swap != 1), "t1": sg & (swap == 1), "t2": ~sg & (truth < 8), "t3": truth >= 8}
            im = np.floor((m.astype(np.float64) - 1.70) / 0.005 + 1e-9).astype(np.int64)
            js = np.floor((ms.astype(np.float64) - 1.60) / 0.005 + 1e-9).astype(np.int64)
            inr = (js >= 0) & (js < 110) & (im >= 0) & (im < 70)
            for k, b in enumerate(BINS):
                for cn, cm in cats.items():
                    q = base & inr & (ipt == k) & cm
                    G[f"G|{A.sample}|{b[0]}-{b[1]}|{cn}"] = np.bincount((s[q].astype(np.int64) * 70 + im[q]) * 110 + js[q], minlength=nS * 7700).astype(np.uint32)
        if A.rows:
            q = base & (s >= 2)                                   # ≥ WP90
            ev = EM.lookup(uid[q])
            g = ev["grl"] & (ev["npu"] == 0)
            R.setdefault("sp", []).append(np.full(int(g.sum()), isp, np.int8))
            for nme, v in (("m", m[q][g]), ("ms", ms[q][g]), ("pt", pt[q][g]), ("s", s[q][g]), ("ipt", ipt[q][g]), ("uid", uid[q][g]), ("nsel", ev["nsel"][g]),
                           ("t0", (ev["tag_A"] | ev["tag_C"])[g])):
                R.setdefault(nme, []).append(v)
    return G, {k: np.concatenate(v) for k, v in R.items()}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    for x in ("app", "wps", "sample", "out"):
        ap.add_argument("--" + x, required=True)
    ap.add_argument("--h5", default=""); ap.add_argument("--tg", default=""); ap.add_argument("--rows", action="store_true"); ap.add_argument("--nproc", type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    fs = sorted(glob.glob(os.path.join(a.app, f"app_{a.sample}_*.npz")))
    G, R = {}, {}
    with Pool(a.nproc, initializer=init, initargs=(a,)) as p:
        for g, r in p.imap_unordered(one, fs, chunksize=2):
            for k, v in g.items():
                G[k] = G.get(k, 0) + v.astype(np.uint64)
            for k, v in r.items():
                R.setdefault(k, []).append(v)
    out = {k: v.reshape(-1, 70, 110).astype(np.uint32) for k, v in G.items()}
    out.update({"R_" + k: np.concatenate(v) for k, v in R.items()})
    np.savez_compressed(a.out, ME=ME, SE=SE, nfiles=len(fs), **out)
    print(a.sample, "pliki", len(fs), "klucze", len(G), "wiersze", len(out.get("R_m", [])), f"{time.time() - t0:.0f}s", flush=True)
