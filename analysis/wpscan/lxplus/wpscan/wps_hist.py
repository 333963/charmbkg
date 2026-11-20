#!/usr/bin/env python3
"""WPSCAN (19.11.2026, projekt B) — histogramy masa × przedział score dla skanu punktu pracy (kryteria: WPSCAN_2026-11-19/KRYTERIA.md).
Wejście: pliki app etapu 4 (dn4_apply.py; app_<próbka>_*.npz, logit float16). Jedno przejście po plikach, tylko numpy.
Progi: kwantyle score sygnału MC fd_nd (prawda 0/1, D0 bez swap, okno |m − m0| < win) dla ε = EPS — jak dn4_ana.MC.
Przedział score s = liczba progów ≤ score (0 = poniżej progu ε = max(EPS), len(EPS) = powyżej progu ε = min(EPS)); WP ε_j przechodzą s ≥ j + 1,
gdzie j = indeks w EPS malejąco. Wyjście npz, klucze <gatunek>|<model>|<bin>|...:
  thr [nE]; eps [nE] (rosnące progi ↔ ε malejące: EPS_DESC)
  fd|sig_w, fd|sig_n, fd|refl_w, fd|refl_n [nS, nm]           (sygnał / odbicie K↔π z fd_nd; _w ważone, _n liczności)
  mb|all [nS, nm], mb|sig [nS, nm], mb|refl [nS, nm]          (pseudo-dane mb_nd: ok & truth < 100; prawdziwy sygnał; odbicie)
  <próbka>|<klasa> [nrep + 1, nS, nm] uint32                  (dane; replika 0 = nominał, 1..nrep = wagi Poissona(1) per zdarzenie)
Bez AFP, tylko pojedyncze D.
Użycie: wps_hist.py --app DIR --out OUT.npz [--nrep 30] [--max-files N]"""
import argparse
import glob
import json
import os
import sys
import time
import zlib

import numpy as np

BINS = [(3, 5), (5, 8), (8, 30)]
EPS_DESC = [round(0.95 - 0.05 * i, 2) for i in range(18)]          # 0,95 … 0,10 (progi rosnące)
SP = {"D0": {"m0": 1.86484, "win": 0.055, "mr": (1.70, 2.05)}, "Dplus": {"m0": 1.86966, "win": 0.045, "mr": (1.70, 2.20)}}
TAGS = {"D0": ["D0.xgb.H.sb.n.s1", "D0.tr.H+e04.mc.d5.s1"], "Dplus": ["Dplus.xgb.H.sb.n.s1", "Dplus.tr.H+e04.mc.d5.s1"]}
MBW = 0.001                                                        # bin masy 1 MeV (fit: 5 MeV; szablony: siatka 2 MeV)
SAMPLES = ["fd_nd", "mb_nd", "data_427929", "data_428770"]


def bname(b):
    return f"{b[0]}-{b[1]}"


def files(appdir, sample, maxf):
    fs = sorted(glob.glob(os.path.join(appdir, f"app_{sample}_*.npz")))
    return fs[:maxf] if maxf else fs


def read(f, sp, need):
    z = np.load(f, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    t = meta["tags"].get(sp) or []
    if not t or f"{sp}__m" not in z.files:
        return None
    T = {c: z[f"{sp}__{c}"] for c in need}
    T["tags"] = t
    return T


def classes(sample, cls):
    grl, mb1, j12, mbts = (cls & 1) > 0, (cls & 2) > 0, (cls & 4) > 0, (cls & 8) > 0
    if sample == "data_427929":
        return {"MB427929": mb1 & grl}
    return {"all": np.ones(len(cls), bool), "J12only": j12 & ~mbts, "notJ12only": ~(j12 & ~mbts)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--app", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--nrep", type=int, default=30); ap.add_argument("--max-files", type=int, default=0)
    a = ap.parse_args(argv)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:6.0f}s] {s}", flush=True)  # noqa: E731
    nE, nS = len(EPS_DESC), len(EPS_DESC) + 1
    out, info = {}, {"eps_desc": EPS_DESC, "bins": BINS, "tags": TAGS, "mbw": MBW, "sp": SP, "nrep": a.nrep, "files": {}, "rows": {}}
    for sp, c in SP.items():
        lo, hi = c["mr"]
        nm = int(round((hi - lo) / MBW))

        def mbin(m):
            ib = np.floor((m.astype(np.float64) - lo) / MBW + 1e-9).astype(np.int64)
            return ib, (ib >= 0) & (ib < nm)

        # ---- 1) fd_nd: progi + szablony
        parts = [r for r in (read(f, sp, ("m", "pt", "ok", "truth", "swap", "w", "lg")) for f in files(a.app, "fd_nd", a.max_files)) if r is not None]
        if not parts:
            log(f"{sp}: brak fd_nd — pomijam"); continue
        tags = parts[0]["tags"]
        F = {k: np.concatenate([p[k] for p in parts]) for k in ("m", "pt", "ok", "truth", "swap", "w", "lg")}
        info["files"][f"{sp}|fd_nd"] = len(parts); info["rows"][f"{sp}|fd_nd"] = int(len(F["m"]))
        jt = {t: tags.index(t) for t in TAGS[sp]}
        thr = {}
        ibF, okF = mbin(F["m"])
        for t in TAGS[sp]:
            for b in BINS:
                inb = (F["pt"] >= b[0]) & (F["pt"] < b[1]) & F["ok"] & np.isin(F["truth"], (0, 1))
                sig = inb & ((F["swap"] != 1) if sp == "D0" else True)
                rf = inb & (F["swap"] == 1) if sp == "D0" else np.zeros(len(inb), bool)
                sw = sig & (np.abs(F["m"] - c["m0"]) < c["win"])
                lg = F["lg"][:, jt[t]].astype(np.float64)
                th = np.array([np.quantile(lg[sw], 1 - e) for e in EPS_DESC]) if sw.sum() >= 200 else np.full(nE, np.inf)
                thr[(t, b)] = th
                s = np.searchsorted(th, lg, side="right")
                key = f"{sp}|{t}|{bname(b)}"
                out[key + "|thr"] = th
                W = float(F["w"][sig].astype(np.float64).sum())
                out[key + "|eps"] = np.array([float(F["w"][sig & (s >= j + 1)].astype(np.float64).sum()) / max(W, 1e-30) for j in range(nE)])
                out[key + "|nsig_tot"] = np.array([int(sig.sum()), int(sw.sum()), int((sig & ~okF).sum())])
                for nmk, msk in (("sig", sig), ("refl", rf)):
                    mm = msk & okF
                    idx = s[mm] * nm + ibF[mm]
                    out[f"{key}|fd|{nmk}_w"] = np.bincount(idx, weights=F["w"][mm].astype(np.float64), minlength=nS * nm).reshape(nS, nm)
                    out[f"{key}|fd|{nmk}_n"] = np.bincount(idx, minlength=nS * nm).reshape(nS, nm).astype(np.uint32)
                log(f"{key}: sygnał MC {int(sig.sum())} (okno {int(sw.sum())}), odbicie {int(rf.sum())}, progi {th[0]:.3f} … {th[-1]:.3f}")
        del F, parts
        # ---- 2) mb_nd (pseudo-dane, domknięcie) i dane
        for sample in SAMPLES[1:]:
            fs = files(a.app, sample, a.max_files)
            isd = sample.startswith("data")
            acc, nrow, nf = {}, 0, 0
            for f in fs:
                T = read(f, sp, ("uid", "m", "pt", "ok", "cls", "truth", "swap", "lg"))
                if T is None:
                    continue
                if T["tags"] != tags:
                    raise RuntimeError(f"{f}: modele {T['tags']} ≠ fd_nd {tags}")
                nf += 1; nrow += len(T["m"])
                ib, okm = mbin(T["m"])
                base = T["ok"] & okm
                if isd:
                    rng = np.random.default_rng(zlib.crc32(os.path.basename(f).encode()) + 1000003 * (1 if sp == "D0" else 2))
                    u, inv = np.unique(T["uid"], return_inverse=True)
                    wr = rng.poisson(1.0, (len(u), a.nrep)).astype(np.uint8)
                    cl = classes(sample, T["cls"])
                else:
                    base = base & (T["truth"] < 100)
                    tsig = np.isin(T["truth"], (0, 1)) & ((T["swap"] != 1) if sp == "D0" else True)
                    trf = np.isin(T["truth"], (0, 1)) & (T["swap"] == 1) if sp == "D0" else np.zeros(len(base), bool)
                    cl = {"all": np.ones(len(base), bool), "sig": tsig, "refl": trf}
                for b in BINS:
                    inb = base & (T["pt"] >= b[0]) & (T["pt"] < b[1])
                    if not inb.any():
                        continue
                    for t in TAGS[sp]:
                        s = np.searchsorted(thr[(t, b)], T["lg"][:, jt[t]].astype(np.float64), side="right")
                        for cn, cm in cl.items():
                            mm = inb & cm
                            idx = s[mm] * nm + ib[mm]
                            key = f"{sp}|{t}|{bname(b)}|{sample[5:] if isd else 'mb'}|{cn}"
                            if isd:
                                h = acc.setdefault(key, np.zeros((a.nrep + 1, nS * nm), np.uint32))
                                h[0] += np.bincount(idx, minlength=nS * nm).astype(np.uint32)
                                wm = wr[inv[mm]]
                                for r in range(a.nrep):
                                    h[r + 1] += np.bincount(idx, weights=wm[:, r], minlength=nS * nm).astype(np.uint32)
                            else:
                                h = acc.setdefault(key, np.zeros(nS * nm, np.uint32))
                                h += np.bincount(idx, minlength=nS * nm).astype(np.uint32)
            for key, h in acc.items():
                out[key] = h.reshape((-1, nS, nm)) if isd else h.reshape(nS, nm)
            info["files"][f"{sp}|{sample}"] = nf; info["rows"][f"{sp}|{sample}"] = int(nrow)
            log(f"{sp} {sample}: {nf} plików (z {len(fs)}), {nrow} wierszy")
    info["sec"] = round(time.time() - t0, 1)
    out["meta"] = np.array(json.dumps(info))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    log(f"zapisano {a.out} ({os.path.getsize(a.out) / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
