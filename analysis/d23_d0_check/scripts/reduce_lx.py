#!/usr/bin/env python3
"""Weryfikacja D0 2-3 GeV: redukcja plikow app (dnet4) na lxplus -> maly npz. TYLKO ODCZYT z EOS.
Uzycie: reduce_lx.py MODE(d23|d150) SAMPLE OUT.npz   (d23: pT 2-3; d150: pT 3-5)
Zapis: hist2D (m 5 MeV x lg 0.05) wszystkich kandydatow klasy (per model, MC per kategoria prawdy) + wiersze po luznym cieciu
(+ score blizniaka K<->pi z tego samego zdarzenia, jesli jest w tabeli)."""
import glob, json, sys, os
import numpy as np
from multiprocessing import Pool
APP = "/eos/user/k/klyczak/charm_q/dnet4/app"
TREE, NET = "D0.xgb.H.sb.n.s1", "D0.tr.H+e04.mc.d5.s1"
ME = np.arange(1.70, 2.05 + 1e-9, 0.005)
LE = np.arange(-12.0, 14.0 + 1e-9, 0.05)
mode, sample = sys.argv[1], sys.argv[2]
PT = (2.0, 3.0) if mode == "d23" else (3.0, 5.0)
LOOSE = {"d23": (1.7646, 2.1875), "d150": (1.0, 1.4)}[mode]      # (drzewo, siec): d23 = progi WP50 z dnet4_d23.json
isdata = sample.startswith("data")

def one(f):
    z = np.load(f, allow_pickle=False)
    if "D0__m" not in z.files:
        return None
    tags = json.loads(str(z["meta"]))["tags"]["D0"]
    jt, jn = tags.index(TREE), tags.index(NET)
    uid, m, pt, ok, cls, truth, swap, w, ms = (z["D0__" + c] for c in ("uid", "m", "pt", "ok", "cls", "truth", "swap", "w", "m_swap"))
    lg = z["D0__lg"][:, [jt, jn]].astype(np.float32)
    sel = ok & (pt >= PT[0]) & (pt < PT[1])
    if sample == "data_427929":
        sel &= ((cls & 1) > 0) & ((cls & 2) > 0)
    uid, m, pt, cls, truth, swap, w, ms, lg = (x[sel] for x in (uid, m, pt, cls, truth, swap, w, ms, lg))
    n = len(m)
    # blizniak: ten sam uid, to samo pT, m_blizniaka = m_swap
    o = np.lexsort((pt, uid))
    tw = np.full((n, 2), np.nan, np.float32)
    has = np.zeros(n, bool)
    uo, po, mo, mso, lgo = uid[o], pt[o], m[o], ms[o], lg[o]
    two = np.full((n, 2), np.nan, np.float32); haso = np.zeros(n, bool)
    for k in (1, -1, 2, -2, 3, -3):
        a = np.arange(n); b = a + k
        v = (b >= 0) & (b < n)
        a, b = a[v], b[v]
        g = (uo[a] == uo[b]) & (np.abs(po[a] - po[b]) < 2e-3) & (np.abs(mo[b] - mso[a]) < 3e-4) & ~haso[a]
        two[a[g]] = lgo[b[g]]; haso[a[g]] = True
    tw[o] = two; has[o] = haso
    if isdata:
        cat = np.zeros(n, np.int8)
    else:   # 0 sygnal poprawna hipoteza, 1 sygnal zamieniony, 2 truth 2..7, 3 truth>=8
        cat = np.where(np.isin(truth, (0, 1)), np.where(swap == 1, 1, 0), np.where(truth >= 8, 3, 2)).astype(np.int8)
    ncat = 1 if isdata else 4
    H = np.zeros((ncat, 2, len(ME) - 1, len(LE) - 1))
    ww = w.astype(np.float64) if not isdata else np.ones(n)
    for c in range(ncat):
        q = cat == c
        for j in range(2):
            H[c, j] = np.histogram2d(m[q], np.clip(lg[q, j], LE[0] + 1e-3, LE[-1] - 1e-3), [ME, LE], weights=ww[q])[0]
    keep = (lg[:, 0] >= LOOSE[0]) | (lg[:, 1] >= LOOSE[1])
    if not isdata:
        keep |= cat < 3
    R = dict(uid=uid[keep], m=m[keep], ms=ms[keep], pt=pt[keep], cls=cls[keep], truth=truth[keep], swap=swap[keep], w=w[keep], lg=lg[keep],
             tw=tw[keep], has=has[keep], cat=cat[keep])
    return H, R, n, int(has.sum())

if __name__ == "__main__":
    fs = sorted(glob.glob(f"{APP}/{mode}/app_{sample}_*.npz"))
    with Pool(6) as p:
        out = [x for x in p.imap(one, fs, chunksize=4) if x is not None]
    H = sum(x[0] for x in out)
    R = {k: np.concatenate([x[1][k] for x in out]) for k in out[0][1]}
    ntot, ntw = sum(x[2] for x in out), sum(x[3] for x in out)
    np.savez_compressed(sys.argv[3], H=H, ME=ME, LE=LE, ntot=ntot, ntw=ntw, nfiles=len(fs), **R)
    print(mode, sample, "pliki", len(fs), "kandydaci", ntot, "z blizniakiem", ntw, "wiersze zapisane", len(R["m"]), flush=True)
