#!/usr/bin/env python3
"""TAGGED krok 3 (PRZED odczytem): kontrola złączenia, fity inkluzywne (kształt), pasma boczne regionu tagowanego (P_tag,SB, rząd tła), zabawki.
Okno sygnału w zdarzeniach z tagiem NIE jest czytane (tg_lib.load_cand(blind=True) usuwa wiersze). Wyjście: out/prep.json, tabele/prep.md."""
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tg_lib as L  # noqa: E402

# pozycje: (klucz, run, region, gatunek, biny pT, WP, selekcja ξ, rola)
ITEMS = []
for sp in L.SPN:
    ITEMS += [(f"{sp}|427929|MBTS|3-30|0.7|T2", "427929", "MBTS", sp, L.BINS, 0.70, "T2", "main" if sp == "D0" else "K1"),
              (f"{sp}|427929|MBTS|5-30|0.7|T2", "427929", "MBTS", sp, L.BINS[1:], 0.70, "T2", "pt5")]
    ITEMS += [(f"{sp}|427929|MBTS|3-30|{e}|T2", "427929", "MBTS", sp, L.BINS, e, "T2", "K3") for e in (0.50, 0.80, 0.90)]
    ITEMS += [(f"{sp}|428770|{r}|3-30|0.7|T2", "428770", r, sp, L.BINS, 0.70, "T2", "reg") for r in ("R1", "R2", "R3", "J12")]
    ITEMS += [(f"{sp}|427929|MBTS|3-30|0.7|CR", "427929", "MBTS", sp, L.BINS, 0.70, "CR", "K5"), (f"{sp}|428770|J12|3-30|0.7|CR", "428770", "J12", sp, L.BINS, 0.70, "CR", "K5")]
    ITEMS += [(f"{sp}|427929|MBAFP|3-30|0.7|T2", "427929", "MBAFP", sp, L.BINS, 0.70, "T2", "K6")]
    ITEMS += [(f"{sp}|428770|R2|{b}|0.7|T2", "428770", "R2", sp, [b], 0.70, "T2", "R2bin") for b in L.BINS]
ITEMS += [(f"D0|427929|MBTS_{s}|3-30|0.7|T2", "427929", f"MBTS_{s}", "D0", L.BINS, 0.70, "T2", "K7") for s in "AC"]

_cache = {}


def data(run, blind):
    k = (run, blind)
    if k not in _cache:
        _cache[k] = (incl(run), L.load_cand(run, blind=blind))
    return _cache[k]


_ifit = {}
_incl = {}


def incl(run):
    if run not in _incl:
        _incl[run] = L.load_incl(run)[0]
    return _incl[run]


def ifit(run, sp, cls, e):
    k = (run, sp, cls, e)
    if k not in _ifit:
        _ifit[k] = L.incl_fit(incl(run), sp, cls, e)
    return _ifit[k]


def build(item, blind=True, rng=None):
    """→ dict: S (kształt), H (hist tagowany 5 MeV; blind: okno puste), msb (maska pasm), fits inkluzywne, liczniki pasm."""
    key, run, reg, sp, bins, e, xsel, role = item
    I, C = data(run, blind)
    cls = L.INCL_CLS[(run, reg)]
    fits = ifit(run, sp, cls, e)
    isp = L.SPN.index(sp)
    ib = [L.BINS.index(b) for b in bins]
    sel = (C["sp"] == isp) & C["base"] & C[xsel] & L.region_mask(C, run, reg) & (C["s"] >= L.jwp(e) + 1) & np.isin(C["ipt"], ib)
    H = L.hist5(C["m"][sel], sp, rng)
    wm = L.winmask(sp, rng)
    S = L.sig_shape(sp, fits, rng=rng, bins=bins)
    # pasma: inkluzywne (ta sama klasa) i tagowane, pełny zakres nominalny
    wm0 = L.winmask(sp)
    Hi = sum(L.V.full5(L.incl_h(I, sp, b, cls, e)) for b in bins)
    nsb_t = float(L.hist5(C["m"][sel], sp)[~wm0].sum())
    nsb_i = float(Hi[~wm0].sum())
    wpt_sb = [float((sel & (C["ipt"] == i) & (np.abs(C["m"] - L.WIN[sp][0]) >= L.WIN[sp][1])).sum()) for i in ib]
    Ninc = sum(f["N"] for f in fits if f["bin"] in bins)
    eN = 1.3 * math.sqrt(sum(f["err"] ** 2 for f in fits if f["bin"] in bins))
    return {"key": key, "run": run, "reg": reg, "sp": sp, "bins": bins, "eps": e, "xsel": xsel, "role": role, "S": S, "H": H, "wm": wm, "fits": fits, "sel": sel,
            "nsb_tag": nsb_t, "nsb_incl": nsb_i, "P_tag_SB": nsb_t / max(nsb_i, 1.0), "wpt_sb": wpt_sb, "N_incl": Ninc, "N_incl_err": eN, "Hincl": Hi, "cls": cls}


def main():
    t0 = time.time()
    os.makedirs(L.OUT, exist_ok=True)
    if os.path.exists(os.path.join(L.OUT, "UNBLINDED")):
        sys.exit("STOP: region już odsłonięty — prep nie jest powtarzany")
    R = {"check": {}, "items": {}, "time": time.strftime("%F %T")}
    # ---- kontrola złączenia: histogramy bez weta npu = wps_hist.npz
    L.V._init()
    for run, cls, src in (("427929", "MB", "427929|MB427929"), ("428770", "all", "428770|all")):
        I, C = data(run, True)
        mx = 0
        for sp in L.SPN:
            for bn in L.BINS:
                a = I[f"{sp}|{bn}|{cls}"].astype(np.int64)
                b = L.V.Z[f"{sp}|{L.V.TAG[sp]}|{bn}|{src}"][0].astype(np.int64)
                mx = max(mx, int(np.abs(a - b).sum()))
        ev = {}
        for isp, sp in enumerate(L.SPN):
            m = (C["sp"] == isp) & C["T2"]
            ev[sp] = int(len(np.unique(C["event_number"][m])))
        R["check"][run] = {"sum_abs_diff_vs_wps": mx, "T2_events_with_cand_outside_window": ev,
                           "found_bits": (float(C["found_bits"].mean()) if "found_bits" in C else None)}
        print(f"[{time.time() - t0:5.0f}s] {run}: |Δ| wobec wps_hist = {mx}; zdarzenia T2 z kandydatem (poza oknem) {ev}", flush=True)
    # ---- pozycje
    jobs = []
    for item in ITEMS:
        B = build(item, blind=True)
        msb = ~B["wm"]
        n, dq, vb = L.bkg_order(B["S"], B["H"], msb)
        npred = B["N_incl"] * B["P_tag_SB"]
        r = {k: B[k] for k in ("run", "reg", "sp", "bins", "eps", "xsel", "role", "nsb_tag", "nsb_incl", "P_tag_SB", "wpt_sb", "N_incl", "N_incl_err", "cls")}
        r.update({"nbkg": n, "dq_order": dq, "bkg_sb": [float(t) for t in vb], "N_pred_SBanchor": npred,
                  "bkg_in_window_extrap": float((L.TFit(B["S"], len(B["H"]), n, mask=msb).bsh(vb[2:]) * vb[1])[B["wm"]].sum()),
                  "incl": [{k: f[k] for k in ("bin", "N", "err", "dmu", "sc", "valid", "chi2", "ndf", "frefl", "kap")} for f in B["fits"]]})
        R["items"][item[0]] = r
        if item[7] == "main":
            pts, nt = [0, 10, 20, 40], 1000
        elif item[7] == "R2bin" and B["nsb_tag"] < 2000:
            pts, nt = [], 0
        else:
            pts, nt = [0, max(round(npred), 1)], 300
        for p in pts:
            if item[7] == "main" and p == 0:
                nt = 2000
            for c in range(8):
                jobs.append((item[0], p, (B["S"], len(B["H"]), n, vb, msb, float(p), nt // 8 + (1 if c < nt % 8 else 0), hash((item[0], p, c)) % (2 ** 31), "exp", None)))
        print(f"[{time.time() - t0:5.0f}s] {item[0]}: N_incl {B['N_incl']:.0f} ± {B['N_incl_err']:.0f}, pasma tag {B['nsb_tag']:.0f} / incl {B['nsb_incl']:.0f}, rząd tła {n} (dq {dq:.1f}), N_pred(pasma) {npred:.1f}", flush=True)
    if os.environ.get("TG_NOTOYS") != "1":
        import zlib
        jobs = [(k, p, a[:7] + (zlib.crc32(f"{k}|{p}|{i}".encode()),) + a[8:]) for i, (k, p, a) in enumerate(jobs)]
        with Pool(7) as pool:
            res = pool.map(L.toy_worker, [j[2] for j in jobs], chunksize=1)
        acc = {}
        for (k, p, _), T in zip(jobs, res):
            acc.setdefault((k, p), []).append(T)
        null = {}
        for (k, p), Ts in acc.items():
            T = np.concatenate(Ts)
            R["items"][k].setdefault("toys", []).append(L.toy_summary(T, p))
            if p == 0:
                null[k] = np.sort(T[:, 4])
        np.savez_compressed(os.path.join(L.OUT, "toys_null_q0.npz"), **{k.replace("|", "__"): v for k, v in null.items()})
    R["sec"] = round(time.time() - t0, 1)
    json.dump(R, open(os.path.join(L.OUT, "prep.json"), "w"), indent=1, ensure_ascii=False)
    print(f"[{time.time() - t0:5.0f}s] zapisano out/prep.json")


if __name__ == "__main__":
    main()
