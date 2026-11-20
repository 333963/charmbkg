#!/usr/bin/env python3
"""REG428 r1 (20.11.2026, projekt B; lxplus, tylko odczyt wejść, wynik w /tmp/klyczak_reg428) — run 428770 per region triggera. Tylko bity triggera.
  bits   : ROOT v2.0 data_428770 (drzewo CharmAnalysis) → bits/b_<fnum>.npz: bits u2 [wpis], lb u2 [wpis], ev u8 [wpis] + ps.json (preskale per LB)
  hist   : pliki app etapu 4 (score drzewa H per kandydat) + mapa bitów (event_uid = task << 36 | fnum << 22 | wpis) → hist/h_<próbka>_<część>.npz
  dstar  : tabela D* (dn_<próbka>.npz, wiersze p_) + score drzewa H (dnet4/sc) + bity → dstar.npz (wiersze w oknie masy D0; fit lokalnie)
  merge  : suma hist/h_* → reg_hist.npz; bsum: liczniki zdarzeń per bit / LB → bits_summary.json
Bity (u2): 0 GRL, 1 HLT_noalg_L1AFP_A_OR_C, 2 HLT_j20_L1AFP_A_OR_C, 3 HLT_noalg_L1J12, 4 TBP L1_J12, 5 TBP L1_AFP_A_OR_C, 6 HLT_j20_L1J12,
  7 HLT_mb_sptrk_pt4_L1AFP_A_OR_C, 8 HLT_j20f_L1AFP_A_OR_C, 9 HLT_noalg_L1AFP_A_OR_C_J12, 10 HLT_j20_L1AFP_A_OR_C_J12, 11 HLT_j20_L1RD0_FILLED, 12 HLT_j20_L1MBTS_2.
Użycie: r1_extract.py bits|hist|dstar|merge|bsum [--nproc 2] [--sample data_428770] [--parts 0-19]"""
import argparse
import glob
import json
import os
import re
import sys
import time
import zlib
from multiprocessing import Pool

import numpy as np

E = "/eos/user/k/klyczak"
RAW = E + "/charm_raw_v20"
APP = E + "/charm_q/dnet4/app/d150"
WPSF = E + "/charm_q/wpscan/wps_hist.npz"
DN = E + "/charm_q/dnet/data"
SC = E + "/charm_q/dnet4/sc"
WD = "/tmp/klyczak_reg428"
CH = [("P", "HLT_noalg_L1AFP_A_OR_C"), ("P", "HLT_j20_L1AFP_A_OR_C"), ("P", "HLT_noalg_L1J12"), ("T", "L1_J12"), ("T", "L1_AFP_A_OR_C"),
      ("P", "HLT_j20_L1J12"), ("P", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"), ("P", "HLT_j20f_L1AFP_A_OR_C"), ("P", "HLT_noalg_L1AFP_A_OR_C_J12"),
      ("P", "HLT_j20_L1AFP_A_OR_C_J12"), ("P", "HLT_j20_L1RD0_FILLED"), ("P", "HLT_j20_L1MBTS_2")]          # bit = indeks + 1
PSC = ["HLT_noalg_L1AFP_A_OR_C", "HLT_j20_L1AFP_A_OR_C", "HLT_noalg_L1J12"]
PRE = {"P": "EventInfo_trigPassed_", "T": "EventInfo_trigTBP_", "S": "EventInfo_trigPS_"}
BINS = [(3, 5), (5, 8), (8, 30)]
FINE = [3, 4, 5, 6.5, 8, 10, 12, 16, 30]
EPS_DESC = [round(0.95 - 0.05 * i, 2) for i in range(18)]
J70 = EPS_DESC.index(0.70)                                           # WP70 ⇔ s ≥ J70 + 1
SP = {"D0": {"m0": 1.86484, "win": 0.055, "mr": (1.70, 2.05)}, "Dplus": {"m0": 1.86966, "win": 0.045, "mr": (1.70, 2.20)}}
TAG = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}
MBW = 0.001
NREP = 30
REPCLS = ["allGRL", "R1", "R2", "R3x", "rest", "J12", "j20AFP", "J12hi", "R3xhi", "R2hi"]
LBHI = 16                                                            # LB 4–16: preskala L1_AFP = 2


def log(s):
    print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)


def res(keys, name):
    for n in (name, name + "_NOSYS"):
        if n in keys:
            return n
    return None


# ------------------------------------------------------------------------------------------------ bits
def bits_one(f):
    import uproot
    fnum = int(re.search(r"\._(\d{6})\.", os.path.basename(f)).group(1))
    task = int(os.path.basename(f).split(".")[2])
    out = f"{WD}/bits/b_{fnum:06d}.npz"
    if os.path.isfile(out):
        return fnum, None
    tr = uproot.open(f)["CharmAnalysis"]
    keys = set(tr.keys())
    want = {"grl": res(keys, "EventInfo_passGRL_MinBiasAFP"), "lb": res(keys, "EventInfo_lumiBlock"), "ev": res(keys, "EventInfo_eventNumber"),
            "run": res(keys, "EventInfo_runNumber")}
    miss = []
    for i, (k, c) in enumerate(CH):
        b = res(keys, PRE[k] + c)
        if b:
            want[f"b{i + 1}"] = b
        else:
            miss.append(k + ":" + c)
    for c in PSC:
        b = res(keys, PRE["S"] + c)
        if b:
            want["ps_" + c] = b
    a = tr.arrays(list(want.values()), library="np")
    g = lambda n: a[want[n]]  # noqa: E731
    n = len(g("lb"))
    bits = (g("grl") == 1).astype(np.uint16)
    for i in range(len(CH)):
        if f"b{i + 1}" in want:
            bits |= (g(f"b{i + 1}") == 1).astype(np.uint16) << np.uint16(i + 1)
    lb = g("lb").astype(np.uint16)
    ps = {}
    for c in PSC:
        if "ps_" + c in want:
            v = g("ps_" + c).astype(np.float64)
            for l in np.unique(lb):
                ps.setdefault(c, {})[int(l)] = sorted(set(np.round(v[lb == l], 3).tolist()))[:4]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out + ".tmp.npz", bits=bits, lb=lb, ev=g("ev").astype(np.uint64), task=np.uint64(task), run=np.unique(g("run")),
             info=np.array(json.dumps({"n": int(n), "miss": miss, "ps": ps, "task": task})))
    os.replace(out + ".tmp.npz", out)
    return fnum, {"n": int(n), "miss": miss, "ps": ps, "task": task}


def run_bits(a):
    fs = sorted(glob.glob(f"{RAW}/{a.sample}/*.root"))
    log(f"bits: {len(fs)} plików ROOT")
    info = {}
    pf = f"{WD}/bits/ps_{a.sample}.json"
    if os.path.isfile(pf):
        info = json.load(open(pf))
    with Pool(a.nproc) as P:
        for i, (fnum, r) in enumerate(P.imap_unordered(bits_one, fs)):
            if r is not None:
                info[str(fnum)] = r
            if (i + 1) % 10 == 0:
                log(f"  {i + 1}/{len(fs)}")
    json.dump(info, open(pf, "w"))
    log("bits OK")


_B = {}


def load_bits(fnum):
    if fnum not in _B:
        if len(_B) > 6:
            _B.clear()
        z = np.load(f"{WD}/bits/b_{fnum:06d}.npz")
        _B[fnum] = (z["bits"], z["lb"], z["ev"], int(z["task"]))
    return _B[fnum]


def map_uid(uid):
    """uid → (bits, lb, ev); −1/0 gdy brak. Sprawdza task."""
    uid = uid.astype(np.uint64)
    fn = ((uid >> np.uint64(22)) & np.uint64(0x3FFF)).astype(np.int64)
    en = (uid & np.uint64(0x3FFFFF)).astype(np.int64)
    tk = (uid >> np.uint64(36)).astype(np.int64)
    bits = np.zeros(len(uid), np.uint16); lb = np.zeros(len(uid), np.uint16); ev = np.zeros(len(uid), np.uint64); ok = np.zeros(len(uid), bool)
    for f in np.unique(fn):
        m = fn == f
        try:
            B, L, EV, task = load_bits(int(f))
        except FileNotFoundError:
            continue
        mm = m & (tk == task) & (en < len(B))
        bits[mm], lb[mm], ev[mm], ok[mm] = B[en[mm]], L[en[mm]], EV[en[mm]], True
    return bits, lb, ev, ok


def bit(b, i):
    return (b & np.uint16(1 << i)) > 0


def classes428(b, lb):
    grl, r1, j20, j12 = bit(b, 0), bit(b, 1), bit(b, 2), bit(b, 3)
    hi = lb > LBHI
    R1, R2 = grl & r1, grl & j20 & ~r1
    R3x = grl & j12 & ~r1 & ~j20
    return {"allGRL": grl, "R1": R1, "R2": R2, "R3x": R3x, "rest": grl & ~r1 & ~j20 & ~j12, "J12": grl & j12, "j20AFP": grl & j20,
            "J12hi": grl & j12 & hi, "R3xhi": R3x & hi, "R2hi": R2 & hi,
            "R1_J12": R1 & j12, "R1_j20": R1 & j20, "R1hi": R1 & hi, "R1hi_j20": R1 & j20 & hi, "R1hi_J12": R1 & j12 & hi, "R1_J12tbp": R1 & bit(b, 4),
            "R1_j20J12": R1 & bit(b, 6), "R1_j20f": R1 & bit(b, 8),
            "J12_AFPtbp": grl & j12 & bit(b, 5), "J12_j20AFP": grl & j12 & j20, "J12_R1": grl & j12 & r1, "J12_j20J12": grl & j12 & bit(b, 6),
            "j20AFP_J12": grl & j20 & j12, "allGRL_AFPtbp": grl & bit(b, 5), "nonGRL": ~grl}


def classes427(cls):
    grl, mb1, j12 = (cls & 1) > 0, (cls & 2) > 0, (cls & 4) > 0
    return {"MB": mb1 & grl, "MB_J12": mb1 & grl & j12}


# ------------------------------------------------------------------------------------------------ hist
def hist_part(t):
    sample, part = t
    out = f"{WD}/hist/h_{sample}_{part}.npz"
    if os.path.isfile(out):
        return out, "jest"
    W = np.load(WPSF, allow_pickle=False)
    fs = sorted(glob.glob(f"{APP}/app_{sample}_{part}_*.npz"))
    is428 = sample == "data_428770"
    acc, val = {}, {"rows": 0, "nomap": 0, "grl_mismatch": 0, "j12_mismatch": 0, "files": len(fs)}
    for f in fs:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        for sp, c in SP.items():
            if f"{sp}__m" not in z.files:
                continue
            jt = meta["tags"][sp].index(TAG[sp])
            lo, hi = c["mr"]
            nm = int(round((hi - lo) / MBW))
            m, pt, ok, cls, uid = z[f"{sp}__m"], z[f"{sp}__pt"], z[f"{sp}__ok"], z[f"{sp}__cls"], z[f"{sp}__uid"]
            lg = z[f"{sp}__lg"][:, jt].astype(np.float64)
            ib = np.floor((m.astype(np.float64) - lo) / MBW + 1e-9).astype(np.int64)
            base = ok & (ib >= 0) & (ib < nm)
            val["rows"] += int(len(m))
            if is428:
                b, lb, _, mok = map_uid(uid)
                val["nomap"] += int((~mok).sum())
                val["grl_mismatch"] += int((mok & (bit(b, 0) != ((cls & 1) > 0))).sum())
                val["j12_mismatch"] += int((mok & (bit(b, 3) != ((cls & 4) > 0))).sum())
                base &= mok
                cl = classes428(b, lb)
                rng = np.random.default_rng(zlib.crc32(os.path.basename(f).encode()) + 1000003 * (1 if sp == "D0" else 2))
                u, inv = np.unique(uid, return_inverse=True)
                wr = rng.poisson(1.0, (len(u), NREP)).astype(np.uint8)
            else:
                cl = classes427(cls)
            nS = len(EPS_DESC) + 1
            ifn = np.digitize(pt, FINE) - 1
            sfine = np.zeros(len(m), np.int64)
            for b_ in BINS:
                bn = f"{b_[0]}-{b_[1]}"
                inb = (pt >= b_[0]) & (pt < b_[1])
                thr = W[f"{sp}|{TAG[sp]}|{bn}|thr"]
                s = np.searchsorted(thr, lg, side="right")
                sfine[inb] = s[inb]
                sel = base & inb
                for cn, cm in cl.items():
                    mm = sel & cm
                    key = f"{sp}|{bn}|{sample[5:]}|{cn}"
                    h = acc.setdefault(key, np.zeros(nS * nm, np.uint32))
                    if mm.any():
                        h += np.bincount(s[mm] * nm + ib[mm], minlength=nS * nm).astype(np.uint32)
                    if is428 and cn in REPCLS:
                        hr = acc.setdefault(key + "|rep70", np.zeros((NREP, nm), np.uint32))
                        mp = mm & (s >= J70 + 1)
                        if mp.any():
                            wm = wr[inv[mp]]
                            for r in range(NREP):
                                hr[r] += np.bincount(ib[mp], weights=wm[:, r], minlength=nm).astype(np.uint32)
            nf = len(FINE) - 1
            self_ = base & (ifn >= 0) & (ifn < nf)
            pas = (sfine >= J70 + 1).astype(np.int64)
            for cn, cm in cl.items():
                mm = self_ & cm
                key = f"{sp}|fine|{sample[5:]}|{cn}"
                h = acc.setdefault(key, np.zeros(nf * 2 * nm, np.uint32))
                if mm.any():
                    h += np.bincount((ifn[mm] * 2 + pas[mm]) * nm + ib[mm], minlength=nf * 2 * nm).astype(np.uint32)
    o = {}
    for k, h in acc.items():
        sp = k.split("|")[0]
        nm = int(round((SP[sp]["mr"][1] - SP[sp]["mr"][0]) / MBW))
        o[k] = h if h.ndim == 2 else h.reshape((-1, nm)) if "|fine|" not in k else h.reshape((len(FINE) - 1, 2, nm))
    o["val"] = np.array(json.dumps(val))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez_compressed(out + ".tmp.npz", **o)
    os.replace(out + ".tmp.npz", out)
    return out, val


def run_hist(a):
    lo, hi = (int(x) for x in a.parts.split("-"))
    tasks = [(a.sample, p) for p in range(lo, hi + 1)]
    with Pool(a.nproc) as P:
        for out, v in P.imap_unordered(hist_part, tasks):
            log(f"{out}: {v}")


def run_merge(a):
    acc, vals = {}, {}
    for f in sorted(glob.glob(f"{WD}/hist/h_data_*_*.npz")):
        if ".tmp" in f:
            continue
        z = np.load(f)
        for k in z.files:
            if k == "val":
                vals[os.path.basename(f)] = json.loads(str(z[k])); continue
            acc[k] = acc[k] + z[k] if k in acc else z[k].copy()
    W = np.load(WPSF, allow_pickle=False)
    for sp in SP:
        for b_ in BINS:
            k = f"{sp}|{TAG[sp]}|{b_[0]}-{b_[1]}|thr"
            acc["thr|" + k] = W[k]
    acc["meta"] = np.array(json.dumps({"val": vals, "eps_desc": EPS_DESC, "fine": FINE, "bins": BINS, "sp": SP, "tag": TAG, "mbw": MBW, "nrep": NREP,
                                       "repcls": REPCLS, "lbhi": LBHI, "bits": [k + ":" + c for k, c in CH]}))
    np.savez_compressed(f"{WD}/reg_hist.npz", **acc)
    tot = {k: sum(v[k] for v in vals.values()) for k in ("rows", "nomap", "grl_mismatch", "j12_mismatch", "files")}
    log(f"merge: {len(vals)} części, {tot}, {os.path.getsize(WD + '/reg_hist.npz') / 1e6:.1f} MB")


# ------------------------------------------------------------------------------------------------ bsum
def run_bsum(a):
    names = ["GRL"] + [c for _, c in CH]
    tot, perlb, n = {}, {}, 0
    combos = {}
    for f in sorted(glob.glob(f"{WD}/bits/b_*.npz")):
        z = np.load(f)
        b, lb = z["bits"], z["lb"]
        n += len(b)
        g = bit(b, 0)
        for i, nmn in enumerate(names):
            tot[nmn] = tot.get(nmn, 0) + int((g & bit(b, i)).sum())
        cl = classes428(b, lb)
        for cn, cm in cl.items():
            combos[cn] = combos.get(cn, 0) + int(cm.sum())
        for cn in ("allGRL", "R1", "j20AFP", "J12", "R2", "R3x"):
            u, c = np.unique(lb[cl[cn]], return_counts=True)
            d = perlb.setdefault(cn, {})
            for l, k in zip(u.tolist(), c.tolist()):
                d[l] = d.get(l, 0) + k
    ps = {}
    for f in sorted(glob.glob(f"{WD}/bits/b_*.npz")):
        z = np.load(f)
        if "info" in z.files:
            ps[f] = json.loads(str(z["info"]))
    psm = {}
    for r in ps.values():
        for c, d in r["ps"].items():
            for l, v in d.items():
                psm.setdefault(c, {}).setdefault(l, set()).update(v)
    psm = {c: {l: sorted(v) for l, v in sorted(d.items(), key=lambda x: int(x[0]))} for c, d in psm.items()}
    miss = sorted({m for r in ps.values() for m in r["miss"]})
    json.dump({"n_events": n, "grl_and_bit": tot, "classes_events": combos, "perlb": perlb, "ps": psm, "missing_branches": miss, "files": len(ps)},
              open(f"{WD}/bits_summary.json", "w"))
    log(f"bsum: zdarzeń {n}, {combos}, brak gałęzi {miss}")


# ------------------------------------------------------------------------------------------------ dstar
def run_dstar(a):
    import zipfile
    o = {}
    m0, win = SP["D0"]["m0"], SP["D0"]["win"]
    for s in ("data_428770", "data_427929", "fd_nd"):
        z = np.load(f"{DN}/dn_{s}.npz", allow_pickle=False)
        sc = np.load(f"{SC}/{TAG['D0']}.npz")[f"{s}__p"]
        g = lambda c: z["D0__p_" + c]  # noqa: E731
        m, pt = g("m"), g("pt")
        sel = (np.abs(m - m0) < win) & (pt >= 3) & (pt < 30) & np.isfinite(sc)
        if s == "fd_nd":
            sel &= np.isin(g("truth"), (0, 1)) & np.isin(g("truth_ds"), (0, 1))
        i = np.flatnonzero(sel)
        log(f"dstar {s}: wierszy p_ {len(m)}, w oknie {len(i)}")
        o[f"{s}|dm"], o[f"{s}|pt"], o[f"{s}|sc"], o[f"{s}|m"] = g("dm")[i].astype(np.float32), pt[i].astype(np.float32), sc[i].astype(np.float32), m[i].astype(np.float32)
        if s == "fd_nd":
            continue
        uid = g("uid")[i]
        o[f"{s}|uid"] = uid
        tp = lambda c: (g("e_tp_" + c)[i] == 1)  # noqa: E731
        grl = g("e_grl")[i] == 1
        o[f"{s}|cls"] = (grl * 1 + tp("HLT_mb_mbts_L1MBTS_1") * 2 + tp("HLT_noalg_L1J12") * 4 + tp("HLT_noalg_L1AFP_A_OR_C") * 8).astype(np.uint8)
        if s == "data_428770":
            b, lb, ev, mok = map_uid(uid)
            o[f"{s}|bits"], o[f"{s}|lb"], o[f"{s}|mapok"] = b, lb, mok
            evt = g("event")[i].astype(np.uint64)
            chk = {"n": int(len(i)), "nomap": int((~mok).sum()), "event_mismatch": int((mok & (ev != evt)).sum()),
                   "grl_mismatch": int((mok & (bit(b, 0) != grl)).sum()), "j12_mismatch": int((mok & (bit(b, 3) != tp("HLT_noalg_L1J12"))).sum()),
                   "afp_mismatch": int((mok & (bit(b, 1) != tp("HLT_noalg_L1AFP_A_OR_C"))).sum()), "lb_mismatch": int((mok & (lb != g("e_lb")[i])).sum())}
            log(f"dstar walidacja mapy uid: {chk}")
            o["check"] = np.array(json.dumps(chk))
    np.savez_compressed(f"{WD}/dstar.npz", **o)
    log(f"zapisano dstar.npz {os.path.getsize(WD + '/dstar.npz') / 1e6:.1f} MB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["bits", "hist", "dstar", "merge", "bsum"])
    ap.add_argument("--nproc", type=int, default=2); ap.add_argument("--sample", default="data_428770"); ap.add_argument("--parts", default="0-19")
    a = ap.parse_args()
    os.makedirs(WD, exist_ok=True)
    {"bits": run_bits, "hist": run_hist, "dstar": run_dstar, "merge": run_merge, "bsum": run_bsum}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
