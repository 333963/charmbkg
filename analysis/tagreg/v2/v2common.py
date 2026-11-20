"""TAGREG v2 — wspólne: wczytanie bits_<run>.parquet, rozpakowanie bitów, regiony, rodziny. Tylko liczniki."""
import json, os, re
import numpy as np, pandas as pd, pyarrow.parquet as pq
D = os.path.expanduser("~/Documents/charm/charm_work/TAGREG_2026-11-18/runs/lxplus")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RUNS = (428770, 427929)
TRIG11 = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
          "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2"]
FAM = ["AFPnoalg", "AFPnoalg-inne", "AFP+tor", "AFP+dzet", "dzet", "MBTS", "pozostale", "noneHLT"]

def family(c):
    if c == "HLT_noalg_L1AFP_A_OR_C": return 0
    if c in ("HLT_noalg_L1AFP_A_OR_C_J12", "HLT_noalg_L1MU5VF_AFP_A_OR_C"): return 3
    if c.startswith("HLT_noalg_L1AFP_"): return 1
    if re.match(r"HLT_mb_sptrk_pt\d_L1AFP_", c): return 2
    if re.match(r"HLT_(j20f?|2j10_20detaAA)_L1AFP_", c): return 3
    if re.match(r"HLT_(j\d+f?|2j\d+(_20detaAA)?)_L1", c) or re.match(r"HLT_noalg_L1J\d+$", c): return 4
    if "L1MBTS" in c: return 5
    return 6

def unpack(W, k):
    M = np.zeros((W.shape[0], k), bool)
    for j in range(k):
        M[:, j] = (W[:, j // 64] >> np.uint64(j % 64)) & np.uint64(1)
    return M

def load(run):
    ch = json.load(open(f"{D}/chains_{run}.json"))
    tb = pq.read_table(f"{D}/bits_{run}.parquet").to_pandas()
    names, l1 = ch["chains"], ch["l1tbp"]
    P = unpack(np.stack([tb[c].to_numpy() for c in ("trig_w0", "trig_w1")], 1), len(names))
    T = unpack(tb[["tbp_w0"]].to_numpy(), len(l1))
    ih = [j for j, c in enumerate(names) if c.startswith("HLT_")]
    il = [j for j, c in enumerate(names) if c.startswith("L1_")]
    hlt = [names[j] for j in ih]
    l1av = [names[j] for j in il]
    assert l1av == l1, "L1 po preskali vs TBP: inne listy"
    ev = tb.drop(columns=[c for c in tb.columns if c.startswith(("trig_w", "tbp_w"))])
    a, c = ev.tag_A.to_numpy() > 0, ev.tag_C.to_numpy() > 0
    ok = (~a | (ev.xi_cal_A.to_numpy() <= ev.xi_A.to_numpy() + 0.02)) & (~c | (ev.xi_cal_C.to_numpy() <= ev.xi_C.to_numpy() + 0.02))
    T1 = ev.npu.to_numpy() == 0
    reg = {"T0": np.ones(len(ev), bool), "T1": T1, "T2": T1 & ok}
    cand = {"D0": ev.n_D0.to_numpy() > 0, "Dp": ev.n_Dp.to_numpy() > 0, "Dst": ev.n_Dst.to_numpy() > 0}
    cand["D"] = cand["D0"] | cand["Dp"] | cand["Dst"]
    return dict(ev=ev, H=P[:, ih], hlt=hlt, AV=P[:, il], TBP=T, l1=l1, reg=reg, cand=cand, a=a, c=c,
                fam=np.array([family(x) for x in hlt]))

def famclass(H, fam):
    """Rozłączna rodzina zdarzenia: pierwsza (najniższy indeks) rodzina z bitem; 7 = noneHLT."""
    out = np.full(H.shape[0], 7)
    for f in range(6, -1, -1):
        out[H[:, fam == f].any(1)] = f
    return out

def lbblocks(lb, T1, n=8):
    q = np.quantile(lb[T1], np.linspace(0, 1, n + 1)); q[-1] += 1
    e = np.unique(np.round(q).astype(int))
    return e
