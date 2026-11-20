"""TAGREG — wspólne: wczytanie zdarzeń (AFP_mix/runs/lxplus/*.parquet, tylko odczyt) do cache npz w TAGREG/cache.
Żadnych fitów; liczby kandydatów tylko w całym oknie slim albo w pasmach bocznych."""
import os, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.dirname(HERE)
IN = os.path.join(W, "AFP_mix", "runs", "lxplus")
BT = os.path.join(W, "B_trig", "runs", "lxplus")
CACHE = os.path.join(HERE, "cache"); OUT = os.path.join(HERE, "out")
RUNS = (427929, 428770)
SIDES = ("A", "C")
DELTA = 0.02
MBTS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
J12 = "HLT_noalg_L1J12"; AN = "HLT_noalg_L1AFP_A_OR_C"; AS = "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"; RD0 = "HLT_noalg_L1RD0_FILLED"
CH = MBTS + [J12, AN, AS, RD0]
L1 = ["L1_MBTS_1", "L1_MBTS_2", "L1_J12", "L1_AFP_A_OR_C", "L1_RD0_FILLED"]
SPECIES = {"Dplus": dict(code=0, var="Dmeson_mass", sb=((1.75, 1.80), (1.94, 2.00))),
           "D0": dict(code=3, var="Dmeson_mass", sb=((1.75, 1.80), (1.93, 2.00))),
           "Dstar": dict(code=1, var="Dmeson_delta_mass", sb=((0.150, 0.160),))}
PT = ((2.0, 3.0), (3.0, 5.0), (5.0, 8.0), (8.0, 30.0))
SBC = [f"{sp}|{lo:g}-{hi:g}" for sp in SPECIES for lo, hi in PT]   # kolumny liczby kandydatów w pasmach bocznych
NSEL_EDGES = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 80, 10000])
DISJ = ["AFPnoalg", "AFPsptrk4", "AFPboth", "AFP&J12", "AFP&MBTS", "J12only!A", "MBTS!A", "none"]


def load(run, force=False):
    os.makedirs(CACHE, exist_ok=True)
    fn = os.path.join(CACHE, f"ev_{run}.npz")
    if os.path.exists(fn) and not force:
        z = np.load(fn, allow_pickle=False); return {k: z[k] for k in z.files}
    ev = pd.read_parquet(os.path.join(IN, f"data_{run}_events.parquet"))
    n = len(ev); evno = ev["event_number"].to_numpy()
    assert len(np.unique(evno)) == n, "event_number nieunikalny"
    ca = pd.read_parquet(os.path.join(IN, f"data_{run}_cands.parquet"), columns=["event_number", "Dmeson_species", "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_pt", "Dmeson_Lxy"])
    order = np.argsort(evno); pos = np.searchsorted(evno[order], ca["event_number"].to_numpy()); crow = order[pos]
    assert (evno[crow] == ca["event_number"].to_numpy()).all()
    spec = ca["Dmeson_species"].to_numpy(); pt = ca["Dmeson_pt"].to_numpy()
    out = {"evno": evno.astype(np.int64), "lb": ev["lumiBlock"].to_numpy().astype(np.int32), "mu": ev["actual_mu"].to_numpy().astype(np.float32),
           "grl": ev["passGRL_MinBiasAFP"].to_numpy() == 1, "npu": ev["nPileupVertices"].to_numpy().astype(np.int16),
           "nsel": ev["nSelBL"].to_numpy().astype(np.int32), "n_cand": ev["n_cand"].to_numpy(), "n_cand_sel": ev["n_cand_sel"].to_numpy()}
    for s in SIDES:
        out[f"xical_{s}"] = ev[f"xi_cal_{s}"].to_numpy().astype(np.float64)
        out[f"np_{s}"] = ev[f"afp_nProtons_{s}"].to_numpy(); out[f"nw_{s}"] = ev[f"afp_nInWindow_{s}"].to_numpy()
        out[f"xi1_{s}"] = ev[f"afp_xi1_{s}"].to_numpy().astype(np.float64); out[f"tag_{s}"] = ev[f"afp_tag_{s}"].to_numpy() == 1
    for c in CH:
        out["T_" + c] = ev["trigPassed_" + c].to_numpy() == 1
    for c in L1:
        out["TBP_" + c] = ev["trigTBP_" + c].to_numpy() == 1
    # kandydaci: dowolny (slim), D0 pT>=3 (nasza luźna selekcja), pasma boczne per gatunek × pT
    hasD = np.zeros(n, bool); hasD[crow[np.isin(spec, [0, 1, 3])]] = True
    out["hasD"] = hasD
    d0 = (spec == 3) & (pt >= 3.0)
    h = np.zeros(n, bool); h[crow[d0]] = True; out["hasD0"] = h
    out["d0_lxy_min"] = np.array(float(ca["Dmeson_Lxy"].to_numpy()[spec == 3].min()))
    out["pt_min"] = np.array([float(pt[spec == d["code"]].min()) for d in SPECIES.values()])
    sbc = np.zeros((n, len(SBC)), np.uint16); anysb = np.zeros(n, bool); ci = 0
    for sp, d in SPECIES.items():
        ms = spec == d["code"]; v = ca[d["var"]].to_numpy(); sb = np.zeros(len(ca), bool)
        for lo, hi in d["sb"]:
            sb |= ms & (v >= lo) & (v <= hi)
        for lo, hi in PT:
            m = sb & (pt >= lo) & (pt < hi)
            sbc[:, ci] = np.minimum(np.bincount(crow[m], minlength=n), 65535); ci += 1
        hs = np.zeros(n, bool); hs[crow[sb & (pt >= 3.0)]] = True; out[f"hasSB_{sp}"] = hs; anysb |= hs
    out["sbc"] = sbc; out["hasSB"] = anysb
    np.savez(fn, **out)
    return out


def classes(d):
    M = np.zeros(len(d["lb"]), bool)
    for c in MBTS:
        M |= d["T_" + c]
    J, an, as_ = d["T_" + J12], d["T_" + AN], d["T_" + AS]; A = an | as_
    cl = {"AFPnoalg": an & ~as_ & ~J & ~M, "AFPsptrk4": as_ & ~an & ~J & ~M, "AFPboth": an & as_ & ~J & ~M,
          "AFP&J12": A & J & ~M, "AFP&MBTS": A & M, "J12only!A": J & ~M & ~A, "MBTS!A": M & ~A, "none": ~J & ~M & ~A}
    return cl, dict(A=A, J=J, M=M, an=an, as_=as_)


def regions(d):
    """T0 = GRL ∧ tag(A∨C); T1 = T0 ∧ npu==0; T2 = T1 ∧ ξ_cal-pass na każdej tagowanej stronie."""
    tag = d["tag_A"] | d["tag_C"]
    ok = np.ones(len(tag), bool)
    for s in SIDES:
        ok &= ~d[f"tag_{s}"] | (d[f"xical_{s}"] <= d[f"xi1_{s}"] + DELTA)
    T0 = d["grl"] & tag; T1 = T0 & (d["npu"] == 0); T2 = T1 & ok
    return T0, T1, T2


def trigps(run):
    m = json.load(open(os.path.join(IN, f"data_{run}.meta.json")))
    return m["trigps_columns"], np.array(m["trigps"])


def jdump(o, path):
    def conv(x):
        if isinstance(x, (np.floating, np.integer)): return x.item()
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.bool_): return bool(x)
        raise TypeError(type(x))
    json.dump(o, open(path, "w"), indent=1, default=conv, ensure_ascii=False)
