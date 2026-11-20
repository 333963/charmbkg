#!/usr/bin/env python3
"""MCTMPL — histogramy do testu szablonów tła z MC (zadanie #7 planu prac; projekt B, 17.11.2026).

Wejście: h5 v2.0 (compound `events`), wyjście: JEDEN mały npz z histogramami (liczniki int64) — żadnych tabel kandydatów.
Selekcja kandydatów (nasza, nie CharmNet): D0 (decayType 10) L_xy ≥ --lxy-d0 (0,25 mm), D+ (decayType 4) L_xy ≥ --lxy-dp
(ntupla: 0,5 mm), 3 ≤ pT < 30 GeV, |η| < 2,2, okna masy jak w ficie TRExFitter głównego łańcucha: D0 1,70–2,05, D+ 1,70–2,20 GeV.
Score = nasz klasyfikator H (modele D0_H / Dplus_H, OOF po foldach); oś `sb` = liczba przekroczonych progów
(0: < WP95, 1: [WP95, WP80), 2: [WP80, WP50), 3: ≥ WP50) — progi z wp_thresholds.json (tryb --thr-out je wyznacza z fd_nd).
Dane: wymagany passGRL; ZDARZENIA REGIONU SYGNAŁU GŁÓWNEGO ŁAŃCUCHA SĄ ODRZUCANE I NIGDZIE NIE LICZONE
  (tag AFP po którejkolwiek stronie ∧ nPileupVertices == 0 ∧ ξ_cal ≤ ξ_AFP + 0,02 po każdej tagowanej stronie).
Oś K: dane — klasy trygera (nierozłączne, lista w meta); MC — kategoria prawdy `fc` (24 wartości, FC_NAMES; grupy głównego łańcucha (TRExFitter) GROUPS).
Oś `cell` (49) = bin nSelBL × bin pT toru wiodącego — do przeważenia MC do aktywności klasy trygera w analizie.
Histogramy: patrz SPEC (nazwa → osie). Użycie:
  mctmpl_hist.py --models DIR(models) --wp wp_thresholds.json --tag data|mb|fd --run 428770 --out OUT.npz FILE.h5 [...]
  mctmpl_hist.py --models DIR --tag fd --thr-out wp_thresholds.json FILE.h5 [...]     (progi WP z sygnału prawdy)
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import yaml

VERSION = "1.0"
PT_EDGES = np.array([3.0, 5.0, 8.0, 30.0])
PT_NAMES = ["3-5", "5-8", "8-30"]
WPS = ["0.95", "0.80", "0.50"]
ETA_MAX = 2.2
M0 = {"D0": 1.86484, "Dplus": 1.86966}
MWIN = {"D0": (1.70, 2.05, 70), "Dplus": (1.70, 2.20, 100)}          # 5 MeV; fit TRExFitter głównego łańcucha: 10 MeV (rebin ×2)
SIGWIN = {"D0": 0.055, "Dplus": 0.045}                               # okno do progów WP (jak config charmtag)
MREG = (0.040, 0.075, 0.135)                                         # SR |Δm| < 40 MeV; SB 75–135 MeV (lo, hi)
NSEL_EDGES = np.array([0, 20, 30, 40, 50, 60, 75, 1e9])
LEAD_EDGES = np.array([0, 3, 4, 5, 7, 10, 15, 1e9])
NCELL = (len(NSEL_EDGES) - 1) * (len(LEAD_EDGES) - 1)
MK, MPI, MP, MPHI, DM0 = 0.493677, 0.13957039, 0.93827209, 1.019461, 0.1454258
MBTS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
J12, AN, AS = "HLT_noalg_L1J12", "HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"
TRIG = MBTS + [J12, AN, AS]
CLASSES = {428770: ["all", "J12only", "MBTSbased", "AFPstream", "AFPonly", "none"], 427929: ["all", "MB"]}
FC_NAMES = ["SigC", "SigB", "ReflSwap", "ReflOther", "CM_411", "CM_421", "CM_431", "CM_413", "CM_423", "CM_4122", "CM_other",
            "DMM_C", "DMM_B", "PC_411", "PC_421", "PC_423", "PC_413", "PC_431", "PC_4122", "PC_other",
            "MatchedBL", "OtherNoB", "OtherB", "Unknown"]
# grupy = próbki TRExFitter głównego łańcucha (HistogramManager.cxx:555-556; Reflection = MatchedC/B z zamienionymi K/π, tylko D0)
GROUPS = {"Sig": [0, 1], "Reflection": [2, 3], "CharmMismatch": [4, 5, 6, 7, 8, 9, 10], "DecayModeMismatch": [11, 12],
          "PartComb": [13, 14, 15, 16, 17, 18, 19], "Comb": [20, 21, 22, 23]}
DS_NAMES = ["Dstar_ok", "Dstar_swapped", "other", "unused"]          # oś K histogramów ds_* w MC
EV = {"uid": "event_uid", "grl": "passGRL_MinBiasAFP", "nsel": "nSelBL", "lead": "leadTrackPtBL", "npu": "nPileupVertices",
      "tagA": "afp_tag_A", "tagC": "afp_tag_C", "xiA": "afp_xi1_A", "xiC": "afp_xi1_C", "xcA": "xi_cal_A", "xcC": "xi_cal_C"}
RAWX = {"truth": "Dmeson_truth_category", "hs": "Dmeson_hypo_swap", "Dtpdg": "Dmeson_truthPdgId", "dmass": "Dmeson_delta_mass",
        "chan": "Dmeson_truth_channel"}
for _i in range(3):
    RAWX.update({f"pdg{_i}": f"Dmeson_daughter{_i}_pdgId", f"tpdg{_i}": f"Dmeson_daughter{_i}_truthPdgId",
                 f"trk{_i}": f"Dmeson_daughter{_i}_trackId"})
NB = {"imp": 40, "ct": 40, "dsdm": 80, "m10": 35, "dpdm": 100, "kkpi": 80, "mkk": 60, "pkpi": 60}
# nazwa → (gatunek, osie po [pt, sb, K])
SPEC = {
    "d0_m": ("cell", "m"), "d0_imp": ("mreg", "imp"), "d0_ct": ("mreg", "ct"),
    "ds_right": ("dsdm", "m10"), "ds_twin": ("dsdm", "m10"),
    "dp_m": ("cell", "veto", "m"), "dp_imp": ("mreg", "imp"), "dp_ct": ("mreg", "ct"), "dp_dm": ("dpdm",),
    "dp_mkk": ("mkk",), "dp_kkpi": ("phi", "kkpi"), "dp_mphi": ("phi", "m"), "dp_pkpi": ("pkpi",), "dp_mds": ("dswin", "m"),
}
AXES = {  # opis osi (lo, hi, n) — zapisywany w meta
    "d0_m": MWIN["D0"], "dp_m": MWIN["Dplus"], "imp": (0.0, 0.2, 40), "ct_D0": (0.0, 1.0, 40), "ct_Dplus": (0.0, 4.0, 40),
    "dsdm": (0.1396, 0.1596, 80), "m10": (1.70, 2.05, 35), "dpdm": (0.135, 0.185, 100), "mkk": (0.98, 1.10, 60),
    "kkpi": (1.75, 2.15, 80), "pkpi": (2.15, 2.45, 60),
    "phi": "0: |m(KK) − 1,0195| < 10 MeV; 1: 1,035 < m(KK) < 1,055", "dswin": "0: 143,5 < Δm_min < 147,5 MeV; 1: 150–158 MeV",
    "veto": "0: Δm_min < 160 MeV (odrzucone wetem D*); 1: przechodzi weto", "mreg": "0: SR |m−m0| < 40 MeV; 1: SB dolne; 2: SB górne (75–135 MeV)",
}


def dims(name, K):
    sp = "D0" if name.startswith(("d0_", "ds_")) else "Dplus"
    d = {"cell": NCELL, "m": MWIN[sp][2], "mreg": 3, "veto": 2, "phi": 2, "dswin": 2}
    return (3, 4, K) + tuple(d[a] if a in d else NB[a] for a in SPEC[name])


def bidx(x, lo, hi, n, clip=False):
    """Indeks binu i maska ważności; clip=True → przepełnienie do skrajnych binów."""
    x = np.asarray(x, float)
    i = np.floor((x - lo) / (hi - lo) * n).astype(np.int64)
    ok = np.isfinite(x) & ((i >= 0) & (i < n) if not clip else True)
    return np.clip(i, 0, n - 1), ok


def minv(parts):
    """Masa niezmiennicza: parts = [(pt, eta, phi, m), ...] (GeV)."""
    E = px = py = pz = 0.0
    for pt, eta, phi, m in parts:
        x, y, z = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        E = E + np.sqrt(x * x + y * y + z * z + m * m)
        px, py, pz = px + x, py + y, pz + z
    return np.sqrt(np.maximum(E * E - px * px - py * py - pz * pz, 0.0))


def fine_cat(sp, t, hs, pdg, tpdg, Dtpdg, chan):
    """Kategoria prawdy fc (indeks FC_NAMES). Refleksja wg definicji w charmpp (DSelection.cxx:584-594):
    MatchedC/B i |pdgId hipotezy| ≠ |truthPdgId| dla którejś z dwóch pierwszych córek (tylko D0)."""
    n = len(t)
    fc = np.full(n, 23, np.int64)
    sig = (t == 0) | (t == 1)
    fc[sig] = t[sig]
    if sp == "D0":
        sw = sig & ((np.abs(pdg[0]) != np.abs(tpdg[0])) | (np.abs(pdg[1]) != np.abs(tpdg[1])))
        fc[sw] = np.where(hs[sw] == 1, 2, 3)
    a = np.abs(Dtpdg)
    cm = (t == 2) | (t == 3)
    fc[cm] = 10
    for code, pid in ((4, 411), (5, 421), (6, 431), (7, 413), (8, 423), (9, 4122)):
        fc[cm & (a == pid)] = code
    fc[t == 6], fc[t == 7] = 11, 12
    pc = (t >= 8) & (t <= 11)
    fc[pc] = 19
    if pc.any():
        par = pd.Series(chan[pc]).astype(str).str.extract(r"(\d+)'?$")[0].fillna("0").astype(np.int64).to_numpy()
        sub = np.full(par.size, 19, np.int64)
        for code, pid in ((13, 411), (14, 421), (15, 423), (16, 413), (17, 431), (18, 4122)):
            sub[par == pid] = code
        fc[pc] = sub
    fc[(t == 4) | (t == 5)] = 20
    fc[t == 12] = 21
    fc[(t >= 13) & (t <= 15)] = 22
    return fc


class Hists:
    def __init__(self, K, Kds):
        self.K, self.Kds = K, Kds
        self.h = {n: np.zeros(dims(n, Kds if n.startswith("ds_") else K), np.int64) for n in SPEC}

    def fill(self, name, pt, sb, keys, rest, sel=None):
        """keys: tablica int (MC: fc) albo lista masek bool (dane: klasy). rest: lista (indeks, ważny)."""
        shape = self.h[name].shape
        ok = np.ones(len(pt), bool) if sel is None else sel.copy()
        idx = []
        for i, v in rest:
            idx.append(i)
            ok &= v
        ks = [(k, m) for k, m in enumerate(keys)] if isinstance(keys, list) else [(None, None)]
        for k, m in ks:
            o = ok & m if k is not None else ok
            if not o.any():
                continue
            kk = np.full(int(o.sum()), k, np.int64) if k is not None else keys[o]
            flat = np.ravel_multi_index((pt[o], sb[o], kk) + tuple(i[o] for i in idx), shape)
            self.h[name] += np.bincount(flat, minlength=self.h[name].size).reshape(shape)


def read_file(path, cols, chunk):
    import h5py
    with h5py.File(path, "r") as f:
        ds = f["events"]
        n = ds.shape[0]
        for a in range(0, n, chunk):
            x = ds.fields(cols)[a:min(n, a + chunk)]
            yield pd.DataFrame({c: x[c] for c in cols})


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True, help="katalog z D0_H i Dplus_H")
    ap.add_argument("--wp", default="")
    ap.add_argument("--thr-out", default="", help="tryb progów: zapisz kwantyle score sygnału prawdy (MC) i zakończ")
    ap.add_argument("--out", default="")
    ap.add_argument("--tag", required=True, help="data|mb|fd")
    ap.add_argument("--run", type=int, default=0, help="dane: 428770|427929 (klasy trygera)")
    ap.add_argument("--chunk", type=int, default=1_500_000)
    ap.add_argument("--lxy-d0", type=float, default=0.25)
    ap.add_argument("--lxy-dp", type=float, default=0.50)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    is_mc = a.tag != "data"
    if not is_mc and a.run not in CLASSES:
        sys.exit("BŁĄD: dane wymagają --run 428770|427929")
    if not a.thr_out and not (a.wp and a.out):
        sys.exit("BŁĄD: podaj --wp i --out (albo --thr-out)")

    from charmtag import features as F
    from charmtag.io import resolve_columns
    from charmtag.model import Ensemble

    MD = {"D0": os.path.join(a.models, "D0_H"), "Dplus": os.path.join(a.models, "Dplus_H")}
    cfg = yaml.safe_load(open(os.path.join(MD["D0"], "config_used.yaml")))
    units = cfg.get("units") or {}
    upt = float(units.get("d_pt", 1.0))
    ens = {s: Ensemble(MD[s]) for s in MD}
    spc = {s: cfg["species"][s] for s in MD}
    thr = json.load(open(a.wp)) if a.wp else None
    cls_names = FC_NAMES if is_mc else CLASSES[a.run]
    H = Hists(len(cls_names), 4 if is_mc else len(cls_names))
    thr_scores = {s: [[] for _ in PT_NAMES] for s in MD}
    cnt = {"rows": 0, "D0_sel": 0, "Dplus_sel": 0, "Dstar_rows": 0, "ds_right": 0, "ds_twin": 0}
    t0 = time.time()

    for path in a.files:
        import h5py
        with h5py.File(path, "r") as f:
            names = list(f["events"].dtype.names)
        cm = resolve_columns(names, cfg["columns"], 3)
        dt_col = cm.get("decay_type")
        rawx = {k: c for k, c in RAWX.items() if is_mc or k not in ("chan",)}
        miss = [c for c in list(rawx.values()) + list(EV.values()) if c not in names]
        if miss:
            sys.exit(f"BŁĄD: {path}: brak kolumn {miss} (to nie jest h5 v2.0?)")
        trig = {c: f"trigPassed_{c}" for c in TRIG if f"trigPassed_{c}" in names}
        if not is_mc and len(trig) != len(TRIG):
            sys.exit(f"BŁĄD: {path}: brak bitów trygera {sorted(set(TRIG) - set(trig))}")
        need = sorted({c for c in F.raw_columns(cm, 3) + [dt_col] + list(rawx.values()) + list(EV.values()) + list(trig.values()) if c})
        for raw in read_file(path, need, a.chunk):
            if a.max_rows:
                raw = raw.iloc[:a.max_rows]
            cnt["rows"] += len(raw)
            e = {k: raw[c].to_numpy() for k, c in EV.items()}
            if is_mc:
                keep = np.ones(len(raw), bool)
                clsm = None
            else:
                tA, tC = e["tagA"] == 1, e["tagC"] == 1
                xi_ok = (~tA | (e["xcA"] <= e["xiA"] + 0.02)) & (~tC | (e["xcC"] <= e["xiC"] + 0.02))
                sr = (tA | tC) & (e["npu"] == 0) & xi_ok                     # region sygnału głównego łańcucha — NIE patrzymy
                keep = (e["grl"] == 1) & ~sr
                tp = {c: raw[trig[c]].to_numpy() == 1 for c in TRIG}
                M = np.zeros(len(raw), bool)
                for c in MBTS:
                    M |= tp[c]
                J, A = tp[J12], tp[AN] | tp[AS]
                if a.run == 428770:
                    clsm = [np.ones(len(raw), bool), J & ~M, M, A, A & ~J & ~M, ~J & ~M & ~A]
                else:
                    clsm = [np.ones(len(raw), bool), tp["HLT_mb_mbts_L1MBTS_1"]]
            cell_all = (np.clip(np.searchsorted(NSEL_EDGES, e["nsel"], "right") - 1, 0, len(NSEL_EDGES) - 2) * (len(LEAD_EDGES) - 1)
                        + np.clip(np.searchsorted(LEAD_EDGES, e["lead"], "right") - 1, 0, len(LEAD_EDGES) - 2))
            dt = raw[dt_col].to_numpy()
            d0_tab = None
            for sp in ("D0", "Dplus"):
                nd = int(spc[sp]["n_daughters"])
                r_i = np.flatnonzero((dt == int(spc[sp]["decay_type"])) & keep)
                if r_i.size == 0:
                    continue
                r = raw.iloc[r_i].reset_index(drop=True)
                L = F.to_logical(r, cm, nd)
                m, pt = L["m"].to_numpy(float), L["pt"].to_numpy(float)
                lo, hi, nb = MWIN[sp]
                ok = ((L["Lxy"].to_numpy(float) >= (a.lxy_d0 if sp == "D0" else a.lxy_dp)) & (pt >= 3) & (pt < 30)
                      & (np.abs(L["eta"].to_numpy(float)) < ETA_MAX) & (m >= lo) & (m < hi))
                if not ok.any():
                    continue
                r, L, r_i = r[ok].reset_index(drop=True), L[ok].reset_index(drop=True), r_i[ok]
                L = F.compute(L, spc[sp], units)
                score = np.asarray(ens[sp].score(L)[0], float)
                m, pt = L["m"].to_numpy(float), L["pt"].to_numpy(float)
                pb = np.clip(np.searchsorted(PT_EDGES, pt, "right") - 1, 0, 2)
                cnt[f"{sp}_sel"] += len(r)
                if is_mc:
                    t = r[RAWX["truth"]].to_numpy()
                    keys = fine_cat(sp, t, r[RAWX["hs"]].to_numpy(), [r[RAWX["pdg0"]].to_numpy(), r[RAWX["pdg1"]].to_numpy()],
                                    [r[RAWX["tpdg0"]].to_numpy(), r[RAWX["tpdg1"]].to_numpy()], r[RAWX["Dtpdg"]].to_numpy(),
                                    r[RAWX["chan"]].to_numpy())
                else:
                    keys = [c[r_i] for c in clsm]
                dm = np.full(len(r), np.nan)
                if sp == "Dplus":
                    dm = np.fmin(m - L["m_kpi1"].to_numpy(float), m - L["m_kpi2"].to_numpy(float))
                    veto = ~(dm < 0.16)
                if a.thr_out:
                    s = (keys <= 1) & (np.abs(m - M0[sp]) < SIGWIN[sp]) & (veto if sp == "Dplus" else True)
                    for b in range(3):
                        thr_scores[sp][b].append(score[s & (pb == b)])
                    continue
                th = np.array([[thr[sp][p][w] for w in WPS] for p in PT_NAMES])          # [pt, wp]
                sb = (score[:, None] >= th[pb]).sum(axis=1)
                cell = cell_all[r_i]
                dmm = m - M0[sp]
                mreg = np.where(np.abs(dmm) < MREG[0], 0, np.where((dmm < -MREG[1]) & (dmm > -MREG[2]), 1,
                                np.where((dmm > MREG[1]) & (dmm < MREG[2]), 2, -1)))
                rg = (np.clip(mreg, 0, 2), mreg >= 0)
                imp = bidx(np.abs(L["impact"].to_numpy(float)), *AXES["imp"], clip=True)
                ct = bidx(L["Lxy"].to_numpy(float) * M0[sp] / pt, *AXES[f"ct_{sp}"], clip=True)
                mi = bidx(m, lo, hi, nb)
                one = np.ones(len(r), bool)
                if sp == "D0":
                    H.fill("d0_m", pb, sb, keys, [(cell, one), mi])
                    H.fill("d0_imp", pb, sb, keys, [rg, imp])
                    H.fill("d0_ct", pb, sb, keys, [rg, ct])
                    ta, tb = r[RAWX["trk0"]].to_numpy().astype(np.int64), r[RAWX["trk1"]].to_numpy().astype(np.int64)
                    d0_tab = pd.DataFrame({"uid": e["uid"][r_i], "lo": np.minimum(ta, tb), "hi": np.maximum(ta, tb), "k": ta,
                                           "pb": pb, "sb": sb, "m": m})
                else:
                    H.fill("dp_m", pb, sb, keys, [(cell, one), (veto.astype(np.int64), one), mi])
                    H.fill("dp_imp", pb, sb, keys, [rg, imp], sel=veto)
                    H.fill("dp_ct", pb, sb, keys, [rg, ct], sel=veto)
                    H.fill("dp_dm", pb, sb, keys, [bidx(dm, *AXES["dpdm"])])
                    dsw = np.where((dm > 0.1435) & (dm < 0.1475), 0, np.where((dm > 0.150) & (dm < 0.158), 1, -1))
                    H.fill("dp_mds", pb, sb, keys, [(np.clip(dsw, 0, 1), dsw >= 0), mi])
                    P = [(L[f"d_pt{i}"].to_numpy(float), L[f"d_eta{i}"].to_numpy(float), L[f"d_phi{i}"].to_numpy(float)) for i in range(3)]
                    for j, o in ((1, 2), (2, 1)):                 # córka j (π tego samego znaku) jako K albo p
                        mkk = minv([P[0] + (MK,), P[j] + (MK,)])
                        kkpi = minv([P[0] + (MK,), P[j] + (MK,), P[o] + (MPI,)])
                        pkpi = minv([P[0] + (MK,), P[j] + (MP,), P[o] + (MPI,)])
                        ph = np.where(np.abs(mkk - MPHI) < 0.010, 0, np.where((mkk > 1.035) & (mkk < 1.055), 1, -1))
                        pr = (np.clip(ph, 0, 1), ph >= 0)
                        H.fill("dp_mkk", pb, sb, keys, [bidx(mkk, *AXES["mkk"])])
                        H.fill("dp_kkpi", pb, sb, keys, [pr, bidx(kkpi, *AXES["kkpi"])])
                        H.fill("dp_mphi", pb, sb, keys, [pr, mi])
                        H.fill("dp_pkpi", pb, sb, keys, [bidx(pkpi, *AXES["pkpi"])])
            # ---- tag D*: wiersze decayType 1 ↔ wiersze D0 tej samej pary torów (ta sama hipoteza = „right”, odwrotna = „twin”)
            if a.thr_out or d0_tab is None:
                continue
            s_i = np.flatnonzero((dt == 1) & keep)
            if s_i.size == 0:
                continue
            rs = raw.iloc[s_i]
            ta, tb = rs[RAWX["trk0"]].to_numpy().astype(np.int64), rs[RAWX["trk1"]].to_numpy().astype(np.int64)
            dms = rs[RAWX["dmass"]].to_numpy().astype(float)
            S = pd.DataFrame({"uid": e["uid"][s_i], "lo": np.minimum(ta, tb), "hi": np.maximum(ta, tb), "k": ta, "dm": dms,
                              "md0": rs[cm.get("m")].to_numpy().astype(float) - dms, "si": np.arange(s_i.size)})
            if is_mc:
                t, hs = rs[RAWX["truth"]].to_numpy(), rs[RAWX["hs"]].to_numpy()
                S["kds"] = np.where(((t == 0) | (t == 1)) & (hs == 0), 0, np.where(((t == 0) | (t == 1)) & (hs == 1), 1, 2))
            S = S[(S["dm"] >= AXES["dsdm"][0]) & (S["dm"] < AXES["dsdm"][1])]
            cnt["Dstar_rows"] += len(S)
            J2 = S.merge(d0_tab, on=["uid", "lo", "hi"], suffixes=("", "_d0"))
            for name, sel in (("ds_right", J2["k"] == J2["k_d0"]),
                              ("ds_twin", (J2["k"] != J2["k_d0"]) & (np.abs(J2["md0"] - M0["D0"]) < MREG[0]))):
                X = J2[sel.to_numpy()]
                if len(X) == 0:
                    continue
                cnt[name] += len(X)
                si = X["si"].to_numpy()
                keys = X["kds"].to_numpy().astype(np.int64) if is_mc else [c[s_i][si] for c in clsm]
                mm = X["md0"].to_numpy() if name == "ds_right" else X["m"].to_numpy()
                H.fill(name, X["pb"].to_numpy(), X["sb"].to_numpy(), keys, [bidx(X["dm"].to_numpy(), *AXES["dsdm"]), bidx(mm, *AXES["m10"])])
        print(f"[{a.tag}] {os.path.basename(path)}: {cnt}  {time.time() - t0:.0f}s", flush=True)

    if a.thr_out:
        out = {"_def": "kwantyle score sygnału prawdy (fc SigC/SigB, |m−m0| < okno, D+: po wecie D*), per bin pT; pliki: "
                       + ",".join(os.path.basename(p) for p in a.files), "_n": {}}
        for sp in MD:
            out[sp] = {}
            for b, p in enumerate(PT_NAMES):
                s = np.concatenate(thr_scores[sp][b]) if thr_scores[sp][b] else np.zeros(0)
                out["_n"][f"{sp}:{p}"] = int(s.size)
                out[sp][p] = {w: float(np.quantile(s, 1 - float(w))) if s.size else float("nan") for w in WPS}
        json.dump(out, open(a.thr_out, "w"), indent=1)
        print(json.dumps(out, indent=1))
        return 0
    meta = {"mctmpl": VERSION, "tag": a.tag, "run": a.run, "files": [os.path.basename(p) for p in a.files], "counts": cnt,
            "keys": cls_names, "keys_ds": DS_NAMES if is_mc else cls_names, "groups": GROUPS, "pt": PT_NAMES, "wps": WPS,
            "sb": "liczba przekroczonych progów 0..3", "thresholds": {s: thr[s] for s in MD}, "spec": {k: list(v) for k, v in SPEC.items()},
            "axes": {k: (list(v) if not isinstance(v, str) else v) for k, v in AXES.items()}, "m0": M0, "lxy": [a.lxy_d0, a.lxy_dp],
            "nsel_edges": NSEL_EDGES.tolist(), "lead_edges": LEAD_EDGES.tolist(), "eta_max": ETA_MAX,
            "blind": "dane: GRL ∧ ¬(tag AFP ∧ nPileupVertices==0 ∧ ξ_cal ≤ ξ_AFP+0,02); zdarzenia regionu sygnału nieliczone"}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, meta=np.array(json.dumps(meta)), **H.h)
    print(f"[{a.tag}] zapisano {a.out}: {os.path.getsize(a.out) / 1e6:.2f} MB, {time.time() - t0:.0f}s, {cnt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
