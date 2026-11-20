"""DLZ_dplus/DATA: wspólne definicje (cechy, foldy, progi WP, ekstrakcja z h5) — używane lokalnie i w paczce lxplus Q_lxplus/dlz_data.

Cechy `_cn` liczone IDENTYCZNIE jak w FEAT_charmnet/scripts/extract.py (ten sam wzór), z h5 v1.2 lub v2.0.
"""
import os
import sys

import numpy as np
import pandas as pd

CHARMBKG = os.environ.get("CHARMBKG_SRC") or ("/home/k/Projects/charm/charm-bkg" if os.path.isdir("/home/k/Projects/charm/charm-bkg") else os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
sys.path.insert(0, CHARMBKG)
from charmbkg import config as C, features as F  # noqa: E402
from charmbkg.samples import colmap, dstar_veto_mask, fold_of  # noqa: E402

cfg = C.load()
SP = cfg["species"]["Dplus"]
M0 = SP["m0"]
H = F.feature_list(cfg, "H")
DALITZ = ["mKpi1_over_mD_cn", "mKpi2_over_mD_cn"]
PTRATIO = ["dMinPt_cn", "dMinPtRatio_cn", "dMedPtRatio_cn", "dMaxPtRatio_cn"]
DLZ8 = DALITZ + PTRATIO + ["maxdR_cn", "absImpact_cn"]
VARIANTS = {"H": H, "HDLZ": H + DLZ8, "Hdal": H + DALITZ}
BINS = [(3, 5), (5, 8), (8, 30)]
WPS = [0.95, 0.90, 0.80, 0.70, 0.50]
SID = {"fd_nd": 1, "mb_nd": 2}
EXTRA = ["Dmeson_maxdR", "Dmeson_impact", "Dmeson_mKpi1", "Dmeson_mKpi2"]
V20 = ["passGRL_MinBiasAFP", "trigPassed_HLT_mb_mbts_L1MBTS_1", "trigPassed_HLT_noalg_L1J12", "nSelBL", "average_mu"]


def add_cn(L, raw):
    """8 cech H+DLZ; L = ramka logiczna po F.compute (pT w GeV), raw = surowe kolumny h5."""
    nd = 3
    L["absImpact_cn"] = np.abs(raw["Dmeson_impact"].to_numpy()).astype(np.float32)
    pt = np.column_stack([L[f"d_pt{i}"] for i in range(nd)])
    ps = np.sort(pt, axis=1)
    L["dMinPt_cn"] = ps[:, 0]
    L["dMinPtRatio_cn"] = ps[:, 0] / L["pt"]
    L["dMedPtRatio_cn"] = ps[:, 1] / L["pt"]
    L["dMaxPtRatio_cn"] = ps[:, -1] / L["pt"]
    L["maxdR_cn"] = raw["Dmeson_maxdR"].to_numpy()
    L["mKpi1_over_mD_cn"] = raw["Dmeson_mKpi1"].to_numpy() * 1e-3 / L["m"]
    L["mKpi2_over_mD_cn"] = raw["Dmeson_mKpi2"].to_numpy() * 1e-3 / L["m"]
    return L


def read_h5(path, lxy_min=0.25):
    """Kandydaci D+ (decay_type, pT 3–30, okno masy, L_xy ≥ lxy_min) z jednego h5 → ramka z cechami H + _cn + dstar_ok + kolumny v2.0."""
    import h5py
    nd = 3
    cm = colmap(cfg, path, nd)
    cols = list(dict.fromkeys(F.raw_columns(cm, nd) + EXTRA + V20))
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        cols = [c for c in cols if c in names]
        # najpierw tanie kolumny do selekcji
        pre = ds.fields([cm["decay_type"], cm["pt"], cm["m"], cm["event"]])[:]
        evn = np.asarray(pre[cm["event"]]).astype(np.int64)
        blk_all = np.r_[0, np.cumsum(np.diff(evn) != 0)] if len(evn) else np.zeros(0, np.int64)
        sel = (pre[cm["decay_type"]] == SP["decay_type"])
        ptv, mv = pre[cm["pt"]], pre[cm["m"]]
        scale = 1e-3 if np.nanmedian(ptv[sel][:1000]) > 500 else 1.0  # MeV → GeV (jak F.compute/units)
        sel &= (ptv * scale >= 3) & (ptv * scale < 30) & (mv * scale >= SP["mass_range"][0]) & (mv * scale < SP["mass_range"][1])
        idx = np.flatnonzero(sel)
        arr = ds.fields(cols)[:][idx] if len(idx) else ds.fields(cols)[:0]
    raw = pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
    L = F.to_logical(raw, cm, nd)
    L["weight"] = L["weight"].fillna(1.0)
    L = F.compute(L, {**SP, "n_daughters": nd}, cfg.get("units"))
    L = add_cn(L, raw)
    L["dstar_ok"] = dstar_veto_mask(L, SP)
    import re
    mfi = re.search(r"_(\d{6})_", os.path.basename(path))
    L["evkey"] = (int(mfi.group(1)) if mfi else 0) * 10_000_000 + blk_all[idx]  # MC: blok kolejnych wierszy o tym samym event_number
    for c in V20:
        if c in raw:
            L[c] = raw[c].to_numpy()
    keep = (L["Lxy"].to_numpy() >= lxy_min) & (L["m"].to_numpy() >= SP["mass_range"][0]) & (L["m"].to_numpy() < SP["mass_range"][1]) \
        & (L["pt"].to_numpy() >= 3) & (L["pt"].to_numpy() < 30)
    return L[keep].reset_index(drop=True)


def data_fold(L):
    return fold_of(L["run"].to_numpy(np.int64), L["event"].to_numpy().astype(np.int64), 5)


def thresholds(sS, wS, ptS, wps=WPS):
    """progi per bin pT z sygnału MC OOF (ważone): dict (lo,hi,wp) → próg."""
    t = {}
    for lo, hi in BINS:
        a = (ptS >= lo) & (ptS < hi)
        s, w = sS[a], wS[a]
        o = np.argsort(-s)
        c = np.cumsum(w[o]) / w.sum()
        for wp in wps:
            t[(lo, hi, wp)] = float(s[o][min(np.searchsorted(c, wp), len(o) - 1)])
    return t


def ptbin(pt):
    return np.digitize(pt, [3, 5, 8, 30]) - 1


def load_models(mdir, variant, seed):
    import xgboost as xgb
    out = []
    for k in range(5):
        b = xgb.Booster()
        b.load_model(os.path.join(mdir, f"{variant}_s{seed}", f"fold{k}.json"))
        out.append(b)
    return out


def predict(models, L, feats, fold, nthread=3):
    import xgboost as xgb
    s = np.full(len(L), np.nan, np.float32)
    X = L[feats].to_numpy(np.float32)
    for k, b in enumerate(models):
        sel = fold == k
        if sel.any():
            b.set_param({"nthread": nthread})
            s[sel] = b.predict(xgb.DMatrix(X[sel]), iteration_range=(0, b.best_iteration + 1))
    return s
