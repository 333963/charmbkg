"""DNET (17.11.2026, projekt B) — wspólne definicje: zestawy cech, konfiguracje modeli (bloki jobów), ramki z h5 v2.0.

Cechy liczone przez charmbkg (migawka src/, jak dlz_data): H = 18 cech XGB H (charmbkg.features); D0: cechy CharmNet (`_cn`, wzory
FEAT_charmnet/scripts/extract.py) i nowe (`_nf`, wzory DLZ_d0/scripts/d0_feats.py — przepisane 1:1, bez zależności od ścieżek FEAT);
D+: DLZ8 (dlzlib.add_cn). Izolacja: bezwzględna (ptcone40/100 [GeV]) i względna (ptcone/pT). Tylko pojedyncze D (bez par D–D).
"""
import importlib.util
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("CHARMBKG_SRC") or os.path.join(HERE, "src"))
from charmbkg import config as C, features as F  # noqa: E402
from charmbkg.kinematics import four_vec, swapped_mass_sp  # noqa: E402
from charmbkg.samples import colmap, dstar_veto_mask  # noqa: E402

_bb = importlib.util.spec_from_file_location("bb", os.path.join(HERE, "ext", "bethe_bloch_charmpy.py"))
BB = importlib.util.module_from_spec(_bb)
_bb.loader.exec_module(BB)

cfg = C.load()
H = F.feature_list(cfg, "H")
SPECIES = ("D0", "Dplus")
BINS = [(3, 5), (5, 8), (8, 30)]
WPS = [0.80, 0.70, 0.50]
FAR = {"D0": 0.100, "Dplus": 0.110}            # dalekie SB danych (tło treningowe „sb” i ε_B danych), jak FEAT_charmnet
REFL_VETO = 0.025                              # D0: |m_swap − m0| > 25 MeV w SB danych
SIGNAL_CATS, BKG_MIN = (0, 1), 8               # truth 0/1 = sygnał; ≥ 8 (i < 100) = tło kombinatoryczne

# ---- zestawy cech
CN_D0 = ["absImpact_cn", "d0AbsZ0_cn", "d1AbsZ0_cn", "dMinPtRatio_cn", "dMaxPtRatio_cn", "apAlpha_cn"]
KIN2 = ["ctsK_nf", "ptAsym_nf"]
KINPT = ["dMinPt_nf", "dMaxPt_nf", "dMinPtRatio_cn", "dMaxPtRatio_cn"]
IP = ["d0prod_nf", "d0sigK_nf", "d0sigPi_nf", "d0sigMin_nf", "d0sigMax_nf", "d0chi2_nf", "dz0_nf", "absImpact_cn", "d0AbsZ0_cn", "d1AbsZ0_cn"]
TIME = ["ctau_nf"]
PID = ["bbK_K_nf", "bbK_pi_nf", "bbPi_pi_nf", "bbPi_K_nf", "pPi_nf", "nptsPi_nf", "dedxPi_nf"]
NF_D0 = KIN2 + ["dMinPt_nf", "dMaxPt_nf", "dR_nf"] + [f for f in IP if f.endswith("_nf")] + TIME + PID
ISOABS = ["ptcone40_iso", "ptcone100_iso"]
ISOREL = ["ptcone40rel_iso", "ptcone100rel_iso"]
DALITZ = ["mKpi1_over_mD_cn", "mKpi2_over_mD_cn"]
DLZ8 = DALITZ + ["dMinPt_cn", "dMinPtRatio_cn", "dMedPtRatio_cn", "dMaxPtRatio_cn", "maxdR_cn", "absImpact_cn"]
ALLF = {"D0": H + CN_D0 + NF_D0 + ISOABS + ISOREL, "Dplus": H + DLZ8 + ISOABS + ISOREL}
_KIN = KIN2 + KINPT + ["dR_nf", "apAlpha_cn"]
SETS = {
    "D0": {
        "H": H,                                 # odniesienie (18 cech, bez izolacji)
        "HX": H + IP + TIME,                    # + parametr zderzenia / czas (nie rekonstruują masy) = DLZ_d0 „NOKIN”
        "HK": H + IP + TIME + _KIN,             # + kinematyka rozpadu (rekonstruuje masę → potrzebna dekorelacja)
        "HKP": H + IP + TIME + _KIN + PID,      # + dE/dx córek względem Bethe–Blocha
        "HKia": H + IP + TIME + _KIN + ISOABS,  # + izolacja bezwzględna (jak CharmNet)
        "HKir": H + IP + TIME + _KIN + ISOREL,  # + izolacja względna
    },
    "Dplus": {"H": H, "HDLZ": H + DLZ8, "HDLZia": H + DLZ8 + ISOABS, "HDLZir": H + DLZ8 + ISOREL},
}
DAU = ["pt", "eta", "phi", "d0", "z0", "dedx", "npts", "charge"]      # blok córek Dd [n, nd, 8] (etap 2)

EVENT_COLS = {"grl": "passGRL_MinBiasAFP", "npu": "nPileupVertices", "nSelBL": "nSelBL", "lead_pt": "leadTrackPtBL", "lb": "lumiBlock",
              "afpA": "afp_tag_A", "afpC": "afp_tag_C", "ntrk500": "nTracksPt500Eta25"}
TRIG_CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2",
               "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12",
               "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2"]
MBTS_CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
               "HLT_mb_sptrk_L1MBTS_1", "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2"]
J12 = "HLT_noalg_L1J12"
EXTRA = ["Dmeson_ptcone40", "Dmeson_ptcone100", "Dmeson_hypo_swap", "Dmeson_maxdR", "Dmeson_mKpi1", "Dmeson_mKpi2", "event_uid",
         "event_entry"] + [f"Dmeson_daughter{i}_charge" for i in range(3)]


# ---- konfiguracje modeli: tag = <gatunek>.<model>.<cechy>.<tło>.<dekorelacja><λ>.s<seed>
#   model: xgb (hiperparametry H) | m64 (MLP 3×64 ReLU, przepis CharmNet) | r128 (MLP rezydualny 128, SiLU + LayerNorm)
#   tło:   mc (mb truth ≥ 8, cały zakres mas) | sb (dalekie SB danych 428770) | sbmc (BCE: SB danych; dekorelacja: tło MC w całym zakresie mas)
#   dek.:  d<λ> = DisCo (dCor², ważony) | m<λ> = MoDe[1] (dopuszcza liniową ε_B(m), gładki indykator) | n = brak
def _c(sp, model, feats, bkg, dec="n", seed=1, full=False):
    return {"sp": sp, "model": model, "feats": feats, "bkg": bkg, "dec": dec, "seed": seed, "full": full,
            "tag": f"{sp}.{model}.{feats}.{bkg}.{dec}.s{seed}"}


LAM = os.environ.get("DN_LAM", "5,25").split(",")      # λ nominalne (≥ 3 seedy); skan: DN_LAM_SCAN
LAM_SCAN = os.environ.get("DN_LAM_SCAN", "1,5,25,100").split(",")
BLOCKS = {
    # A: drzewa (odniesienie + nowe cechy + kontrola „drzewo z kinematyką rzeźbi”)
    "xgbA": [_c("D0", "xgb", "H", "sb", seed=s) for s in (1, 2, 3)] + [_c("D0", "xgb", "H", "mc", seed=s) for s in (1, 2, 3)],
    "xgbB": [_c("D0", "xgb", f, b) for f in ("HX", "HK") for b in ("mc", "sb")] + [_c("D0", "xgb", "HKP", "mc")],
    # B: sama zmiana modelu (cechy H)
    "mlpB1": [_c("D0", "m64", "H", "mc", "n"), _c("D0", "m64", "H", "mc", f"d{LAM[0]}", full=True), _c("D0", "m64", "H", "sb", f"d{LAM[0]}")],
    "mlpB2": [_c("D0", "m64", "H", "mc", f"d{LAM[0]}", seed=s) for s in (2, 3)] + [_c("D0", "r128", "H", "mc", f"d{LAM[0]}")],
    # C: H + kinematyka rozpadu z dekorelacją; skan λ
    "mlpC1": [_c("D0", "r128", "HK", "mc", "n")] + [_c("D0", "r128", "HK", "mc", f"d{lam}", full=True) for lam in LAM_SCAN],
    "mlpC2": [_c("D0", "r128", "HK", "mc", f"d{lam}", seed=s) for lam in LAM for s in (2, 3)],
    "mlpC3": [_c("D0", "r128", "HK", "mc", f"m{LAM[0]}"), _c("D0", "m64", "HK", "mc", f"d{LAM[0]}"), _c("D0", "r128", "HX", "mc", f"d{LAM[0]}"),
              _c("D0", "r128", "HKP", "mc", f"d{LAM[0]}")],
    "mlpC4": [_c("D0", "r128", "HK", "sbmc", f"d{LAM[0]}", full=True), _c("D0", "r128", "HK", "sb", f"d{LAM[0]}"),
              _c("D0", "r128", "HKia", "mc", f"d{LAM[0]}"), _c("D0", "r128", "HKir", "mc", f"d{LAM[0]}")],
    # D+
    "dpA": [_c("Dplus", "xgb", f, b) for f in ("H", "HDLZ") for b in ("mc", "sb")] + [_c("Dplus", "m64", "H", "mc", f"d{LAM[0]}")],
    "dpC": [_c("Dplus", "r128", "HDLZ", "mc", "n")] + [_c("Dplus", "r128", "HDLZ", "mc", f"d{LAM[0]}", seed=s, full=(s == 1)) for s in (1, 2, 3)]
           + [_c("Dplus", "r128", "HDLZ", "sbmc", f"d{LAM[0]}")],
}
REF = {"D0": "D0.xgb.H.sb.n.s1", "Dplus": "Dplus.xgb.H.sb.n.s1"}       # „obecne drzewo XGB H” (tło z SB danych, jak nominał)


def all_configs():
    return [c for b in BLOCKS.values() for c in b]


# ---- hash zdarzeń (jak cnclone)
def mix64(x, salt=0):
    with np.errstate(over="ignore"):
        x = np.asarray(x, np.uint64) + np.uint64(0x9E3779B97F4A7C15) * np.uint64(salt + 1)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return x ^ (x >> np.uint64(31))


def unit_hash(x, salt=0):
    return (mix64(x, salt) >> np.uint64(11)).astype(np.float64) / float(1 << 53)


def fold_of(uid, k, salt=7):
    return (mix64(uid, salt) % np.uint64(k)).astype(np.int8)


# ---- cechy
def _ap_alpha(L):
    a = four_vec(L["d_pt0"], L["d_eta0"], L["d_phi0"], BB.M_K)[:, 1:]
    b = four_vec(L["d_pt1"], L["d_eta1"], L["d_phi1"], BB.M_PI)[:, 1:]
    d = a + b
    u = d / np.linalg.norm(d, axis=1, keepdims=True)
    la, lb = (a * u).sum(1) ** 2, (b * u).sum(1) ** 2
    return np.where(la + lb > 0, (la - lb) / np.where(la + lb > 0, la + lb, 1.0), 0.0)


def add_d0(L, sp):
    """Cechy `_cn` (FEAT extract.py) i `_nf` (DLZ_d0 d0_feats.add; dE/dx surowe, bez przesunięcia MC). Córka 0 = K, 1 = π."""
    A_D0, B_D0, M0 = 0.012, 0.075, 1.86484
    pt = [L[f"d_pt{i}"].to_numpy(float) for i in (0, 1)]
    eta = [np.clip(L[f"d_eta{i}"].to_numpy(float), -5, 5) for i in (0, 1)]
    phi = [L[f"d_phi{i}"].to_numpy(float) for i in (0, 1)]
    d0 = [L[f"d_d0{i}"].to_numpy(float) for i in (0, 1)]
    z0 = [L[f"d_z0{i}"].to_numpy(float) for i in (0, 1)]
    D = L["pt"].to_numpy(float)
    o = {}
    o["absImpact_cn"] = np.abs(L["impact"].to_numpy(float))
    o["d0AbsZ0_cn"], o["d1AbsZ0_cn"] = np.abs(z0[0]), np.abs(z0[1])
    o["dMinPtRatio_cn"], o["dMaxPtRatio_cn"] = np.minimum(pt[0], pt[1]) / D, np.maximum(pt[0], pt[1]) / D
    o["apAlpha_cn"] = _ap_alpha(L)
    o["ctsK_nf"] = L["cos_theta_star"].to_numpy(float)
    o["ptAsym_nf"] = (pt[0] - pt[1]) / (pt[0] + pt[1])
    o["dMinPt_nf"], o["dMaxPt_nf"] = np.minimum(pt[0], pt[1]), np.maximum(pt[0], pt[1])
    dphi = np.abs(phi[0] - phi[1])
    dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    o["dR_nf"] = np.hypot(eta[0] - eta[1], dphi)
    o["d0prod_nf"] = d0[0] * d0[1]
    sig = [np.hypot(A_D0, B_D0 / (pt[i] * np.sqrt(1 / np.cosh(eta[i])))) for i in (0, 1)]
    s = [np.abs(d0[i]) / sig[i] for i in (0, 1)]
    o["d0sigK_nf"], o["d0sigPi_nf"] = s[0], s[1]
    o["d0sigMin_nf"], o["d0sigMax_nf"] = np.minimum(s[0], s[1]), np.maximum(s[0], s[1])
    o["d0chi2_nf"] = s[0] ** 2 + s[1] ** 2
    o["dz0_nf"] = np.abs(z0[0] - z0[1])
    o["ctau_nf"] = L["Lxy"].to_numpy(float) * M0 / D
    p = [pt[i] * np.cosh(eta[i]) for i in (0, 1)]
    for i, nm in ((0, "K"), (1, "Pi")):
        dedx = L[f"d_dedx{i}"].to_numpy(float)
        npts = L[f"d_npts{i}"].to_numpy(float)
        v = (npts > 2) & (dedx > 0.1) & (dedx < 8.0) & (p[i] > 0.1) & (p[i] < 50)
        for hyp, mass in (("K", BB.M_K), ("pi", BB.M_PI)):
            x = np.full(len(L), np.nan)
            if v.any():
                x[v] = BB.distance(mass, p[i][v], dedx[v])
            o[f"bb{nm}_{hyp}_nf"] = x
    o["pPi_nf"] = np.log(p[1])
    o["nptsPi_nf"] = L["d_npts1"].to_numpy(float)
    o["dedxPi_nf"] = L["d_dedx1"].to_numpy(float)
    for k, v in o.items():
        L[k] = np.asarray(v, np.float32)
    L["m_swap"] = swapped_mass_sp(L, sp)
    L["ok"] = True
    return L


def add_dplus(L, raw, sp):
    pt = np.column_stack([L[f"d_pt{i}"] for i in range(3)])
    ps = np.sort(pt, axis=1)
    D = L["pt"].to_numpy(float)
    L["absImpact_cn"] = np.abs(L["impact"].to_numpy(float)).astype(np.float32)
    L["dMinPt_cn"] = ps[:, 0]
    L["dMinPtRatio_cn"], L["dMedPtRatio_cn"], L["dMaxPtRatio_cn"] = ps[:, 0] / D, ps[:, 1] / D, ps[:, 2] / D
    L["maxdR_cn"] = raw["Dmeson_maxdR"].to_numpy()
    L["mKpi1_over_mD_cn"] = raw["Dmeson_mKpi1"].to_numpy() * 1e-3 / L["m"]
    L["mKpi2_over_mD_cn"] = raw["Dmeson_mKpi2"].to_numpy() * 1e-3 / L["m"]
    L["m_swap"] = np.nan
    L["ok"] = dstar_veto_mask(L, sp)
    return L


def ftag_of(path):
    m = re.search(r"_(\d{6})_[A-Za-z]", os.path.basename(path))
    return int(m.group(1)) if m else 0


def read_species(path, spn, lxy_min, extra_cols=()):
    """Kandydaci gatunku z jednego h5 v2.0 (pT 3–30, zakres mas, L_xy ≥ lxy_min) → (ramka L z cechami ALLF + kolumny pomocnicze, wiersze h5)."""
    import h5py
    import pandas as pd
    sp = cfg["species"][spn]
    nd = int(sp["n_daughters"])
    cm = colmap(cfg, path, nd)
    cols = list(dict.fromkeys(F.raw_columns(cm, nd) + EXTRA + list(extra_cols)))
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        cols = [c for c in cols if c in names]
        pre = ds.fields([cm["decay_type"], cm["pt"], cm["m"], cm["Lxy"]])[:]
        ptv, mv = pre[cm["pt"]], pre[cm["m"]]
        sel = (pre[cm["decay_type"]] == sp["decay_type"]) & (ptv >= 3) & (ptv < 30) & (mv >= sp["mass_range"][0]) & (mv < sp["mass_range"][1]) \
            & (pre[cm["Lxy"]] >= lxy_min)
        idx = np.flatnonzero(sel)
        arr = ds.fields(cols)[:][idx] if len(idx) else ds.fields(cols)[:0]
    raw = pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
    L = F.to_logical(raw, cm, nd)
    L["weight"] = L["weight"].fillna(1.0)
    L = F.compute(L, {**sp, "n_daughters": nd}, cfg.get("units"))
    L = add_d0(L, sp) if spn == "D0" else add_dplus(L, raw, sp)
    pc40, pc100 = raw["Dmeson_ptcone40"].to_numpy(float), raw["Dmeson_ptcone100"].to_numpy(float)
    L["ptcone40_iso"], L["ptcone100_iso"] = pc40, pc100
    L["ptcone40rel_iso"], L["ptcone100rel_iso"] = pc40 / L["pt"].to_numpy(float), pc100 / L["pt"].to_numpy(float)
    L["uid"] = raw["event_uid"].to_numpy().astype(np.uint64)
    L["entry"] = raw["event_entry"].to_numpy().astype(np.int64) if "event_entry" in raw else -1
    L["swap"] = raw["Dmeson_hypo_swap"].to_numpy().astype(np.int8) if "Dmeson_hypo_swap" in raw else -1
    for i in range(nd):
        L[f"d_charge{i}"] = raw[f"Dmeson_daughter{i}_charge"].to_numpy(float) if f"Dmeson_daughter{i}_charge" in raw else np.nan
    for c in extra_cols:
        if c in raw:
            L["x_" + c] = raw[c].to_numpy()
    return L, idx


def table(L, rows, spn, tag, ftag, prefix, ev, tp):
    """Ramka → dict tablic numpy z prefiksem (a_/p_)."""
    nd = int(cfg["species"][spn]["n_daughters"])
    f4 = lambda c: L[c].to_numpy(np.float32)  # noqa: E731
    with np.errstate(invalid="ignore"):
        X = np.stack([L[c].to_numpy(np.float64) for c in ALLF[spn]], axis=1).astype(np.float32)
    Dd = np.stack([np.stack([L[{"pt": "d_pt", "eta": "d_eta", "phi": "d_phi", "d0": "d_d0", "z0": "d_z0", "dedx": "d_dedx", "npts": "d_npts",
                                "charge": "d_charge"}[v] + str(i)].to_numpy(np.float32) for v in DAU], axis=1) for i in range(nd)], axis=1)
    o = {"uid": L["uid"].to_numpy(np.uint64), "run": L["run"].to_numpy(np.int64), "event": L["event"].to_numpy().astype(np.int64),
         "row": np.asarray(rows, np.int64), "ftag": np.full(len(L), ftag, np.int32), "entry": L["entry"].to_numpy(np.int64),
         "m": L["m"].to_numpy(np.float32), "pt": f4("pt"), "eta": f4("eta"), "phi": f4("phi"), "w": f4("weight"),
         "truth": (L["truth"].to_numpy(float).astype(np.int16) if tag != "data" else np.full(len(L), -1, np.int16)),
         "swap": L["swap"].to_numpy(np.int8), "m_swap": f4("m_swap"), "ok": L["ok"].to_numpy(bool), "X": X, "Dd": Dd}
    for k, c in ev.items():
        if "x_" + c in L:
            o["e_" + k] = L["x_" + c].to_numpy(np.float32 if k == "lead_pt" else np.int32)
    for c, col in tp.items():
        if "x_" + col in L:
            o["e_tp_" + c] = L["x_" + col].to_numpy(np.int8)
    return {prefix + k: v for k, v in o.items()}


def auc(score, y):
    from scipy.stats import rankdata
    y = np.asarray(y, bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(score)
    return float((r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))
