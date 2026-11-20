"""FEAT_charmnet: ekstrakcja kandydatów D0/D+ z h5 v1.2 → ramki z 18 cechami H (charmbkg.features) + nowe cechy CharmNet.

Użycie: python extract.py SAMPLE_DIR OUT.pkl.gz SPECIES   (SPECIES: D0 | Dplus)
Klucz zdarzenia `evkey` = file_idx·10⁷ + numer bloku kolejnych wierszy o tym samym event_number
(fd: event_number powtarza się między plikami AOD wewnątrz jednej części; zweryfikowane: bloki = unikalne (segment, event)).
Charmbkg używany jako biblioteka (bez edycji).
"""
import glob
import importlib.util
import os
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.environ.get("CHARMBKG_SRC", "/home/k/Projects/charm/charm-bkg"))
from charmbkg import config as C, features as F  # noqa: E402
from charmbkg.samples import colmap, dstar_veto_mask  # noqa: E402
from charmbkg.kinematics import swapped_mass_sp, four_vec  # noqa: E402

_bb = importlib.util.spec_from_file_location("bb_charmpy", os.path.join(os.path.dirname(__file__), "ext", "bethe_bloch_charmpy.py"))
BB = importlib.util.module_from_spec(_bb)
_bb.loader.exec_module(BB)

EXTRA = ["Dmeson_ptcone40", "Dmeson_ptcone60", "Dmeson_ptcone80", "Dmeson_ptcone100", "Dmeson_maxdR",
         "Dmeson_hypo_swap", "Dmeson_impact"]
DEDX_MC_SHIFT = -0.023  # CharmNet: MC pixel dE/dx przesunięte przed DistPerp (dane bez zmian)


def distperp(p, dedx, npts, mass):
    out = np.full(len(p), np.nan)
    v = (npts > 2) & (dedx > 0.1) & (dedx < 8.0) & (p > 0.1) & (p < 50.0)
    if v.any():
        out[v] = BB.distance_perp(mass, p[v], dedx[v])
    return out


def ap_alpha(L):
    a = four_vec(L["d_pt0"], L["d_eta0"], L["d_phi0"], BB.M_K)[:, 1:]
    b = four_vec(L["d_pt1"], L["d_eta1"], L["d_phi1"], BB.M_PI)[:, 1:]
    d = a + b
    u = d / np.linalg.norm(d, axis=1, keepdims=True)
    la, lb = (a * u).sum(1) ** 2, (b * u).sum(1) ** 2
    return np.where(la + lb > 0, (la - lb) / (la + lb), 0.0)


def main(sdir, out, species):
    cfg = C.load()
    sp = cfg["species"][species]
    nd = int(sp["n_daughters"])
    is_data = "data" in os.path.basename(os.path.normpath(sdir))
    files = sorted(glob.glob(os.path.join(sdir, "*_none.h5")))
    parts = []
    for fi, path in enumerate(files):
        cm = colmap(cfg, path, nd)
        cols = list(dict.fromkeys(F.raw_columns(cm, nd) + EXTRA + ["Dmeson_mKpi1", "Dmeson_mKpi2"]))
        with h5py.File(path, "r") as f:
            ds = f["events"]
            names = set(ds.dtype.names)
            cols = [c for c in cols if c in names]
            arr = ds.fields(cols)[:]
        raw = pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
        ev = raw[cm["event"]].to_numpy().astype(np.int64)
        blk = np.r_[0, np.cumsum(np.diff(ev) != 0)]
        sel = (raw[cm["decay_type"]].to_numpy() == sp["decay_type"]) & (raw[cm["pt"]].to_numpy() >= 3) & (raw[cm["pt"]].to_numpy() < 30)
        m = raw[cm["m"]].to_numpy()
        sel &= (m >= sp["mass_range"][0]) & (m < sp["mass_range"][1])
        raw, blk = raw[sel].reset_index(drop=True), blk[sel]
        L = F.to_logical(raw, cm, nd)
        L["weight"] = L["weight"].fillna(1.0)
        L = F.compute(L, {**sp, "n_daughters": nd}, cfg.get("units"))
        L["evkey"] = fi * 10_000_000 + blk
        L["file"] = fi
        for c in EXTRA:
            if c in raw:
                L[c.replace("Dmeson_", "")] = raw[c].to_numpy()
        # --- nowe cechy CharmNet (nazwy *_cn) ---
        L["ptcone40_cn"], L["ptcone100_cn"] = L["ptcone40"], L["ptcone100"]
        L["ptcone60_cn"], L["ptcone80_cn"] = L["ptcone60"], L["ptcone80"]
        L["absImpact_cn"] = np.abs(L["impact"])
        pt = np.column_stack([L[f"d_pt{i}"] for i in range(nd)])  # GeV (po compute)
        ps = np.sort(pt, axis=1)
        L["dMinPt_cn"] = ps[:, 0]
        L["dMinPtRatio_cn"] = ps[:, 0] / L["pt"]
        L["dMaxPtRatio_cn"] = ps[:, -1] / L["pt"]
        if nd == 3:
            L["dMedPtRatio_cn"] = ps[:, 1] / L["pt"]
        az = np.abs(np.column_stack([L[f"d_z0{i}"] for i in range(nd)]))
        L["dMinAbsZ0_cn"] = az.min(1)
        for i in range(nd):
            L[f"d{i}AbsZ0_cn"] = az[:, i]
        if species == "D0":
            L["apAlpha_cn"] = ap_alpha(L)
            L["m_swap"] = swapped_mass_sp(L, sp)
        else:
            L["maxdR_cn"] = L["maxdR"]
            # mKpi w MeV, m w GeV → /1000 (poprawka M5 jednostek)
            L["mKpi1_over_mD_cn"] = raw["Dmeson_mKpi1"].to_numpy() * 1e-3 / L["m"]
            L["mKpi2_over_mD_cn"] = raw["Dmeson_mKpi2"].to_numpy() * 1e-3 / L["m"]
            L["dstar_ok"] = dstar_veto_mask(L, sp)
            masses = (BB.M_K, BB.M_PI, BB.M_PI)
            for i in range(nd):
                dedx = L[f"d_dedx{i}"].to_numpy(float) + (0.0 if is_data else DEDX_MC_SHIFT)
                pp = L[f"d_pt{i}"].to_numpy(float) * np.cosh(np.clip(L[f"d_eta{i}"].to_numpy(float), -10, 10))
                dp = distperp(pp, dedx, L[f"d_npts{i}"].to_numpy(float), masses[i])
                L[f"d{i}_dEdxValid_cn"] = np.isfinite(dp).astype(np.float32)
                L[f"d{i}_dEdxDistPerp_cn"] = np.nan_to_num(dp, nan=0.0)
        for c in L.columns:
            if L[c].dtype == np.float64 and c not in ("m", "m_swap"):
                L[c] = L[c].astype(np.float32)
        parts.append(L)
        print(f"{path}: {len(L)} kand.", flush=True)
    L = pd.concat(parts, ignore_index=True)
    L["sample"] = os.path.basename(os.path.normpath(sdir))
    L.to_pickle(out)
    print("zapisano", out, len(L))


if __name__ == "__main__":
    main(*sys.argv[1:4])
