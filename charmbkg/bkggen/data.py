"""bkggen: wczytanie kandydatów (h5 v1.1 przez charmbkg albo podpróbka npz z extract_sub.py) → ramka logiczna + cechy.

Konwencja jak w charmbkg: m, pt w GeV; d_pt* przeliczone do GeV (units.d_pt); cechy z features.compute.
"""
from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

from .. import features as F
from ..io import resolve_columns
from ..samples import mark_swaps

from ..constants import M_K, M_PI  # noqa: F401
from .. import kinematics as _kin
BKG_TRUTH_MIN = 8  # truth_category ≥ 8 = tło kombinatoryczne (MC)


def _npz_frame(paths):
    parts = []
    for p in paths:
        z = np.load(p)
        parts.append(pd.DataFrame({k: z[k] for k in z.files}))
    return pd.concat(parts, ignore_index=True)


def load_npz(cfg: dict, pattern: str, species: str = "D0", decay_type: int | None = None, nd: int | None = None,
             tag: str = "") -> pd.DataFrame:
    """Podpróbka npz (kolumny surowe h5) → ramka logiczna z cechami charmbkg."""
    paths = sorted(glob.glob(os.path.expanduser(pattern)))
    if not paths:
        raise FileNotFoundError(pattern)
    raw = _npz_frame(paths)
    sp = cfg["species"][species]
    dt = int(decay_type if decay_type is not None else sp["decay_type"])
    nd = int(nd or sp["n_daughters"])
    raw = raw[raw["Dmeson_decayType"].to_numpy() == dt].reset_index(drop=True)
    cm = resolve_columns(list(raw.columns), cfg["columns"], nd)
    L = F.to_logical(raw, cm, nd)
    L["weight"] = L["weight"].fillna(1.0)
    L = F.compute(L, {**sp, "n_daughters": nd}, cfg.get("units"))
    L["row"] = raw["_row"].to_numpy() if "_row" in raw else np.arange(len(raw))
    for extra in ("Dmeson_soft_pion_d0", "Dmeson_soft_pion_dR", "Dmeson_maxdR"):
        if extra in raw:
            L[extra.replace("Dmeson_", "")] = raw[extra].to_numpy()
    if nd == 2 and dt == sp["decay_type"]:
        mark_swaps(L, sp, cfg)
    L["sample"] = tag
    return L


def load_h5(cfg: dict, key: str, species: str = "D0", decay_type=None, nd=None, select=None, keep_frac=1.0):
    """Pełne h5 przez charmbkg (Mac/lxplus)."""
    from ..samples import load_species
    return load_species(cfg, key, species, cache=False, decay_type=decay_type, nd=nd, select=select,
                        keep_frac=keep_frac, tagx="bkggen")


def four_vec(pt, eta, phi, m):
    """→ charmbkg.kinematics.four_vec (alias zgodności)."""
    return _kin.four_vec(pt, eta, phi, m)


def inv_mass(*p4):
    """→ charmbkg.kinematics.inv_mass (alias zgodności)."""
    return _kin.inv_mass(*p4)


def swapped_mass(L: pd.DataFrame, kaon_index=0, pion_index=1) -> np.ndarray:
    """D0: masa z zamienionymi hipotezami K↔π. Implementacja: charmbkg.kinematics."""
    return _kin.swapped_mass(L, kaon_index, pion_index)


def nominal_mass(L, kaon_index=0, pion_index=1):
    return _kin.nominal_mass(L, kaon_index, pion_index)


def is_bkg_truth(L):
    t = L["truth"].to_numpy(dtype=float)
    return (t >= BKG_TRUTH_MIN) & (t < 100)


def is_sig_truth(L, cfg):
    return L["truth"].isin(cfg["signal"]["truth_values"]).to_numpy()


def in_range(L, sp, pt_lo=3.0, pt_hi=30.0):
    lo, hi = sp["mass_range"]
    m, pt = L["m"].to_numpy(), L["pt"].to_numpy()
    return (m >= lo) & (m < hi) & (pt >= pt_lo) & (pt < pt_hi)
