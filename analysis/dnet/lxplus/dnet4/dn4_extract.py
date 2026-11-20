#!/usr/bin/env python3
"""DNET etap 4 — ekstrakcja h5 v2.0 → tabela dn_<próbka>.npz jak ../dnet1/dn_extract.py (ten sam kod, te same cechy), z dwiema różnicami
ustawianymi zmiennymi środowiska (kod etapu 1 nie jest zmieniany — podmieniana jest tylko funkcja wyboru wierszy w pamięci):
  DN4_PT="lo,hi"            okno pT kandydata [GeV) zamiast twardego 3–30 (bin 2–3 GeV: produkcja dpt2to3)
  DN4_MR_D0="lo,hi", DN4_MR_Dplus="lo,hi"   zakres mas zamiast cfg mass_range (D0: szerszy zakres do wariantu zakresu fitu)
Argumenty: jak dn_extract.py (--tag data|fd|mb --out OUT.npz [--frac-d0 f] [--frac-dp f] [--lxy-min x] [--species ...] PLIK.h5 ...)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dn4_lib as L4  # noqa: E402,F401  (sys.path → ../dnet1)
import dn_feat as DF  # noqa: E402
from charmbkg import features as F  # noqa: E402
from charmbkg.samples import colmap  # noqa: E402


def _rng(name, default):
    v = os.environ.get(name, "")
    return tuple(float(x) for x in v.split(",")) if v else tuple(default)


def read_species_win(path, spn, lxy_min, extra_cols=()):
    """Kopia dn_feat.read_species z oknem pT i zakresem mas ze środowiska (reszta 1:1)."""
    import h5py
    import pandas as pd
    cfg = DF.cfg
    sp = cfg["species"][spn]
    nd = int(sp["n_daughters"])
    plo, phi = _rng("DN4_PT", (3.0, 30.0))
    mlo, mhi = _rng("DN4_MR_" + spn, sp["mass_range"])
    cm = colmap(cfg, path, nd)
    cols = list(dict.fromkeys(F.raw_columns(cm, nd) + DF.EXTRA + list(extra_cols)))
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        cols = [c for c in cols if c in names]
        pre = ds.fields([cm["decay_type"], cm["pt"], cm["m"], cm["Lxy"]])[:]
        ptv, mv = pre[cm["pt"]], pre[cm["m"]]
        sel = (pre[cm["decay_type"]] == sp["decay_type"]) & (ptv >= plo) & (ptv < phi) & (mv >= mlo) & (mv < mhi) & (pre[cm["Lxy"]] >= lxy_min)
        idx = np.flatnonzero(sel)
        arr = ds.fields(cols)[:][idx] if len(idx) else ds.fields(cols)[:0]
    raw = pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
    L = F.to_logical(raw, cm, nd)
    L["weight"] = L["weight"].fillna(1.0)
    L = F.compute(L, {**sp, "n_daughters": nd}, cfg.get("units"))
    L = DF.add_d0(L, sp) if spn == "D0" else DF.add_dplus(L, raw, sp)
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
    if len(L):
        pt = L["pt"].to_numpy(float)
        if pt.min() < plo - 1e-3 or pt.max() >= phi + 1e-3:
            raise RuntimeError(f"{path} {spn}: pT poza oknem [{plo}, {phi}): {pt.min():.3f}–{pt.max():.3f}")
    return L, idx


def main(argv=None):
    DF.read_species = read_species_win
    import dn_extract as DX
    print(f"dn4_extract: pT {_rng('DN4_PT', (3.0, 30.0))}, zakres mas D0 {_rng('DN4_MR_D0', DF.cfg['species']['D0']['mass_range'])}, "
          f"D+ {_rng('DN4_MR_Dplus', DF.cfg['species']['Dplus']['mass_range'])}", flush=True)
    return DX.main(argv)


if __name__ == "__main__":
    sys.exit(main())
