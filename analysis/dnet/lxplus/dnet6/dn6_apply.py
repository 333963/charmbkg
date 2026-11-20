#!/usr/bin/env python3
"""DNET etap 6 — dn4_apply.py z maską wolnego pionu dla modeli, których meta ma pole „mask” (kod etapu 4 bez zmian: podmiana w pamięci
funkcji wczytania wierszy i liczenia score). Modele bez maski (drzewo H, pełna sieć) liczone dokładnie jak w etapie 4. Argumenty jak dn4_apply.py."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_Q = os.path.dirname(HERE)
sys.path.insert(0, os.environ.get("DN4") or os.path.join(_Q, "dnet4"))
sys.path.insert(0, HERE)
import dn4_apply as A4  # noqa: E402
import dn4_lib as L4  # noqa: E402
import dn6_mask as MK  # noqa: E402

_REG, _CUR = {}, {}
_load_one0, _load_models0, _score_rows0 = A4.load_one, L4.load_models, L4.score_rows


def _load_one(ddir, edir, s, sp):
    r = _load_one0(ddir, edir, s, sp)
    if r is not None and sp == "D0":
        T = r[0]
        Em, nrm = MK.apply(T["E"], T["Dd"], T["pt"], T["eta"], T["phi"], T["m"])
        _REG.clear(); _REG[id(T["E"])] = Em
        print(f"[{s}] maska D*: wierszy z usuniętym torem {float((nrm > 0).mean()):.4f} (Δm < {MK.DM_MAX}, znak {MK.SIGN})", flush=True)
    return r


def _load_models(mroot, tag):
    meta, ms = _load_models0(mroot, tag)
    _CUR["mask"] = meta.get("mask")
    if _CUR["mask"] and (abs(_CUR["mask"]["dm_max"] - MK.DM_MAX) > 1e-9 or int(_CUR["mask"]["sign"]) != MK.SIGN):
        raise RuntimeError(f"{tag}: maska modelu {_CUR['mask']} ≠ bieżąca (DM_MAX {MK.DM_MAX}, znak {MK.SIGN})")
    return meta, ms


def _score_rows(kind, models, uid, X, Dd=None, E=None, eta=None, phi=None, **kw):
    if kind == "tr" and _CUR.get("mask"):
        E = _REG[id(E)]                               # brak wpisu = błąd (maska wymagana, a wiersze nie przeszły przez _load_one)
    return _score_rows0(kind, models, uid, X, Dd, E, eta, phi, **kw)


A4.load_one, L4.load_models, L4.score_rows = _load_one, _load_models, _score_rows

if __name__ == "__main__":
    sys.exit(A4.main())
