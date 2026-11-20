"""B4 — SF_rel z B2 (D0 z tagu D*) jako **systematyka ε_S** w ``charmtag sf`` (projekt A) (N) i w parach (U), nie nominał.

Plik: ``inputs/sf_rel_b2.json`` (obok ``charmbkg/data/``; ścieżka nadpisywalna ``sf.sf_rel_json``).
API:
* ``sf_rel(species, pt_bin, wp) -> (val, unc)`` — D⁺ = transfer z D0 ⊕ 5 % (N);
* ``eps_s_syst(eps_S, species, pt_bin, wp) -> dict`` — nominał bez zmian, ``syst_rel`` = unc ⊕ |val−1|·0 (centralne
  przesunięcie NIE jest stosowane; raportowane jako ``central_shift_rel``);
* ``annotate_rows(rows, species, key_bin="bin", key_wp="wp", key_eps="eps_S")`` — dopisuje pola ``sf_rel_b2`` do wierszy
  wyniku ``sf``/raportu (addytywnie).
Hak (3 linie, domyślnie wyłączony): w ``pipeline.cmd_sf`` (lub tam, gdzie powstaje tabela SF per bin×WP):
``if cfg.get("sf", {}).get("sf_rel_syst"): from charmbkg import sf_rel; rows = sf_rel.annotate_rows(rows, species)``.
"""

from __future__ import annotations

import json
import math
import os
import re
from functools import lru_cache

_DEFAULT_PATHS = (
    os.path.join(os.path.dirname(__file__), "data", "sf_rel_b2.json"),
    os.path.join(os.path.dirname(os.path.dirname(__file__)), "inputs", "sf_rel_b2.json"),
)


def _norm_species(sp: str) -> str:
    s = str(sp).replace("+", "plus").replace("⁺", "plus")
    return "Dplus" if s.lower().startswith("dp") else "D0"


def _norm_bin(b) -> str:
    """'8-30', '8_30', (8, 30), '8–30', 'pt8_30' → '8-30'."""
    if isinstance(b, (tuple, list)) and len(b) == 2:
        lo, hi = b
        return f"{round(float(lo))}-{round(float(hi))}"
    s = str(b).replace("–", "-").replace("_", "-")
    m = re.findall(r"\d+(?:\.\d+)?", s)
    if len(m) >= 2:
        return f"{round(float(m[0]))}-{round(float(m[1]))}"
    return s


def _norm_wp(wp) -> str:
    s = re.sub(r"[^\d.]", "", str(wp))
    v = float(s)
    if v <= 1.0:
        v *= 100.0
    return str(round(v))


@lru_cache(maxsize=4)
def load(path: str | None = None) -> dict:
    cands = [path] if path else list(_DEFAULT_PATHS)
    for p in cands:
        if p and os.path.exists(p):
            with open(p) as f:
                js = json.load(f)
            js["_path"] = p
            return js
    raise FileNotFoundError(f"sf_rel_b2.json nie znaleziony: {cands}")


def sf_rel(species: str, pt_bin, wp, path: str | None = None) -> tuple[float, float]:
    js = load(path)
    sp = _norm_species(species)
    b = _norm_bin(pt_bin)
    w = _norm_wp(wp)
    extra = 0.0
    tab = js.get(sp)
    if tab is None or "transfer_from" in tab:
        extra = float((tab or {}).get("extra_unc", 0.05))
        tab = js[(tab or {}).get("transfer_from", "D0")]
    if b not in tab:
        raise KeyError(f"bin {b!r} nie w sf_rel_b2.json ({sorted(k for k in tab if '-' in k)})")
    row = tab[b]
    if w not in row:
        # WP95 = referencja (1.0, 0); inne WP → najbliższy znany
        if w == "95":
            return 1.0, math.hypot(0.0, extra)
        known = sorted((abs(int(k) - int(w)), k) for k in row if k.isdigit())
        w = known[0][1]
    val, unc = float(row[w][0]), float(row[w][1])
    return val, math.hypot(unc, extra)


def eps_s_syst(eps_S: float, species: str, pt_bin, wp, path: str | None = None) -> dict:
    val, unc = sf_rel(species, pt_bin, wp, path)
    return {
        "eps_S_nominal": float(eps_S),
        "sf_rel": val,
        "sf_rel_unc": unc,
        "syst_rel": unc,
        "eps_S_syst_abs": float(eps_S) * unc,
        "central_shift_rel": val - 1.0,
        "central_shift_applied": False,
        "source": "B2 sf_rel_b2.json",
    }


def annotate_rows(rows, species: str, key_bin: str = "bin", key_wp: str = "wp", key_eps: str = "eps_S", path=None):
    """Addytywnie dopisuje ``sf_rel_b2`` do każdego wiersza (dict) z binem, WP i ε_S. Wiersze bez klucza pomijane."""
    for r in rows:
        try:
            r["sf_rel_b2"] = eps_s_syst(float(r[key_eps]), species, r[key_bin], r[key_wp], path)
        except (KeyError, TypeError, ValueError, FileNotFoundError) as e:
            r["sf_rel_b2"] = {"skipped": str(e)}
    return rows
