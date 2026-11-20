#!/usr/bin/env python3
"""TEMPLATE/common.py — wspólne wejście dla make_template.py i jump_tests.py (projekt B, AFP-MIX v2, WSTĘPNE).

Selekcja zdarzeń i kandydatów jest KOPIĄ logiki ``AFP_mix/run_data.py::analyse_sample`` (v2) — funkcje pomocnicze importowane
stamtąd, metoda z ``charmbkg.afp_mix`` (bez zmian w obu). Nic tu nie liczy uzysków sygnału ani fitów masy.
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
AFPMIX = os.path.dirname(HERE)
INDIR = os.path.join(AFPMIX, "runs", "lxplus")
sys.path.insert(0, AFPMIX)
sys.path.insert(0, os.path.expanduser("~/Projects/charm/charm-bkg"))
import run_data as RD            # noqa: E402
from charmbkg import afp_mix as M  # noqa: E402

SAMPLES = ("data_427929", "data_428770")
EV_COLS = ["event_number", "lumiBlock", "actual_mu", "passGRL_MinBiasAFP", "nPileupVertices", "nSelBL", "xi_cal_A", "xi_cal_C",
           "afp_nProtons_A", "afp_nProtons_C", "afp_xi1_A", "afp_xi1_C", "afp_xi2_A", "afp_xi2_C", "afp_tag_A", "afp_tag_C"]


def load(s, with_cands=True):
    """→ dict: Dn (dawcy = wszystkie wybrane zdarzenia, nazwy wewnętrzne + hasD, trig, b1), T (biorcy = wybrane z ≥1 kandydatem),
    cands (DataFrame kandydatów w biorcach: t = indeks biorcy, species, m, dm, pt), meta, cfg."""
    meta = json.load(open(os.path.join(INDIR, f"{s}.meta.json")))
    import pyarrow.parquet as pq
    names = pq.ParquetFile(os.path.join(INDIR, f"{s}_events.parquet")).schema_arrow.names
    ev = pd.read_parquet(os.path.join(INDIR, f"{s}_events.parquet"), columns=EV_COLS + [c for c in names if c.startswith("trigPassed_")])
    ev["run_number"] = np.uint32(meta["run_number"])
    n0 = len(ev)
    grl = ev["passGRL_MinBiasAFP"].to_numpy() == 1
    nonafp, afp, _, _ = RD.trigger_classes(ev)
    npu0 = ev["nPileupVertices"].to_numpy() == 0
    sel = grl & nonafp & npu0
    E = RD.standardize_events(ev)
    tcl, mbm = RD.trig_class_and_mb(ev)
    evno = ev["event_number"].to_numpy()
    del ev
    ca = pd.read_parquet(os.path.join(INDIR, f"{s}_cands.parquet"), columns=["event_number", "Dmeson_species", "Dmeson_mass", "Dmeson_delta_mass", "Dmeson_pt"])
    spec = ca["Dmeson_species"].to_numpy()
    crow = pd.Index(evno).get_indexer(ca["event_number"].to_numpy())
    assert (crow >= 0).all()
    csel = np.isin(spec, [sp["code"] for sp in RD.SPECIES.values()])
    has = np.zeros(n0, bool)
    has[crow[csel]] = True
    E["hasD"] = has
    E["trig"] = tcl
    mb_frac = float(mbm[sel].mean())
    E["b1"] = mbm & (~has if mb_frac >= 0.5 else True)
    t_rows = np.flatnonzero(sel & has)
    Dn = E.loc[sel].reset_index(drop=True)
    T = E.iloc[t_rows].reset_index(drop=True)
    del E
    out = {"sample": s, "run": int(meta["run_number"]), "meta": meta, "Dn": Dn, "T": T, "cfg": M.MixConfig(n_min_pool=50, tag_def="nside")}
    if with_cands:
        row2t = np.full(n0, -1, np.int64)
        row2t[t_rows] = np.arange(len(t_rows))
        ct = row2t[crow]
        ok = (ct >= 0) & csel
        out["cands"] = pd.DataFrame({"t": ct[ok], "species": spec[ok], "m": ca["Dmeson_mass"].to_numpy()[ok],
                                     "dm": ca["Dmeson_delta_mass"].to_numpy()[ok], "pt": ca["Dmeson_pt"].to_numpy()[ok]})
    return out


def jsonable(o):
    return RD.jsonable(o)
