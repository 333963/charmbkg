"""charmbkg.s8salad — SALAD (reweighting MC tła → dane w SB, interpolacja do SR) w charmbkg (etap S8).

Dwa zastosowania, oba addytywne (domyślnie wyłączone):
  1) tło treningowe:   training.bkg_source: synthetic, bkg_options.synthetic_model: salad
                       → `training_background_salad(cfg, species)` (wołane z bkgsrc.bkg_synthetic): PRAWDZIWE zdarzenia
                       MC tła (incl_mc truth≥8) w |m−m0| < sideband_min z wagą SALAD w kolumnie `weight`.
  2) systematyka kształtu tła po cięciu:  `charmbkg s8salad -c CFG --species Dplus [--variant H] [--kind sk_mlp]`
                       → ε_B(SR) z MC ważonego SALAD vs ε_B(SB danych) vs ε_B(MC surowe) @WP, per bin pT
                       → <out_dir>/syst/bkg_shape_salad_<species>.json (format jak Y: species/source/definition/entries)
                       oraz C2ST dane vs model w paśmie (inner_band) jako kontrola jakości.
Reweighter: charmbkg.bkggen.salad8.SaladReweighter8 (sklearn; kind sk_mlp|hgb|hgb_nom; 'mlp' = torch S4 gdy jest).
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import numpy as np
import pandas as pd

log = logging.getLogger("charmbkg.s8salad")

DEFAULTS = {"salad_kind": "sk_mlp", "salad_max_n": 400_000, "salad_mlp_iter": 60, "salad_mc_truth_min": 8,
            "salad_region": "inner",
            # S5MEM (05.11): pełna lumi D0 = dziesiątki mln kandydatów → 65 GB RAM na lxplus (HELD, kod 34).
            # Losowa, deterministyczna podpróbka per chunk (samples.load_species keep_frac) do salad_max_rows wierszy
            # na próbkę; ε_B to stosunki → podpróbka jednorodna nie zmienia wartości, stat. 3 M ≫ potrzeb (closure S8 na 0,5 M).
            "salad_max_rows": 3_000_000}   # inner: |m−m0|<sideband_min (SR+pasmo) | sr: |m−m0|<signal_window


def _opts(cfg):
    o = dict(DEFAULTS)
    o.update({k: v for k, v in ((cfg.get("training") or {}).get("bkg_options") or {}).items() if k.startswith("salad_")})
    return o


def _feats_present(df, feats):
    return [f for f in feats if f in df and np.isfinite(df[f].to_numpy(dtype=float)).any()]


def fit_reweighter(cfg, species, MC: pd.DataFrame, SBd: pd.DataFrame, feats, kind=None, seed=None):
    """SALAD: MC (cały zakres mas, kolumny cech + m) vs dane SB. Zwraca (rw, feats_użyte)."""
    from .bkggen.salad8 import SaladReweighter8
    sp = cfg["species"][species]
    o = _opts(cfg)
    kind = kind or o["salad_kind"]
    seed = int(cfg["io"]["seed"] if seed is None else seed)
    feats = _feats_present(MC, _feats_present(SBd, feats))
    m_mc = MC["m"].to_numpy(float)
    sb_mc = np.abs(m_mc - sp["m0"]) > sp["sideband_min"]
    rw = SaladReweighter8(feats, {"m0": sp["m0"]}, kind=kind, seed=seed, max_iter=int(o["salad_mlp_iter"]))
    rw.fit(MC.loc[sb_mc, feats].to_numpy(float), m_mc[sb_mc], SBd[feats].to_numpy(float), SBd["m"].to_numpy(float),
           max_n=int(o["salad_max_n"]))
    log.info("SALAD %s/%s: cechy %d, MC SB %d, dane SB %d, val AUC %.3f", species, kind, len(feats), sb_mc.sum(), len(SBd), rw.val_auc)
    return rw, feats


def keep_frac_for(cfg, key, species, max_rows=None):
    """Ułamek podpróbki tak, by w zakresie mas zostało ≤ max_rows kandydatów (≤0 → bez limitu). Liczy 2 kolumny (szybko)."""
    from . import samples as SA
    max_rows = int(_opts(cfg)["salad_max_rows"] if max_rows is None else max_rows)
    if max_rows <= 0:
        return 1.0
    try:
        n, _ = SA.count_species(cfg, key, species)
    except Exception as ex:  # pragma: no cover — brak źródła/kolumn: bez limitu, z ostrzeżeniem
        log.warning("s8salad: count_species(%s,%s) nieudane (%s) — bez podpróbkowania", key, species, ex)
        return 1.0
    f = 1.0 if n <= max_rows else max_rows / float(n)
    log.info("s8salad: %s/%s: %d kandydatów w zakresie → keep_frac %.4f (limit %d)", key, species, n, f, max_rows)
    return f


def _mc_all(cfg, species, keep_frac=1.0):
    """MC tła (incl_mc, truth≥min) w całym zakresie mas, po wecie D*."""
    from . import samples as SA
    sp = cfg["species"][species]
    tmin = float(_opts(cfg)["salad_mc_truth_min"])

    def sel(L):
        t = L["truth"].to_numpy(dtype=float)
        return SA.in_range(L, sp).to_numpy() & (t >= tmin) & (t < 100) & SA.dstar_veto_mask(L, sp)
    return SA.load_species(cfg, "incl_mc", species, tagx=f"saladmc{tmin:g}", select=sel, keep_frac=keep_frac)


def _data_sb(cfg, species):
    from .bkgsrc import _data_sb as _dsb
    sp = cfg["species"][species]
    return _dsb(cfg, species, sp["sideband_min"], False, "saladsb")


def training_background_salad(cfg, species, feats=None):
    """Tło treningowe = MC tła w |m−m0|<sideband_min (lub SR) z wagą SALAD (kolumna weight *= w)."""
    from . import features as F
    sp = cfg["species"][species]
    o = _opts(cfg)
    feats = list(feats or F.feature_list(cfg, cfg["training"].get("variant", "H")))
    MC = _mc_all(cfg, species, keep_frac=keep_frac_for(cfg, "incl_mc", species))
    SBd = _data_sb(cfg, species)
    rw, feats_u = fit_reweighter(cfg, species, MC, SBd, feats)
    half = sp["sideband_min"] if o["salad_region"] == "inner" else sp["signal_window"]
    B = MC[np.abs(MC["m"].to_numpy(float) - sp["m0"]) < half].reset_index(drop=True)
    w = rw.weights(B[feats_u].to_numpy(float), B["m"].to_numpy(float))
    B["w_salad"] = w
    B["weight"] = B["weight"].to_numpy(float) * w
    B["sample"] = "salad"
    neff = float(w.sum() ** 2 / (w ** 2).sum() / len(w))
    B.attrs["synthetic"] = {"model": "salad", "kind": rw.kind, "n_mc_region": int(len(B)), "n_sb_data": int(len(SBd)),
                            "val_auc_sb": rw.val_auc, "neff_over_n": neff, "features": feats_u, "region": o["salad_region"]}
    log.info("bkg_source=synthetic/salad %s: %d zdarzeń MC w |m−m0|<%.0f MeV, n_eff/N=%.2f", species, len(B), 1e3 * half, neff)
    return B


# ---------------------------------------------------------------- systematyka kształtu po cięciu
def _eff(score, pt, thr, pt_bins, wps, w=None):
    w = np.ones(len(score)) if w is None else np.asarray(w, float)
    out = {}
    tn = {wp: 0.0 for wp in wps}; tk = {wp: 0.0 for wp in wps}
    for lo, hi in zip(pt_bins[:-1], pt_bins[1:]):
        s = (pt >= lo) & (pt < hi)
        key = f"{lo}-{hi}"
        out[key] = {}
        for wp in wps:
            t = thr[key][f"{wp:.2f}"]
            n, k = w[s].sum(), w[s & (score >= t)].sum()
            out[key][f"{int(wp*100)}"] = float(k / n) if n > 0 else float("nan")
            tn[wp] += n; tk[wp] += k
    out["all"] = {f"{int(wp*100)}": float(tk[wp] / tn[wp]) if tn[wp] > 0 else float("nan") for wp in wps}
    return out


def shape_syst(cfg, species, variant=None, kind=None, out=None):
    from . import model as M
    from . import samples as SA
    from .bkggen.feta import c2st
    from .paths import mdir
    variant = variant or cfg["training"]["variant"]
    sp = cfg["species"][species]
    tr = cfg["training"]
    pt_bins, wps = tr["pt_bins"], tr["working_points"]
    ens = M.Ensemble(mdir(cfg, species, variant))
    meta = ens.meta
    thr = meta["thresholds"]
    feats = list(ens.feats)
    fM, fD = keep_frac_for(cfg, "incl_mc", species), keep_frac_for(cfg, "data_none", species)
    MC = _mc_all(cfg, species, keep_frac=fM)
    D = SA.load_species(cfg, "data_none", species, tagx="saladall",
                        select=lambda L: SA.in_range(L, sp).to_numpy() & SA.dstar_veto_mask(L, sp), keep_frac=fD)
    mD, mM = D["m"].to_numpy(float), MC["m"].to_numpy(float)
    dD, dM = np.abs(mD - sp["m0"]), np.abs(mM - sp["m0"])
    SBd = D[dD > sp["sideband_min"]].reset_index(drop=True)
    rw, feats_u = fit_reweighter(cfg, species, MC, SBd, feats, kind=kind)
    wM = rw.weights(MC[feats_u].to_numpy(float), mM)
    sM, _ = ens.score(MC)
    sD, _ = ens.score(D)
    ptM, ptD = MC["pt"].to_numpy(float), D["pt"].to_numpy(float)
    sr, band, sbm = dM < sp["signal_window"], (dM >= sp["signal_window"]) & (dM < sp["sideband_min"]), dM > sp["sideband_min"]
    bandD = (dD >= sp["signal_window"]) & (dD < sp["sideband_min"])
    res = {
        "species": species, "variant": variant, "kind": rw.kind, "source": "S8: SALAD (MC incl_mc truth≥8 → dane SB) vs SB danych vs MC surowe",
        "definition": {"effB_salad_SR": "ε_B @WP na MC tła w SR ważonym SALAD (w(x,m) uczone w SB)",
                       "effB_mcraw_SR": "ε_B @WP na MC tła w SR bez wag",
                       "effB_data_SB": "ε_B @WP na danych w SB (nominalny proxy charmbkg)",
                       "effB_salad_SB": "ε_B @WP na MC SB ważonym (zamknięcie: ≈ effB_data_SB)",
                       "rel_salad_vs_SB": "effB_salad_SR / effB_data_SB − 1 → systematyka kształtu tła po cięciu",
                       "c2st": "AUC klasyfikatora dane vs model w paśmie inner_band (0,5 = nieodróżnialne)"},
        "val_auc_sb": rw.val_auc, "n_mc": int(len(MC)), "n_data": int(len(D)), "features": feats_u,
        "keep_frac": {"incl_mc": fM, "data_none": fD},
        "neff_over_n_SR": float(wM[sr].sum() ** 2 / (wM[sr] ** 2).sum() / max(sr.sum(), 1)),
        "effB_salad_SR": _eff(sM[sr], ptM[sr], thr, pt_bins, wps, wM[sr]),
        "effB_mcraw_SR": _eff(sM[sr], ptM[sr], thr, pt_bins, wps),
        "effB_data_SB": _eff(sD[dD > sp["sideband_min"]], ptD[dD > sp["sideband_min"]], thr, pt_bins, wps),
        "effB_salad_SB": _eff(sM[sbm], ptM[sbm], thr, pt_bins, wps, wM[sbm]),
        "effB_mcraw_SB": _eff(sM[sbm], ptM[sbm], thr, pt_bins, wps),
    }
    rel = {}
    for key, row in res["effB_salad_SR"].items():
        rel[key] = {wp: (row[wp] / res["effB_data_SB"][key][wp] - 1.0) if res["effB_data_SB"][key][wp] else None for wp in row}
    res["rel_salad_vs_SB"] = rel
    try:
        XD, XM = D.loc[bandD, feats_u].to_numpy(float), MC.loc[band, feats_u].to_numpy(float)
        res["c2st"] = {"salad_band": c2st(XD, XM, None, wM[band]), "mcraw_band": c2st(XD, XM),
                       "sbcopy_band": c2st(XD, SBd[feats_u].to_numpy(float))}
    except Exception as ex:  # pragma: no cover
        res["c2st"] = {"error": str(ex)}
    od = out or os.path.join(cfg["paths"]["out_dir"], "syst")
    os.makedirs(od, exist_ok=True)
    p = os.path.join(od, f"bkg_shape_salad_{species}.json")
    with open(p, "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log.info("s8salad %s: zapisano %s (rel all: %s)", species, p, rel.get("all"))
    return res


def main(argv=None):
    from . import config as C
    ap = argparse.ArgumentParser(prog="charmbkg s8salad", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("s8salad", nargs="?")
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--species", default="Dplus", choices=["D0", "Dplus"])
    ap.add_argument("--variant", default=None)
    ap.add_argument("--kind", default=None, choices=["sk_mlp", "hgb", "hgb_nom", "mlp"])
    ap.add_argument("--out", default=None)
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = C.load(a.config)
    r = shape_syst(cfg, a.species, a.variant, a.kind, a.out)
    print(json.dumps({"rel_salad_vs_SB": r["rel_salad_vs_SB"], "c2st": r["c2st"], "val_auc_sb": r["val_auc_sb"]}, indent=1))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
