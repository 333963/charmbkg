#!/usr/bin/env python3
"""AFP-MIX v2 na danych (projekt B, WSTĘPNE): tło przypadkowego tagu AFP z mieszania protonów — parquet z afp_extract.py / afp_merge.py
(<próbka>_events.parquet, <próbka>_cands.parquet, <próbka>.meta.json) → estymator A (warstwy LB × nSelBL), drugie estymatory
B(i), B(ii), bramki v2, rozbicie błędu → RAPORT_DANE_v2.md.

Metoda i bramki: charmbkg.afp_mix v2 (docs/afp_mix_method.md). Tutaj tylko sterowanie w skali danych:
  * selekcja zdarzeń: GRL ∧ tryger nie-AFP ∧ nPileupVertices == 0 (ZAKRES WAŻNOŚCI: tylko trygery nie-AFP — W2),
  * klasa trygera (wariant puli) i próbka B(i) (łańcuchy quasi-minbias L1MBTS) z kolumn trigPassed_*,
  * afp_mix.run_procedure (pełna statystyka, także w wariantach puli) + tabele per gatunek / pT / okno masy, μ, LB,
  * raportowane są WYŁĄCZNIE frakcje tła przypadkowego i liczności tagów; żadnych uzysków sygnału / przekrojów (bez fitów masy).

Użycie:
  run_data.py [--in DIR] [--out DIR] [--samples a,b]            # dane (domyślnie runs/lxplus → runs/lxplus/out_v2)
  run_data.py --make-toy DIR [--toy-scenario nominal]           # toy w formacie afp_ext (+ <próbka>_truth.parquet)
  run_data.py --in DIR --out DIR --tag-def inwindow             # toy: tag = dokładnie 1 proton W OKNIE (dane: n_side == 1)
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("CHARMBKG", os.path.dirname(os.path.dirname(HERE))))   # korzeń repo charmbkg
from charmbkg import afp_mix as M  # noqa: E402

AFP_CHAINS = ("HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C")
# klasa trygera (wyłączna, wg priorytetu): 3 = L1J12, 2 = j20, 1 = sptrk_pt4_L1MBTS_2, 0 = quasi-minbias L1MBTS (reszta)
TRIG_PRIORITY = (("HLT_noalg_L1J12", 3), ("HLT_j20_L1MBTS_2", 2), ("HLT_mb_sptrk_pt4_L1MBTS_2", 1))
MB_CHAINS = ("HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2")
RANDOM_CHAIN = "HLT_noalg_L1RD0_FILLED"
SPECIES = {   # kod Dmeson_species, zmienna okna, pik, sidebandy [GeV]
    "Dplus": dict(code=0, var="Dmeson_mass", peak=(1.84, 1.90), sb=((1.75, 1.80), (1.94, 2.00))),
    "D0": dict(code=3, var="Dmeson_mass", peak=(1.835, 1.895), sb=((1.75, 1.80), (1.93, 2.00))),
    "Dstar": dict(code=1, var="Dmeson_delta_mass", peak=(0.1440, 0.1470), sb=((0.150, 0.160),)),
}
PT_EDGES = (3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 150.0)
REGIONS = ("AorC", "A", "C")


# ---------------------------------------------------------------------------------------------------------------- wejście
def load_sample(d, s):
    meta = json.load(open(os.path.join(d, f"{s}.meta.json")))
    ev = pd.read_parquet(os.path.join(d, f"{s}_events.parquet"))
    ca = pd.read_parquet(os.path.join(d, f"{s}_cands.parquet"))
    run = np.uint32(meta["run_number"])
    ev["run_number"] = run
    ca["run_number"] = run
    return ev, ca, meta


def trigger_classes(ev):
    """(nonafp, afp, tabela łańcuchów, kolumny). nonafp = przeszedł ≥ 1 zapisany łańcuch bez L1AFP; afp = przeszedł łańcuch L1AFP."""
    cols = [c for c in ev.columns if c.startswith("trigPassed_")]
    n = len(ev)
    nonafp, afp, tab = np.zeros(n, bool), np.zeros(n, bool), {}
    for c in cols:
        p = ev[c].to_numpy() == 1
        tab[c[len("trigPassed_"):]] = float(p.mean()) if n else 0.0
        if c[len("trigPassed_"):] in AFP_CHAINS or "L1AFP" in c:
            afp |= p
        else:
            nonafp |= p
    return nonafp, afp, tab, cols


def trig_class_and_mb(ev):
    """Klasa trygera (int8, wyłączna wg TRIG_PRIORITY) i maska „przeszedł łańcuch quasi-minbias L1MBTS”."""
    n = len(ev)
    cls = np.zeros(n, np.int8)
    for ch, code in reversed(TRIG_PRIORITY):
        c = "trigPassed_" + ch
        if c in ev:
            cls[ev[c].to_numpy() == 1] = code
    mb = np.zeros(n, bool)
    for ch in MB_CHAINS:
        c = "trigPassed_" + ch
        if c in ev:
            mb |= ev[c].to_numpy() == 1
    return cls, mb


def standardize_events(ev):
    E = M.standardize(ev)                                   # ← afp_mix.COLMAP (nazwy h5 v2.0)
    for s in M.SIDES:
        E[f"np_{s}"] = E[f"np_{s}"].to_numpy().astype(np.int32)
    return E


def qedges(x, nb):
    e = np.unique(np.quantile(x, np.linspace(0, 1, nb + 1)))
    if len(e) < 3:
        e = np.array([x.min(), x.max() + 1e-9]) if len(e) < 2 else e
    e = e.astype(float)
    e[-1] = np.nextafter(e[-1], np.inf) + 1e-12
    return e


def est_small(mx, boot, tw, rel_syst):
    e = M.estimate(mx, boot, tw=tw)
    out = {k: e[k] for k in ("n_tag", "n_acc", "n_acc_err_stat", "f_acc", "f_acc_err_stat", "eps_xical_acc")}
    out["f_acc_err_syst"] = e["f_acc"] * rel_syst if np.isfinite(e["f_acc"]) else float("nan")
    return out


# ---------------------------------------------------------------------------------------------------------------- analiza próbki
def analyse_sample(s, ev, ca, meta, a, log):
    cfg = M.MixConfig(n_min_pool=a.n_min_pool, tag_def=a.tag_def)
    R = {"sample": s, "run_number": int(meta["run_number"]), "config": {"xi_window": [cfg.xi_lo, cfg.xi_hi], "delta": cfg.delta,
         "veto_thr": cfg.veto_thr, "strata": cfg.strata, "act_edges": list(cfg.act_edges), "n_min_pool": cfg.n_min_pool, "n_boot": a.n_boot,
         "tag_def": a.tag_def, "cr_margin": cfg.cr_margin, "cr_xmax": cfg.cr_x, "trigger_mode": a.trigger_mode, "dstar_lxy_min": a.dstar_lxy_min}}
    # ------------------------------------------------------------ selekcja zdarzeń
    n0 = len(ev)
    grl = ev["passGRL_MinBiasAFP"].to_numpy() == 1 if "passGRL_MinBiasAFP" in ev else np.ones(n0, bool)
    nonafp, afp, trig_tab, tcols = trigger_classes(ev)
    if not tcols:
        nonafp = np.ones(n0, bool)
    afp_only = afp & ~nonafp
    no_trig = ~afp & ~nonafp
    if a.trigger_mode == "nonafp":
        tsel = nonafp
    elif a.trigger_mode.startswith("chain:"):
        tsel = ev["trigPassed_" + a.trigger_mode[6:]].to_numpy() == 1
    else:
        raise SystemExit(f"nieznany --trigger-mode {a.trigger_mode}")
    npu0 = ev["nPileupVertices"].to_numpy() == 0
    base = grl & tsel
    E = standardize_events(ev)
    tgl = E["tag_A"].to_numpy().astype(bool) | E["tag_C"].to_numpy().astype(bool)
    tcl, mbm = trig_class_and_mb(ev)
    rd0 = (ev["trigPassed_" + RANDOM_CHAIN].to_numpy() == 1) if "trigPassed_" + RANDOM_CHAIN in ev else np.zeros(n0, bool)
    R["inputs"] = {"n_events": n0, "n_grl": int(grl.sum()), "n_afp_only": int(afp_only.sum()), "n_afp_only_grl": int((afp_only & grl).sum()),
                   "n_afp_only_grl_tagged": int((afp_only & grl & tgl).sum()), "n_sel_tagged": int((base & npu0 & tgl).sum()),
                   "n_no_stored_trigger": int(no_trig.sum()), "n_trig_sel_grl": int(base.sum()), "n_npu0_of_base": int((base & npu0).sum()),
                   "lb_range": [int(E["lb"].min()), int(E["lb"].max())], "n_lb": int(E["lb"].nunique()),
                   "mu_quantiles_5_50_95": [float(x) for x in np.quantile(E["mu"].to_numpy()[base], [0.05, 0.5, 0.95])] if base.any() else [],
                   "trig_pass_frac": trig_tab, "n_random_trigger_RD0": int(rd0.sum()),
                   "afp_only_tag_rate": float(tgl[afp_only].mean()) if afp_only.any() else None, "base_tag_rate": float(tgl[base].mean()) if base.any() else None}
    # ------------------------------------------------------------ kandydaci → zdarzenia
    spec = ca["Dmeson_species"].to_numpy()
    lxy = ca["Dmeson_Lxy"].to_numpy()
    crow = pd.Index(ev["event_number"].to_numpy()).get_indexer(ca["event_number"].to_numpy())
    if (crow < 0).any():
        raise SystemExit(f"{s}: {int((crow < 0).sum())} kandydatów bez zdarzenia w tabeli zdarzeń")
    csel = np.zeros(len(ca), bool)
    cwin = np.zeros(len(ca), np.int8)      # 1 pik, 2 sideband
    for name, sp in SPECIES.items():
        ms = spec == sp["code"]
        if name == "Dstar" and a.dstar_lxy_min > 0:
            ms &= lxy >= a.dstar_lxy_min
        csel |= ms
        v = ca[sp["var"]].to_numpy()
        cwin[ms & (v >= sp["peak"][0]) & (v <= sp["peak"][1])] = 1
        for lo, hi in sp["sb"]:
            cwin[ms & (v >= lo) & (v <= hi)] = 2
    has = np.zeros(n0, bool)
    has[crow[csel]] = True
    sel = base & npu0
    E["hasD"] = has
    E["trig"] = tcl
    mb_frac = float(mbm[sel].mean()) if sel.any() else 0.0
    # próbka B(i): łańcuchy quasi-minbias L1MBTS; gdy są bez preskali (≥ 50 % zdarzeń) — tylko zdarzenia bez D (rozłączne z biorcami)
    b1 = mbm & (~has if mb_frac >= 0.5 else True)
    E["b1"] = b1
    t_rows = np.flatnonzero(sel & has)
    Dn = E.loc[sel].reset_index(drop=True)
    T = E.iloc[t_rows].reset_index(drop=True)
    del E
    row2t = np.full(n0, -1, np.int64)
    row2t[t_rows] = np.arange(len(t_rows))
    ct = row2t[crow]
    nT = len(T)
    tgD, tgT = M.side_tags(Dn, cfg), M.side_tags(T, cfg)
    R["inputs"].update({"n_donors": len(Dn), "n_targets": nT, "n_cands": int(len(ca)), "n_cands_in_targets": int(((ct >= 0) & csel).sum()),
                        "mb_chain_frac_of_selected": mb_frac, "b1_definition": "łańcuchy L1MBTS quasi-minbias" + (" ∧ bez D" if mb_frac >= 0.5 else " (z D i bez D)"),
                        "n_b1": int((b1 & sel).sum()), "trig_class_frac_targets": {int(k): float(v) for k, v in zip(*np.unique(T["trig"], return_counts=True))},
                        "trig_class_frac_donors": {int(k): float(v) for k, v in zip(*np.unique(Dn["trig"], return_counts=True))}})
    R["tag_counts"] = {"donors": {"A": int(tgD["A"][0].sum()), "C": int(tgD["C"][0].sum()), "AandC": int((tgD["A"][0] & tgD["C"][0]).sum()), "n": len(Dn)},
                       "targets": {"A": int(tgT["A"][0].sum()), "C": int(tgT["C"][0].sum()), "AandC": int((tgT["A"][0] & tgT["C"][0]).sum()), "n": nT}}
    log(f"{s}: zdarzeń {n0}, GRL∧tryger {int(base.sum())}, dawcy {len(Dn)}, biorcy {nT}, próbka B(i) {R['inputs']['n_b1']}")
    if nT < 100 or len(Dn) < 1000:
        R["error"] = "za mało biorców/dawców"
        return R

    # ------------------------------------------------------------ procedura v2
    mu_t, lb_t = T["mu"].to_numpy().astype(float), T["lb"].to_numpy().astype(float)
    gc = M.GateConfig(n_boot=a.n_boot, n_pseudo=a.n_pseudo, mu_edges=tuple(qedges(mu_t, 7)))
    t0 = time.time()
    rep = M.run_procedure(None, Dn, cfg, gc, "AorC", targets=T, n_afp_only=0, rng=a.seed, log=lambda m: log(f"{s}: {m}"), keep=True)
    obj = rep.pop("_obj")
    mx, boot = obj["mx"], obj["boot"]
    log(f"{s}: procedura v2 ({time.time() - t0:.0f} s), werdykt {rep['verdict']}")
    g0 = rep["gates"][0]
    g0.update(n_afp_only_excluded=int((afp_only & grl).sum()), n_afp_only_excluded_tagged=R["inputs"]["n_afp_only_grl_tagged"],
              note="zdarzenia wyłącznie z trygera L1AFP wykluczone u biorców i dawców — dla nich wynik NIE obowiązuje (W2)")
    R.update({k: rep[k] for k in ("estimate", "estimates_by_region", "variants", "info_variants", "profile", "method_Bi", "method_Bii",
                                  "syst_rel", "gates", "verdict", "verdict_failed", "verdict_na")})
    sy = rep["syst_rel"]
    rel_syst = math.sqrt(sum(sy[k] ** 2 for k in ("pool", "extrapolation", "half_diff_methods") if np.isfinite(sy[k])))
    R["syst_rel"]["syst_no_stat"] = rel_syst

    # ------------------------------------------------------------ f_acc(strona, μ) i f_acc(strona, LB)
    R["by_mu"], R["by_lb"] = [], []
    for x, nb, key in ((mu_t, a.n_mu_bins, "by_mu"), (lb_t, a.n_lb_bins, "by_lb")):
        ed = qedges(x, nb)
        b = np.digitize(x, ed) - 1
        for i in range(len(ed) - 1):
            mk = (b == i).astype(float)
            if mk.sum() == 0:
                continue
            row = {"lo": float(ed[i]), "hi": float(ed[i + 1]), "mean": float(x[b == i].mean())}
            for reg in REGIONS:
                row[reg] = est_small(mx[reg], boot, mk, rel_syst)
            R[key].append(row)

    # ------------------------------------------------------------ gatunki × pT × okno (wagi = liczba kandydatów w zdarzeniu)
    pt = ca["Dmeson_pt"].to_numpy()
    okc = (ct >= 0) & csel
    R["species"] = {}
    for name, sp in SPECIES.items():
        ms = okc & (spec == sp["code"])
        if name == "Dstar" and a.dstar_lxy_min > 0:
            ms &= lxy >= a.dstar_lxy_min
        rows = []
        bins = [(PT_EDGES[0], PT_EDGES[-1])] + list(zip(PT_EDGES[:-1], PT_EDGES[1:]))
        for lo, hi in bins:
            mb_ = ms & (pt >= lo) & (pt < hi)
            row = {"pt_lo": lo, "pt_hi": hi}
            for wn, wc in (("peak", 1), ("sb", 2)):
                w = np.bincount(ct[mb_ & (cwin == wc)], minlength=nT).astype(float)
                row[wn] = est_small(mx["AorC"], boot, w, rel_syst)
                for reg in ("A", "C"):
                    row[f"{wn}_{reg}"] = est_small(mx[reg], boot, w, rel_syst)
            pk, sb = row["peak"], row["sb"]
            if pk["n_tag"] > 0 and sb["n_tag"] > 0 and np.isfinite(pk["f_acc_err_stat"]) and np.isfinite(sb["f_acc_err_stat"]) and sb["f_acc"] > 0 and pk["f_acc"] > 0:
                r = pk["f_acc"] / sb["f_acc"]
                row["peak_over_sb"] = r
                row["peak_over_sb_err"] = r * math.hypot(pk["f_acc_err_stat"] / pk["f_acc"], sb["f_acc_err_stat"] / sb["f_acc"])
            rows.append(row)
        R["species"][name] = {"windows": {"var": sp["var"], "peak": sp["peak"], "sb": sp["sb"]}, "rows": rows}

    # ------------------------------------------------------------ I6 (informacyjna): η* kandydatów w zdarzeniach z tagiem
    obs_t = mx["AorC"].own
    cobs = okc.copy()
    cobs[cobs] = obs_t[ct[cobs]]
    if cobs.sum() >= 50:
        Cs = pd.DataFrame({"run": np.uint32(meta["run_number"]), "event": ca["event_number"].to_numpy()[cobs], "eta": ca["Dmeson_eta"].to_numpy()[cobs]})
        g = M.info_i6_eta(Cs, T, cfg, gc, R["estimate"]["f_acc"])
        g["selection"] = "wszyscy kandydaci slim"
        R["gates"].append(g)

    # ------------------------------------------------------------ prawda toy (jeśli jest)
    tp = os.path.join(a.indir, f"{s}_truth.parquet")
    if os.path.exists(tp):
        tr = pd.read_parquet(tp).iloc[t_rows].reset_index(drop=True)
        n_true = int(tr["acc_region"].sum())
        inc = R["estimate"]
        R["toy_truth"] = {"n_acc": n_true, "eps_xical_acc": float(n_true / max(tr["acc_tag_nocut"].sum(), 1)),
                          "fbar": float(tr["lost"].sum() / max(tr["own_region"].sum(), 1)), "f_acc": float(n_true / max(inc["n_tag"], 1)),
                          "z_n_acc": float((inc["n_acc"] - n_true) / math.sqrt(max(n_true, 1) + (inc["n_acc"] * sy["total_rel"]) ** 2))}
    if a.save_p:
        out = pd.DataFrame({"event_number": ev["event_number"].to_numpy()[t_rows]})
        for reg in REGIONS:
            out[f"p_tag_{reg}"] = mx[reg].p_tag.astype(np.float32)
            out[f"p_acc_{reg}"] = mx[reg].p_acc().astype(np.float32)
            out[f"obs_{reg}"] = mx[reg].own
        out.to_parquet(os.path.join(a.out, f"{s}_p_events.parquet"), index=False)
    return R


# ---------------------------------------------------------------------------------------------------------------- raport
def f(x, n=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:.{n}f}"


def pm_(x, e, n=4):
    return f"{f(x, n)} ± {f(e, n)}"


def pct(x, n=1):
    return "—" if x is None or not np.isfinite(x) else f"{100 * x:.{n}f} %"


GATE_KEYS = ("coverage", "frac_no_pool", "median_pool", "n_afp_only_excluded", "lr", "ndf", "n_bins", "p", "n_pool_tags", "n_obs", "n_pred", "n_pred_err",
             "ratio", "ratio_err", "z", "n_acc_A", "n_acc_Bi", "n_tagged_donors", "slope_data", "slope_data_err", "slope_mix", "z_slope",
             "A_obs", "A_err", "f_acc_eta", "f_acc_eta_err", "f_acc_mix", "n_cand")


def gate_rows(R):
    out = []
    for g in R["gates"]:
        st = g["status"] + (" (informacyjna)" if g["status"] == "INFO" else "")
        kv = ", ".join(f"{k}={g[k]:.4g}" for k in GATE_KEYS if k in g and isinstance(g[k], (int, float, np.integer, np.floating)))
        if "per_side" in g and g["gate"] != "G1_pool_homogeneity":
            kv += "; strony: " + ", ".join(f"{sd}: {v.get('n_obs', 0):.0f}/{v.get('n_pred', 0):.1f}" for sd, v in g["per_side"].items())
        if "note" in g:
            kv += f"; {g['note']}"
        out.append(f"| {g['gate']} | **{st}** | {kv} | {g['tests']} | {g['not_tests']} |")
    return out


def load_v1(a):
    p = a.v1_result or os.path.join(a.indir, "out", "result.json")
    try:
        return {r["sample"]: r for r in json.load(open(p))["results"]}
    except Exception:
        return {}


def write_report(RR, a, path, dt):
    toy = any("toy_truth" in R for R in RR)
    V1 = load_v1(a)
    L = [f"# AFP-MIX v2 na danych — WSTĘPNE{' (TOY — test formatu, nie dane)' if toy else ''}", "",
         f"**STATUS: WSTĘPNE.** Wygenerowano {time.strftime('%Y-%m-%d %H:%M')} przez `analysis/afp_mix/run_data.py` v2 (czas {dt:.0f} s). Wejście: `{a.indir}`.",
         "Metoda: `charmbkg/afp_mix.py` v2, `docs/afp_mix_method.md` (charmbkg). "
         "Kandydaci **bez score** (luźna selekcja v2.0) — liczby pod pikiem to wszyscy kandydaci w oknie (sygnał + kombinatoryka), bez fitu. "
         "Raport zawiera wyłącznie frakcje tła przypadkowego tagu i liczności tagów; uzysków sygnału dyfrakcyjnego ani przekrojów nie liczono.", "",
         "## 0. Zakres ważności", "",
         "**Wynik dotyczy wyłącznie zdarzeń wybranych trygerami nie-AFP** (GRL ∧ ≥ 1 zapisany łańcuch bez L1AFP ∧ nPileupVertices == 0). "
         "Zdarzenia wybrane wyłącznie przez łańcuchy L1AFP (tryger wymaga protonu) są wykluczone u biorców i dawców; f_acc z tej metody ich NIE opisuje "
         "(inna aktywność, inne ξ_cal, tag wymuszony trygerem). Liczby per próbka poniżej („Wejścia”). Próbek nie łączono.", "",
         "## 1. Wynik: f_acc (region A∨C, po cięciu ξ_cal) z rozbiciem błędu", "",
         "f_acc = N_acc/N_tag, estymator A (mieszanie w warstwach LB × nSelBL × „z D/bez D”, pule jednostronne ξ_cal > 0,14). Składniki błędu względnego N_acc: stat = bootstrap otagowanych dawców; "
         "pula = największa istotna różnica wariantów puli (|Δ| − σ_Δ, pełna statystyka, wspólny bootstrap); ekstrapolacja = √((κ−1)² + σ_κ²) z liniowego modelu R(ln ξ_cal) w puli; "
         "½Δ metod = połowa największej różnicy A–B(i), A–B(ii).", "",
         "| próbka | N_tag | N_acc (A) | **f_acc** | stat | pula | ekstrapolacja | ½Δ metod | razem (wzgl.) | **f_acc ± razem** | v1 f_acc | werdykt bramek v2 |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for R in RR:
        if "estimate" not in R:
            L.append(f"| {R['sample']} | {R.get('error', '?')} |")
            continue
        e, sy = R["estimate"], R["syst_rel"]
        v1 = V1.get(R["sample"], {}).get("inclusive", {}).get("AorC", {}).get("f_acc")
        vd = f"**{R['verdict']}**" + (f" (FAIL: {', '.join(R['verdict_failed'])})" if R["verdict_failed"] else "") + (f"; N/A: {', '.join(R['verdict_na'])}" if R["verdict_na"] else "")
        L.append(f"| {R['sample']} | {e['n_tag']:.0f} | {pm_(e['n_acc'], e['n_acc_err_stat'], 1)} | **{f(e['f_acc'])}** | {pct(sy['stat'])} | {pct(sy['pool'])} | {pct(sy['extrapolation'])} | "
                 f"{pct(sy['half_diff_methods'])} | {pct(sy['total_rel'])} | **{f(e['f_acc'])} ± {f(e['f_acc'] * sy['total_rel'])}** | {f(v1)} | {vd} |")
    L += ["", "Uwaga: składniki „ekstrapolacja” i „½Δ metod (B(ii))” dotyczą tego samego pytania (zależność częstości od ξ_cal) — suma kwadratowa jest konserwatywna.", ""]
    for R in RR:
        s = R["sample"]
        I = R["inputs"]
        L += [f"## {s} (run {R['run_number']})", ""]
        if "error" in R:
            L += [f"BŁĄD: {R['error']}", ""]
            continue
        e, sy = R["estimate"], R["syst_rel"]
        L += ["### Wejścia, selekcja, zakres", "",
              f"- zdarzeń {I['n_events']:,}; GRL {I['n_grl']:,}; LB {I['lb_range'][0]}–{I['lb_range'][1]} ({I['n_lb']} bloków); μ (5/50/95 %) = "
              + "/".join(f"{x:.4f}" for x in I["mu_quantiles_5_50_95"]),
              f"- **zakres (W2):** zdarzenia wyłącznie z łańcuchów L1AFP po GRL: {I['n_afp_only_grl']:,}, w tym z tagiem {I['n_afp_only_grl_tagged']:,} — **wykluczone**; "
              f"zdarzenia z tagiem w wybranej próbce (nie-AFP, npu0): {I['n_sel_tagged']:,} → metoda opisuje "
              f"{pct(I['n_sel_tagged'] / max(I['n_sel_tagged'] + I['n_afp_only_grl_tagged'], 1))} wszystkich tagów po GRL; bez żadnego zapisanego łańcucha: {I['n_no_stored_trigger']:,}",
              f"- po GRL ∧ tryger: {I['n_trig_sel_grl']:,}; z nPileupVertices == 0: {I['n_npu0_of_base']:,}; dawcy {I['n_donors']:,}; biorcy (≥ 1 kandydat slim) {I['n_targets']:,}",
              f"- tryger losowy {RANDOM_CHAIN}: **{I['n_random_trigger_RD0']:,} zdarzeń** → absolutne P(tag) z trygera losowego niedostępne; "
              f"B(i) używa: {I['b1_definition']} ({I['n_b1']:,} zdarzeń; łańcuchy L1MBTS to {pct(I['mb_chain_frac_of_selected'])} wybranych)",
              "- klasy trygera (3 = L1J12, 2 = j20, 1 = sptrk_pt4_L1MBTS_2, 0 = pozostałe L1MBTS): biorcy "
              + ", ".join(f"{k}: {v / max(I['n_targets'], 1):.3f}" for k, v in I["trig_class_frac_targets"].items())
              + "; dawcy " + ", ".join(f"{k}: {v / max(I['n_donors'], 1):.3f}" for k, v in I["trig_class_frac_donors"].items()),
              "- liczności tagów (przed cięciem ξ_cal): dawcy A {A:,} / C {C:,} / A∧C {AandC:,} z {n:,}".format(**R["tag_counts"]["donors"])
              + "; biorcy A {A:,} / C {C:,} / A∧C {AandC:,} z {n:,}".format(**R["tag_counts"]["targets"]), ""]
        b1, b2 = R["method_Bi"], R["method_Bii"]
        v1r = V1.get(s, {})
        L += ["### Estymatory (N_acc, region A∨C)", "", "| estymator | N_acc | N/N_A | uwagi |", "|---|---|---|---|",
              f"| **A** — mieszanie, warstwy LB × nSelBL × „z D”, pule jednostronne ξ_cal > {R['config']['veto_thr']:.2f} (nominał) | {pm_(e['n_acc'], e['n_acc_err_stat'], 1)} | 1 | otagowanych dawców {e['n_tagged_donors']}; ε_ξcal = {f(e['eps_xical_acc'], 3)}; f̄ = {f(e['fbar'])} |",
              f"| informacyjnie: warstwy LB × nSelBL, dawcy z D i bez D | {f(R['info_variants']['donors_D_and_noD']['n_acc'], 1)} | {f(R['info_variants']['donors_D_and_noD']['ratio'], 3)} ± {f(R['info_variants']['donors_D_and_noD']['delta_err'], 3)} | nie w systematyce — efekt dawców bez D (W1; zob. I2) |",
              f"| informacyjnie: tylko warstwy LB, wszyscy dawcy (≈ v1, pule jednostronne) | {f(R['info_variants']['v1_like_LB_only']['n_acc'], 1)} | {f(R['info_variants']['v1_like_LB_only']['ratio'], 3)} ± {f(R['info_variants']['v1_like_LB_only']['delta_err'], 3)} | nie w systematyce — wielkość efektu warstwowania |"]
        if v1r:
            vi = v1r["inclusive"]["AorC"]
            L.append(f"| v1 (raport `out/RAPORT_DANE.md`) | {pm_(vi['n_acc'], vi['n_acc_err_stat'], 1)} | {f(vi['n_acc'] / e['n_acc'], 3)} | f_acc v1 = {f(vi['f_acc'])}; v1 syst. dekl. {pct(v1r.get('syst_rel', {}).get('total_rel'))} |")
        if b1.get("n_acc") is not None:
            L.append(f"| **B(i)** — druga pula: {I['b1_definition']}, weto ξ_cal po obu stronach, nSelBL ≥ {b1.get('nsel_min')}, absolutne P(tag) w {len(b1.get('mu_edges', [])) + 1} przedziałach μ | {pm_(b1['n_acc'], b1['n_acc_err'], 1)} | "
                     f"{f(1 + b1['diff_rel'], 3)} ± {f(b1['diff_rel_err'], 3)} | otagowanych dawców {b1['n_tagged_donors']} z {b1['n_pool']:,}{'' if b1['available'] else ' — ZA MAŁO, nieużyte'} |")
        else:
            L.append("| **B(i)** | — | — | próbka niedostępna |")
        if b2.get("n_acc") is not None:
            L.append(f"| **B(ii)** — A × N_obs/N_pred w obszarze kontrolnym (ξ_cal ≤ {R['config']['cr_xmax']:.2f}, Δ > {0.02 + R['config']['cr_margin']:.3f}) | {pm_(b2['n_acc'], b2['n_acc_err'], 1)} | "
                     f"{f(b2['ratio'], 3)} ± {f(b2['diff_rel_err'], 3)} | CR: obs {b2['n_obs_cr']:.0f}, pred {b2['n_pred_cr']:.1f}{'' if b2['available'] else ' — ZA MAŁO, nieużyte'} |")
        L += ["", f"Rozbicie błędu względnego N_acc: stat {pct(sy['stat'])}, pula {pct(sy['pool'])}, ekstrapolacja {pct(sy['extrapolation'])} "
              f"(κ = {f(sy['kappa'], 3)} ± {f(sy['kappa_err'], 3)}), ½Δ metod {pct(sy['half_diff_methods'])} "
              f"({', '.join(f'{k}: {pct(v)}' for k, v in sy['half_diff_by_method'].items())}) → **razem {pct(sy['total_rel'])}**.", "",
              "Warianty puli (pełna statystyka; błąd = błąd różnicy względem nominału ze wspólnego bootstrapu):", "",
              "| wariant | N_acc / nominał | istotna część |", "|---|---|---|"]
        for k, v in R["variants"].items():
            L.append(f"| {k} | {f(v['ratio'], 3)} ± {f(v['delta_err'], 3)} | {pct(v['significant'])} |")
        L += ["", "Profil częstości tagu w puli po warstwowaniu (R = obs/pred; G1 i ekstrapolacja):", "",
              "| strona | przedziały ξ_cal | tagi k | pred | nachylenie dR/d ln ξ_cal | ⟨ln ξ_cal⟩ pula → użycie | κ |", "|---|---|---|---|---|---|---|"]
        for sd, p in R["profile"].items():
            L.append(f"| {sd} | {', '.join(f'{x:.2g}' for x in p['edges'])} | {', '.join(f'{x:.0f}' for x in p['k'])} | {', '.join(f'{x:.1f}' for x in p['pred'])} | "
                     f"{pm_(p.get('slope'), p.get('slope_err'), 3)} | {f(p.get('u_pool'), 2)} → {f(p.get('u_use'), 2)} | {pm_(p.get('kappa'), p.get('kappa_err'), 3)} |")
        g7 = next(g for g in R["gates"] if g["gate"] == "G7_methods_A_vs_Bii")
        L += ["", "Obszar kontrolny — skan dolnej granicy Δ − δ (spadek R z marginesem = przeciek własnych protonów przez rozdzielczość ξ_cal):", "",
              "| margines | N_obs | N_pred | R |", "|---|---|---|---|"]
        for c in g7.get("margin_scan", []):
            L.append(f"| {c['margin']:.3f} | {c['n_obs']:.0f} | {c['n_pred']:.1f} | {pm_(c['ratio'], c['ratio_err'], 3)} |")
        L += ["", "### Bramki v2", "", "| bramka | status | kluczowe | co testuje | czego NIE testuje |", "|---|---|---|---|---|"] + gate_rows(R) + [""]
        L += ["### f_acc inkluzywnie per region", "", "| region | N_tag | N_acc (A) | f_acc ± stat | ε_ξcal(przyp.) | f̄ | v1 f_acc |", "|---|---|---|---|---|---|---|"]
        for reg in REGIONS:
            x = R["estimates_by_region"][reg]
            v1 = v1r.get("inclusive", {}).get(reg, {}).get("f_acc") if v1r else None
            L.append(f"| {reg} | {x['n_tag']:.0f} | {pm_(x['n_acc'], x['n_acc_err_stat'], 1)} | {pm_(x['f_acc'], x['f_acc_err_stat'])} | {f(x['eps_xical_acc'], 3)} | {f(x['fbar'])} | {f(v1)} |")
        L.append("")
        if "toy_truth" in R:
            t = R["toy_truth"]
            L += [f"**Prawda toy (A∨C):** N_acc {t['n_acc']} (mix {e['n_acc']:.1f}, z = {t['z_n_acc']:.2f}), ε {f(t['eps_xical_acc'])} (mix {f(e['eps_xical_acc'])}), "
                  f"f̄ {f(t['fbar'])} (mix {f(e['fbar'])}), f_acc {f(t['f_acc'])} (mix {f(e['f_acc'])}).", ""]
        for key, lab in (("by_mu", "μ (actual_mu)"), ("by_lb", "lumiBlock")):
            L += [f"### f_acc(strona, {lab}) — błąd tylko stat.", "", f"| {lab} | ⟨x⟩ | N_tag A | f_acc A | N_tag C | f_acc C | N_tag A∨C | f_acc A∨C |", "|---|---|---|---|---|---|---|---|"]
            for r in R[key]:
                L.append(f"| {r['lo']:.4g}–{r['hi']:.4g} | {r['mean']:.4g} | {r['A']['n_tag']:.0f} | {pm_(r['A']['f_acc'], r['A']['f_acc_err_stat'])} | {r['C']['n_tag']:.0f} | "
                         f"{pm_(r['C']['f_acc'], r['C']['f_acc_err_stat'])} | {r['AorC']['n_tag']:.0f} | {pm_(r['AorC']['f_acc'], r['AorC']['f_acc_err_stat'])} |")
            L.append("")
        L += ["### Gatunki × pT: tag przypadkowy pod pikiem i w sidebandach (A∨C; kandydaci, bez fitu)", "",
              f"N_tag = kandydaci w oknie w zdarzeniach z tagiem po cięciu ξ_cal; N_acc = Σ p_acc po kandydatach. f_acc ± stat ± syst (syst wspólna, względna {pct(sy['syst_no_stat'])}: pula ⊕ ekstrapolacja ⊕ ½Δ metod). "
              "Sidebandy NIE są próbką „bez prawdziwych tagów” (kombinatoryka masowa w zdarzeniach dyfrakcyjnych niesie prawdziwy proton) — stosunek pik/SB jest tylko informacyjny.", ""]
        for name, S in R["species"].items():
            w = S["windows"]
            v1rows = {(r["pt_lo"], r["pt_hi"]): r for r in v1r.get("species", {}).get(name, {}).get("rows", [])} if v1r else {}
            L += [f"**{name}** — {w['var']}: pik {w['peak']}, SB {w['sb']}", "",
                  "| pT [GeV] | N_tag pik | N_acc pik | f_acc pik ± stat ± syst | f_acc pik A | f_acc pik C | f_acc SB | pik/SB | v1 f_acc pik |", "|---|---|---|---|---|---|---|---|---|"]
            for r in S["rows"]:
                pk, sb = r["peak"], r["sb"]
                v1 = v1rows.get((r["pt_lo"], r["pt_hi"]), {}).get("peak", {}).get("f_acc")
                L.append(f"| {r['pt_lo']:g}–{r['pt_hi']:g} | {pk['n_tag']:.0f} | {pm_(pk['n_acc'], pk['n_acc_err_stat'], 1)} | {pm_(pk['f_acc'], pk['f_acc_err_stat'])} ± {f(pk['f_acc_err_syst'])} | "
                         f"{pm_(r['peak_A']['f_acc'], r['peak_A']['f_acc_err_stat'])} | {pm_(r['peak_C']['f_acc'], r['peak_C']['f_acc_err_stat'])} | {pm_(sb['f_acc'], sb['f_acc_err_stat'])} | "
                         f"{pm_(r.get('peak_over_sb'), r.get('peak_over_sb_err'), 3)} | {f(v1)} |")
            L.append("")
    L += ["## Zmiany względem v1", "",
          "- Estymator: warstwy (LB × nSelBL × „z D/bez D”) zamiast samych LB; region A∨C z pul jednostronnych (v1: pula dwustronna); overlay „union” na parach (bez iteracyjnej poprawki własnego okna); pełna pula na poziomach zapasowych.",
          "- Błąd: bootstrap B = 200 wspólny dla wariantów; warianty puli na pełnej statystyce z błędem różnicy (v1: podpróbka → szum); "
          "nowe składniki: ekstrapolacja w ξ_cal, połowa różnicy metod A–B(i), A–B(ii).",
          "- Bramki: G1 → Poisson LR po warstwowaniu (odporna na zera); G2 → zbiory rozłączne (klasy trygera); I2 (informacyjna): z D vs bez D per klasa nSelBL; G4 usunięta; G5 → kształt w obszarze kontrolnym poza pulą, "
          "Poisson LR z pseudoeksperymentami, próg nominalny p ≥ 0,01; nowe G7 (obszar kontrolny / „przesuwane weto” = A vs B(ii)), G8 (A vs B(i)), G9 (obszar dwuprotonowy); "
          "G3/G6 → informacyjne I3/I6.", "",
          "## Ograniczenia i sprawy otwarte (WSTĘPNE)", "",
          "- **Zakres (W2):** tylko trygery nie-AFP. Dla próbki z trygera L1AFP potrzebni są biorcy reprezentatywni dla tej klasy (decyzja WP2).",
          "- **Obszar ξ_cal < 0,075 pozostaje nietestowalny bezpośrednio** (prawdziwe i przypadkowe tagi nieodróżnialne). G7 sięga do ξ_cal ≈ 0,075, G9 testuje zdarzenia dyfrakcyjne "
          "po drugiej stronie; toy: skok częstości poniżej 0,07 nie jest wykrywany przez żadną bramkę (docs/afp_mix_method.md, tabela mocy).",
          "- Tryger losowy (RD0_FILLED) niedostępny w zapisanych łańcuchach — B(i) jest drugą pulą z łańcuchów L1MBTS, nie absolutnym P(tag) z pustych przecięć.",
          "- Brak score CharmNet/XGB w h5 v2.0 → kandydaci po luźnej selekcji; f_acc pod pikiem dotyczy mieszanki sygnał + kombinatoryka.",
          "- Fluktuacje mianownika (liczba dawców w komórce) pominięte w bootstrapie (względnie ~ częstość tagu, ≲ 10⁻³).", ""]
    open(path, "w").write("\n".join(L))


# ---------------------------------------------------------------------------------------------------------------- toy → format afp_ext
def make_toy(d, scenario="nominal", seed=1, n_per_lb=2500, n_lb=60):
    """Toy afp_mix_toy w formacie parquet afp_ext (jedna próbka na run) + <próbka>_truth.parquet (wiersz = zdarzenie)."""
    from charmbkg import afp_mix_toy as T
    os.makedirs(d, exist_ok=True)
    c, e, tr = T.generate(T.ToyConfig(seed=seed, n_per_lb=n_per_lb, n_lb=n_lb, **T.SCENARIOS[scenario]))
    rng = np.random.default_rng(seed + 100)
    chains = ["HLT_noalg_mb_L1MBTS_1", "HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C", "HLT_noalg_L1J12", RANDOM_CHAIN]
    out = []
    for run in np.unique(e["run_number"]):
        me, mc = (e["run_number"] == run).to_numpy(), (c["run_number"] == run).to_numpy()
        ee, cc, tt = e[me].reset_index(drop=True), c[mc].reset_index(drop=True), tr[me].reset_index(drop=True)
        ev = pd.DataFrame({"event_number": ee["event_number"].astype("u8"), "lumiBlock": ee["lumiBlock"].astype("u4"), "actual_mu": ee["actual_mu"].astype("f4"),
                           "passGRL_MinBiasAFP": np.int8(1), "nPileupVertices": ee["nPileupVertices"].astype("i2"), "nSelBL": ee["nSelBL"].astype("i2")})
        for s in "AC":
            ev[f"xi_cal_{s}"] = ee[f"xi_cal_{s}"].astype("f8")
            ev[f"afp_nProtons_{s}"] = ee[f"afp_nProtons_{s}"].astype("u1")
            ev[f"afp_xi1_{s}"] = ee[f"afp_xi1_{s}"].astype("f8")
            ev[f"afp_xi2_{s}"] = ee[f"afp_xi2_{s}"].astype("f8")
            ev[f"afp_tag_{s}"] = ee[f"afp_tag_{s}"].astype("u1")
        at = tt["afp_trig"].to_numpy()
        tcl = ee[M.COLMAP["trig"]].to_numpy()
        for ch in chains:
            ev[f"trigPassed_{ch}"] = np.int8(0)
        ev["trigPassed_HLT_noalg_mb_L1MBTS_1"] = (~at & (tcl == 0)).astype("i1")
        ev["trigPassed_HLT_noalg_L1J12"] = (~at & (tcl == 1)).astype("i1")
        ev["trigPassed_HLT_noalg_L1AFP_A_OR_C"] = at.astype("i1")
        n = len(cc)
        sp = rng.choice([0, 3, 1], n, p=[0.5, 0.35, 0.15]).astype("u1")
        m = cc[M.COLMAP["m"]].to_numpy()
        sig = cc["is_signal"].to_numpy()
        mass = np.where(sp == 0, m, m - (1.8696 - 1.8648))
        dm = np.where(sp == 1, np.where(sig, rng.normal(0.1454, 0.0006, n), rng.uniform(0.139, 0.160, n)), 0.0)
        ca = pd.DataFrame({"event_number": cc["event_number"].astype("u8"), "Dmeson_species": sp, "Dmeson_mass": mass.astype("f4"),
                           "Dmeson_delta_mass": dm.astype("f4"), "Dmeson_pt": np.clip(cc["Dmeson_pt"].to_numpy(), 3.0, 149.0).astype("f4"),
                           "Dmeson_eta": cc["Dmeson_eta"].astype("f4"), "Dmeson_Lxy": np.float32(1.0)})
        s = f"toy_{run}"
        ev.to_parquet(os.path.join(d, f"{s}_events.parquet"), index=False)
        ca.to_parquet(os.path.join(d, f"{s}_cands.parquet"), index=False)
        tt.to_parquet(os.path.join(d, f"{s}_truth.parquet"), index=False)
        json.dump({"sample": s, "run_number": int(run), "toy": {"scenario": scenario, "seed": seed, "n_per_lb": n_per_lb, "n_lb": n_lb}},
                  open(os.path.join(d, f"{s}.meta.json"), "w"), indent=1)
        out.append(s)
    return out


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser(description="AFP-MIX v2 na danych (WSTĘPNE): estymator A, B(i), B(ii), bramki v2, rozbicie błędu")
    ap.add_argument("--in", dest="indir", default=os.path.join(HERE, "runs", "lxplus"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--samples", default="", help="lista po przecinku; domyślnie wszystkie *_events.parquet w --in")
    ap.add_argument("--make-toy", default="", help="katalog: zapisz toy w formacie afp_ext i zakończ")
    ap.add_argument("--toy-scenario", default="nominal")
    ap.add_argument("--toy-seed", type=int, default=1)
    ap.add_argument("--toy-n-per-lb", type=int, default=2500)
    ap.add_argument("--tag-def", default="nside", choices=["nside", "inwindow"], help="definicja tagu (dane h5 v2.0: nside; toy: inwindow)")
    ap.add_argument("--trigger-mode", default="nonafp", help="nonafp (≥1 zapisany łańcuch bez L1AFP) | chain:<HLT_...>")
    ap.add_argument("--dstar-lxy-min", type=float, default=0.0)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-pseudo", type=int, default=2000)
    ap.add_argument("--n-min-pool", type=int, default=50)
    ap.add_argument("--n-mu-bins", type=int, default=6)
    ap.add_argument("--n-lb-bins", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--save-p", action="store_true", help="zapisz p_tag/p_acc per zdarzenie (szablon do TREx)")
    ap.add_argument("--v1-result", default="", help="result.json z v1 do porównania (domyślnie <in>/out/result.json)")
    a = ap.parse_args(argv)
    if a.make_toy:
        print("toy:", make_toy(a.make_toy, a.toy_scenario, a.toy_seed, a.toy_n_per_lb), "→", a.make_toy)
        return 0
    a.out = a.out or os.path.join(a.indir, "out_v2")
    os.makedirs(a.out, exist_ok=True)
    samples = [x for x in a.samples.split(",") if x] or sorted(os.path.basename(p)[:-len("_events.parquet")] for p in glob.glob(os.path.join(a.indir, "*_events.parquet")))
    if not samples:
        raise SystemExit(f"brak *_events.parquet w {a.indir}")
    t00 = time.time()

    def log(msg):
        print(f"[{time.time() - t00:6.0f} s] {msg}", flush=True)

    RR = []
    for s in samples:
        ev, ca, meta = load_sample(a.indir, s)
        RR.append(analyse_sample(s, ev, ca, meta, a, log))
        del ev, ca
        json.dump(jsonable({"status": "WSTEPNE", "version": "v2", "args": vars(a), "results": RR}), open(os.path.join(a.out, "result_v2.json"), "w"), indent=1)
    rp = os.path.join(a.out, "RAPORT_DANE_v2.md")
    write_report(RR, a, rp, time.time() - t00)
    log(f"raport: {rp}")
    for R in RR:
        if "verdict" in R:
            log(f"{R['sample']}: werdykt {R['verdict']} {R['verdict_failed']}; f_acc(A∨C) = {R['estimate']['f_acc']:.4f} ± {R['estimate']['f_acc'] * R['syst_rel']['total_rel']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
