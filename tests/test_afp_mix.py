"""AFP-MIX v2: tło przypadkowego tagu AFP z mieszania protonów — jednostki, domknięcie na toy MC, rozbicie błędu
i testy negatywne (scenariusze naruszające założenia: bramki v2 muszą wykryć to, co deklarują w docs/afp_mix_method.md)."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from charmbkg import afp_mix as M
from charmbkg import afp_mix_toy as T

TOY = dict(n_runs=2, n_lb=60, n_per_lb=2500)
GC = dict(n_boot=100, n_pseudo=300)


@pytest.fixture(scope="module")
def toy():
    C, E, tr = T.prepare(T.ToyConfig(seed=11, **TOY))
    tgt = M.events_from_candidates(C)
    hd = tr[tr.hasD].reset_index(drop=True)
    assert len(tgt) == len(hd) and np.array_equal(tgt["event"].to_numpy(), E["event"].to_numpy()[tr.hasD.to_numpy()])
    return C, E, tr, tgt, hd


@pytest.fixture(scope="module")
def report(toy):
    C, E, *_ = toy
    return M.run_procedure(C, E, gc=M.GateConfig(**GC), keep=True)


def _gate(rep, name):
    return next(g for g in rep["gates"] if g["gate"] == name)


def _flat_events(n=300, xical=0.5):
    return pd.DataFrame({"run": 1, "event": np.arange(n), "lb": 1, "mu": 0.02, "nsel": 30, "hasD": True, "np_A": 0, "np_C": 0, "xi1_A": -1.0,
                         "xi1_C": -1.0, "xi2_A": -1.0, "xi2_C": -1.0, "xical_A": xical, "xical_C": xical})


# ---------------------------------------------------------------------------------------------------------------- jednostki
def test_colmap_single_place():
    df = pd.DataFrame({v: [0.0] for v in M.COLMAP.values()}).rename(columns={"actual_mu": "mu_actual"})
    with pytest.raises(KeyError):
        M.standardize(df)
    out = M.standardize(df, colmap={"mu": "mu_actual"})
    assert "mu" in out and "xical_A" in out and "tag_C" in out and "nsel" in out and "m" in out
    assert M.COLMAP["m"] == "Dmeson_mass"      # h5 v2.0 (przegląd D5)


def test_tag_window_both_definitions():
    n = np.array([0, 1, 1, 2, 2, 3, 3])
    xi1 = np.array([-1, 0.05, 0.20, 0.01, 0.04, 0.01, 0.01])
    xi2 = np.array([-1, -1, -1, 0.05, 0.06, 0.05, 0.12])
    tag, anyin = M.tag_window(n, xi1, xi2, 0.03, 0.10)
    assert tag.tolist() == [False, True, False, True, False, False, False]   # n=3, ξ2 ≤ hi → niejednoznaczne
    assert anyin.tolist() == [False, True, False, True, True, True, False]
    tag_ns, _ = M.tag_window(n, xi1, xi2, 0.03, 0.10, "nside")              # dane: dokładnie 1 proton po stronie
    assert tag_ns.tolist() == [False, True, False, False, False, False, False]


def test_region_pass_charmpp():
    t, f = True, False
    assert M.region_pass(np.array([t]), np.array([t]), np.array([t]), np.array([f]), "AorC").tolist() == [False]
    assert M.region_pass(np.array([t]), np.array([f]), np.array([t]), np.array([f]), "AorC").tolist() == [True]
    assert M.region_pass(np.array([f]), np.array([t]), np.array([f]), np.array([t]), "C").tolist() == [True]


def test_self_exclusion_also_with_bootstrap():
    """Biorca będący jedynym otagowanym dawcą w komórce nie daje sobie tagu (v1: pomijane przy wagach bootstrapu, D1)."""
    ev = _flat_events()
    ev.loc[0, ["np_A", "xi1_A"]] = [1, 0.05]
    cfg = M.MixConfig(use_tag_column=False, n_min_pool=10)
    mx = M.mix(ev.iloc[[0, 1]], ev, cfg, "A")
    assert mx.p_tag[0] == 0.0 and math.isclose(mx.p_tag[1], 1 / (len(ev) - 1))
    boot = M.Boot.for_events(ev, cfg, B=50, seed=1)
    reps = mx.boot_p_tag(boot)
    assert np.allclose(reps, boot.get(M.ekey(ev.iloc[[0]]))[0] / (len(ev) - 1))   # tylko biorca 1 dostaje wkład


def test_union_overlay_and_levels():
    """(a) biorca z własnym protonem po stronie nie dostaje tagu przypadkowego (overlay union zamiast poprawki C z v1);
    (b) komórka z < n_min_pool dawców → poziom zapasowy z PEŁNĄ pulą poziomu (przegląd D7)."""
    ev = _flat_events(400)
    ev.loc[:9, ["np_A", "xi1_A"]] = [1, 0.06]                       # 10 otagowanych dawców w puli
    ev["lb"] = np.where(np.arange(400) < 380, 1, 2)                 # LB 2: tylko 20 dawców
    tg = _flat_events(4, xical=0.05)
    tg["event"] += 10_000
    tg.loc[1, ["np_A", "xi1_A"]] = [1, 0.20]                        # własny proton poza oknem (n_side → brak tagu po dołożeniu)
    tg.loc[2, "xical_A"] = 0.11                                     # Δ = 0,05 > δ → nie przechodzi cięcia, ale leży w obszarze kontrolnym
    tg.loc[3, "lb"] = 2
    cfg = M.MixConfig(use_tag_column=False, n_min_pool=50, tag_def="nside")
    mx = M.mix(tg, ev, cfg, "A")
    p = mx.p_acc()
    assert math.isclose(p[0], 10 / 380 / (1 - 10 / 380), rel_tol=1e-6)   # 1/(1 − q): biorcy z już dołożonym protonem wypadają z sumy
    assert p[1] == 0.0 and p[2] == 0.0
    obs, pm = M.control_region(mx, "A", cfg)
    assert not obs.any() and np.unique(mx.pt[pm]).tolist() == [2]
    assert mx.level.tolist() == [0, 0, 0, 1] and mx.n_pool[3] == 400     # LB 2 → poziom μ: cała pula, nie tylko ubogie LB
    assert math.isclose(p[3], 10 / 400 / (1 - 10 / 400), rel_tol=1e-6)


def test_stat_tools_robust_to_zeros():
    """Poisson LR nie wybucha na pustych przedziałach (v1 G5: χ² Neymana ~ 10⁵, v1 G1: przedział z 0 tagów dominował)."""
    exp = np.array([13.7, 6.0, 2.0, 0.5])
    assert M.poisson_lr(exp, exp) == pytest.approx(0)
    assert 0 < M.poisson_lr([14, 0, 0, 0], exp) < 30
    lab = M.merge_groups(exp, 5.0)
    assert lab.tolist() == [0, 1, 1, 1] and np.bincount(lab, exp).min() >= 5
    z = M.count_z(7, 13.7)                       # przegląd W4: P(≤ 7 | 13,7) = 0,036 → 1,8σ, nie −2,55
    assert -2.0 < z < -1.6
    assert M.count_z(120, 100, 5.0) == pytest.approx(20 / math.sqrt(125))
    assert abs(M.count_z(7, 13.7, 5.0)) < abs(z)  # niepewność przewidywania zmniejsza istotność


# ---------------------------------------------------------------------------------------------------------------- domknięcie toy
def test_closure_n_acc_eps_ndiff(report, toy):
    *_, hd = toy
    est, sy = report["estimate"], report["syst_rel"]
    n_true = hd.acc_region.sum()
    z = (est["n_acc"] - n_true) / math.sqrt(n_true + est["n_acc_err_stat"] ** 2)
    assert abs(z) < 3, (est["n_acc"], n_true)
    eps_true = n_true / hd.acc_tag_nocut.sum()
    assert abs(est["eps_xical_acc"] - eps_true) < 3 * eps_true / math.sqrt(n_true), (est["eps_xical_acc"], eps_true)
    n_own = hd.own_region.sum()
    assert abs(est["n_diff"] - n_own) < 3 * est["n_diff_err"], (est["n_diff"], n_own)
    fb_true = hd.lost.sum() / n_own
    assert abs(est["fbar"] - fb_true) < 3 * math.sqrt(max(hd.lost.sum(), 1)) / n_own
    # drugie estymatory zgodne z prawdą w granicach własnych błędów
    for k in ("method_Bi", "method_Bii"):
        b = report[k]
        assert b["available"] and abs(b["n_acc"] - n_true) < 3 * math.hypot(b["n_acc_err"], math.sqrt(n_true)), (k, b["n_acc"], n_true)
    assert set(sy) >= {"stat", "pool", "extrapolation", "half_diff_methods", "total_rel"}
    assert 0 < sy["stat"] < 0.15 and sy["total_rel"] < 0.4


def test_all_gates_pass_nominal(report):
    bad = [g for g in report["gates"] if g["gate"] in M.VERDICT_GATES and g["status"] != "PASS"]
    assert not bad, bad
    assert report["verdict"] == "PASS"
    assert {g["gate"] for g in report["gates"]} >= set(M.VERDICT_GATES)
    assert all(g["tests"] and g["not_tests"] for g in report["gates"])      # każda bramka mówi, co testuje i czego NIE
    assert not any(g["gate"].startswith("G4") for g in report["gates"])     # G4 z v1 usunięta
    assert all(g["status"] == "INFO" for g in report["gates"] if g["gate"].startswith("I"))   # G3/G6 → informacyjne
    assert M.GateConfig().p_min == 0.01                                     # próg G5 nominalny (v1: 0,002 post hoc)


def test_variants_full_stat_shared_bootstrap(report):
    """Warianty puli: pełna statystyka, błąd RÓŻNICY ze wspólnego bootstrapu (W3) — mniejszy niż niezależne błędy."""
    v = report["variants"]
    assert {"margin_0", "margin_0.05", "margin_0.1", "strata_mu", "act_coarse", "plus_trig_class"} <= set(v)
    assert {"donors_D_and_noD", "v1_like_LB_only"} <= set(report["info_variants"])     # informacyjne, nie w systematyce
    s = report["syst_rel"]["stat"]
    assert v["strata_mu"]["delta_err"] < s and v["plus_trig_class"]["delta_err"] < s
    assert abs(v["strata_mu"]["delta"]) < 0.1
    assert report["syst_rel"]["pool"] == pytest.approx(max(x["significant"] for x in v.values()))


def test_xi_distribution_after_cut(report, toy):
    """Rozkład ξ_AFP przypadkowych po cięciu ξ_cal (korelacja z ξ_cal biorcy) — mieszanie vs prawda, strona A."""
    C, E, tr, tgt, hd = toy
    mx = report["_obj"]["mx"]["A"]
    edges = np.array([0.03, 0.045, 0.06, 0.075, 0.09, 0.10])
    hm = M.xi_distribution(mx, "A", edges)
    acc = hd.acc_region_A.to_numpy() & ~hd.own_tag_A.to_numpy()
    xt = tgt["xi1_A"].to_numpy()[acc]
    xt = np.where(xt < 0.03, tgt["xi2_A"].to_numpy()[acc], xt)        # proton tagu = w oknie
    ht, _ = np.histogram(xt, edges)
    assert abs(hm.sum() - ht.sum()) < 3 * math.sqrt(ht.sum()) + 0.1 * ht.sum()
    s = ht.sum() / hm.sum()
    chi2 = float(((ht - s * hm) ** 2 / np.maximum(s * hm, 1)).sum())
    assert chi2 < 20, (ht, s * hm)
    hm0 = M.xi_distribution(mx, "A", edges, after_cut=False)
    mid = 0.5 * (edges[1:] + edges[:-1])
    assert (hm * mid).sum() / hm.sum() > (hm0 * mid).sum() / hm0.sum()   # cięcie przesuwa ξ w górę


def test_peak_vs_sideband_template(report, toy):
    C, E, tr, tgt, hd = toy
    t = M.acc_template(C, tgt, report["_obj"]["mx"]["AorC"].p_acc(), np.linspace(1.7, 2.2, 51),
                       windows={"peak": (1.84, 1.90), "sb_lo": (1.72, 1.80), "sb_hi": (1.94, 2.10)})
    mp = t["mean_p_acc"]
    assert abs(mp["peak"] / mp["sb_lo"] - 1) < 0.15 and abs(mp["sb_hi"] / mp["sb_lo"] - 1) < 0.15
    assert math.isclose(t["hist"].sum(), t["w"][(C["m"] >= 1.7) & (C["m"] < 2.2)].sum(), rel_tol=1e-9)


# ---------------------------------------------------------------------------------------------------------------- testy negatywne
def test_naive_pool_overestimates(toy):
    """Pula bez weta ξ_cal (dawcy niosą własną dyfrakcję): N_acc zawyżone kilkukrotnie."""
    C, E, tr, tgt, hd = toy
    n_naive = M.mix_region(tgt, E, M.MixConfig(donor_mode="all"), "AorC")[1].total()
    assert n_naive > 2.5 * hd.acc_region.sum()


def test_activity_strata_remove_bias():
    """Częstość zależna od nSelBL: bez warstw aktywności (≈ v1) N_acc zaniżone, v2 domyka; bramki nie alarmują."""
    rep, tr = T.run_scenario("rate_vs_nsel", seed=12, gc=M.GateConfig(**GC), **TOY)
    nA, n_true = rep["estimate"]["n_acc"], tr["n_acc"]
    assert rep["info_variants"]["v1_like_LB_only"]["ratio"] < 0.93
    assert abs(nA - n_true) < 3 * math.hypot(rep["estimate"]["n_acc_err_stat"], math.sqrt(n_true))
    assert _gate(rep, "G2_disjoint_closure")["status"] == "PASS"


@pytest.mark.parametrize("scenario", ["pool_noD_x2", "independence_broken"])
def test_pool_composition_removed_by_construction(scenario):
    """W1: częstość zależna od obecności D (dawcy bez D ×2 / zdarzenia z D ×1,6). v1: bias ±30–60 %, wykrywany tylko przez G2.
    v2: klasa „z D/bez D” w warstwach → A domyka; I2 (informacyjna) pokazuje różnicę; B(i) (próbka bez D) odstaje → G8 FAIL."""
    rep, tr = T.run_scenario(scenario, seed=12, gc=M.GateConfig(**GC), **TOY)
    nA, n_true = rep["estimate"]["n_acc"], tr["n_acc"]
    assert abs(nA - n_true) < 3 * math.hypot(rep["estimate"]["n_acc_err_stat"], math.sqrt(n_true)), (nA, n_true)
    assert abs(rep["info_variants"]["donors_D_and_noD"]["ratio"] - 1) > 0.15
    assert _gate(rep, "I2_hasD_vs_noD")["p"] < 1e-3
    assert _gate(rep, "G8_methods_A_vs_Bi")["status"] == "FAIL" and rep["syst_rel"]["half_diff_methods"] > 0.1


@pytest.mark.parametrize("scenario,gate", [("rate_vs_trig", "G2_disjoint_closure"), ("xical_pow_0.3", "G1_pool_homogeneity"),
                                           ("xical_step_0.12_x2", "G7_methods_A_vs_Bii")])
def test_gate_detects_violation(scenario, gate):
    """Scenariusze naruszające założenia → wskazana bramka FAIL i werdykt FAIL (toy pełnego rozmiaru dla mocy G7)."""
    kw = dict(n_runs=2, n_lb=100, n_per_lb=4000) if gate.startswith("G7") else TOY
    rep, tr = T.run_scenario(scenario, seed=12, gc=M.GateConfig(**GC), **kw)
    assert _gate(rep, gate)["status"] == "FAIL", _gate(rep, gate)
    assert rep["verdict"] == "FAIL" and gate in rep["verdict_failed"]
    if scenario == "xical_step_0.12_x2":      # B(ii) koryguje w stronę prawdy, a połowa różnicy metod wchodzi do błędu
        assert abs(rep["method_Bii"]["n_acc"] / tr["n_acc"] - 1) < abs(rep["estimate"]["n_acc"] / tr["n_acc"] - 1)
        assert rep["syst_rel"]["half_diff_methods"] > 0.1
