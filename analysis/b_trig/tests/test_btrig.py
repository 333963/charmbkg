"""Testy B_trig (MC/toy, bez danych). Uruchom: cd charm_work/B_trig && ../../.venv_h5/bin/python -m pytest -q tests"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import btrig_common as C  # noqa: E402
import btrig_emul as EM  # noqa: E402
from btrig_extract import keep_hash, event_cols  # noqa: E402


def test_find_trig_cols_conventions():
    names = ["tp_HLT_noalg_L1J12", "ps_HLT_noalg_L1J12", "trigPassed_HLT_j20_L1MBTS_2", "trigPS_HLT_j20_L1MBTS_2",
             "EventInfo_trigPassed_HLT_mb_sptrk_pt4_L1MBTS_2_NOSYS", "EventInfo_trigPS_HLT_mb_sptrk_pt4_L1MBTS_2"]
    t = C.find_trig_cols(names)
    assert t["HLT_noalg_L1J12"] == {"pass": "tp_HLT_noalg_L1J12", "ps": "ps_HLT_noalg_L1J12"}
    assert t["HLT_j20_L1MBTS_2"]["pass"] == "trigPassed_HLT_j20_L1MBTS_2"
    assert t["HLT_mb_sptrk_pt4_L1MBTS_2"]["pass"].endswith("_NOSYS")
    assert t["HLT_mb_sptrk_pt4_L1MBTS_2"]["ps"] == "EventInfo_trigPS_HLT_mb_sptrk_pt4_L1MBTS_2"


def test_event_cols_h5v20_names():
    names = ["trigPassed_HLT_noalg_L1J12", "trigPS_HLT_noalg_L1J12", "nSelBL", "leadTrackPtBL", "passGRL_MinBiasAFP",
             "nPileupVertices", "actual_mu", "lumiBlock"]
    m = event_cols(names)
    for k in ("tp_HLT_noalg_L1J12", "ps_HLT_noalg_L1J12", "nSelBL", "lead_pt", "grl", "npu", "actual_mu", "lb"):
        assert k in m


def test_classify_disjoint_and_minus_one():
    ev = {"tp_HLT_noalg_L1J12": np.array([1, 1, 0, 0, -1]), "tp_HLT_j20_L1MBTS_2": np.array([0, 1, 1, 0, -1]),
          "tp_HLT_mb_sptrk_pt4_L1MBTS_2": np.array([0, 0, 0, 0, -1])}
    c = C.classify(ev, 428770)
    assert c["J12only"].tolist() == [True, False, False, False, False]
    assert c["MBTSbased"].tolist() == [False, True, True, False, False]
    tot = c["J12only"].astype(int) + c["MBTSbased"] + c["other"]
    assert (tot == 1).all()                       # rozłączne i pokrywają całość
    assert c["other"][4]                          # −1 (MC/brak) ≠ odpalił


def test_classify_427929_grl():
    ev = {"tp_HLT_mb_mbts_L1MBTS_1": np.array([1, 1, 0]), "grl": np.array([1, 0, 1])}
    c = C.classify(ev, 427929)
    assert c["MB427929"].tolist() == [True, False, False]
    assert c["all427929"].tolist() == [True, False, True]


def test_join_event_table():
    T = {"run": np.array([428770, 428770, 428770]), "event": np.array([5, 1, 9]), "tp_X": np.array([1, 0, 1], np.int8)}
    J = C.join_event_table(np.array([428770, 428770, 428770]), np.array([9, 2, 1]), T, ["tp_X"])
    assert J["found"].tolist() == [True, False, True]
    assert J["tp_X"][0] == 1 and J["tp_X"][2] == 0


def test_chi2_and_pulls():
    r = C.chi2_consistency([1.0, 1.0, 1.0], [0.1, 0.2, 0.1])
    assert abs(r["chi2"]) < 1e-12 and r["ndf"] == 2 and r["p"] > 0.99
    r = C.chi2_consistency([0.8, 1.2], [0.1, 0.1])
    assert abs(r["chi2"] - 8.0) < 1e-9
    assert abs(C.pull_subset(0.9, 0.2, 1.0, 0.1) - (-0.1 / np.sqrt(0.03))) < 1e-9
    assert abs(C.pull_indep(0.83, 0.1, 1.0, 0.0) + 1.7) < 1e-9


def test_sigma_scale():
    assert abs(C.sigma_scale(1000, 9000, 0.1, 1000, 9000) - 0.1) < 1e-12
    # 4× S i B → σ/2
    assert abs(C.sigma_scale(1000, 9000, 0.1, 4000, 36000) - 0.05) < 1e-12


def test_keep_hash_deterministic():
    run = np.full(100000, 428770)
    ev = np.arange(100000)
    k1, k2 = keep_hash(run, ev, 0.3), keep_hash(run, ev, 0.3)
    assert (k1 == k2).all() and abs(k1.mean() - 0.3) < 0.01
    assert keep_hash(run, ev, 1.0).all()


def test_effmap_closure_toy():
    """Toy: prawdziwe ε_J12(x) rośnie z lead_pt; mapa z próbki odniesienia odtwarza ε w komórkach (±3σ)."""
    rng = np.random.default_rng(3)
    n = 400000
    nsel = rng.poisson(15, n)
    lead = rng.exponential(2.0, n)
    pj = 1 / (1 + np.exp(-(lead - 4.0)))
    fired = {"J12": rng.random(n) < pj, "sptrk4": rng.random(n) < (lead > 4), "j20MBTS": rng.random(n) < 0.2 * pj,
             "mbts": np.ones(n, bool)}
    em = EM.effmap(nsel, lead, np.ones(n, bool), fired)
    i, j = EM._bin2(nsel, lead)
    for a in range(len(EM.NSEL_EDGES) - 1):
        for b in range(len(EM.LEAD_EDGES) - 1):
            c = (i == a) & (j == b)
            if c.sum() > 2000:
                tru = pj[c].mean()
                assert abs(em["J12"]["eps"][a, b] - tru) < 3 * np.sqrt(tru * (1 - tru) / c.sum()) + 0.02
    # składanie klas z prescale: P(J12only) + P(MBTS) ≤ 1, P(MBTS) ≈ ε_sptrk/50 + ε_j20/50 + 1/5e4
    pjo = EM.class_prob(em, nsel, lead, "J12only")
    pm = EM.class_prob(em, nsel, lead, "MBTSbased")
    assert (pjo + pm <= 1 + 1e-12).all() and pm.max() < 0.05


def test_density_weights_closure():
    rng = np.random.default_rng(5)
    nm, nd = rng.poisson(10, 200000), rng.poisson(14, 200000)
    lm, ld = rng.exponential(1.5, 200000), rng.exponential(2.0, 200000)
    w, _ = EM.density_weights(nm, lm, nd, ld)
    # obsadzenie komórek (nSelBL, lead_pt) po wagach = dane (z konstrukcji); średnia nSelBL zbliża się tylko do szerokości binu
    # (biny nSelBL po 5 → resztkowa różnica ~0,5–1; ograniczenie metody B, PROCEDURA §4)
    i, j = EM._bin2(nm, lm)
    idd, jd = EM._bin2(nd, ld)
    hm = np.zeros((8, 9)); hd = np.zeros((8, 9))
    np.add.at(hm, (i, j), w); np.add.at(hd, (idd, jd), 1.0)
    n_m = np.zeros((8, 9)); np.add.at(n_m, (i, j), 1.0)
    ok = n_m >= 20                                    # komórki z < 20 MC mają wagę 1 (min_mc)
    assert np.allclose(hm[ok] / hm[ok].sum(), hd[ok] / hd[ok].sum(), atol=2e-3)
    assert abs((w * nm).sum() / w.sum() - nd.mean()) < abs(nm.mean() - nd.mean()) / 2


def test_fit_eps_toy_closure():
    """Jeden toy pass/fail D*: ε = 0,79 (SF95 0,83 × ε_MC 0,95), S = 1800, tło ~ próg — fit_passfail odtwarza ε w ±3σ."""
    rng = np.random.default_rng(11)
    shape = {"mu": 145.53, "sc": 0.68, "k": 1.98, "f": 0.66}

    def sig(n):
        core = rng.random(n) < shape["f"]
        return np.where(core, rng.normal(shape["mu"], shape["sc"], n), rng.normal(shape["mu"], shape["sc"] * shape["k"], n))

    def bkg(n):
        u = rng.power(1.5, n) * 30.0                     # ~ u^0.5 na [0, 30] MeV
        return 139.57 + u
    S, eps = 1800, 0.79
    sp = rng.binomial(S, eps)
    dp = np.r_[sig(sp), bkg(20000)] / 1e3
    df_ = np.r_[sig(S - sp), bkg(60000)] / 1e3
    f = C.fit_eps(dp, df_, shape, interval=False)
    assert f["valid"]
    assert abs(f["eps"] - eps) < 3 * f["eps_err"] + 0.02, f
