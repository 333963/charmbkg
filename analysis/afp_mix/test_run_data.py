"""Test run_data.py v2 na toy wyeksportowanym do formatu parquet afp_ext. Uruchomienie (wymaga pyarrow):
  python -m pytest -q analysis/afp_mix/test_run_data.py"""
import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_data as RD  # noqa: E402
from charmbkg import afp_mix as M  # noqa: E402


def _run(tmp, scenario):
    d = str(tmp / scenario)
    RD.make_toy(d, scenario, seed=4, n_per_lb=2500)   # seed 3: G2 z = 3,1 w jednej z 2 próbek (fluktuacja: 1 na 26 biegów nominalnych)
    assert RD.main(["--in", d, "--out", d + "/out", "--tag-def", "inwindow", "--n-boot", "100", "--n-pseudo", "300"]) == 0
    txt = open(d + "/out/RAPORT_DANE_v2.md").read()
    assert "WSTĘPNE" in txt and "Zakres ważności" in txt and "czego NIE testuje" in txt
    return d, json.load(open(d + "/out/result_v2.json"))["results"]


@pytest.fixture(scope="module")
def nominal(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("toy"), "nominal")


def test_closure_and_gates(nominal):
    _, res = nominal
    for R in res:
        assert R["verdict"] == "PASS", R["verdict_failed"]
        t, e, sy = R["toy_truth"], R["estimate"], R["syst_rel"]
        assert abs(t["z_n_acc"]) < 3
        assert abs(e["eps_xical_acc"] - t["eps_xical_acc"]) < 0.03
        assert abs(e["fbar"] - t["fbar"]) < 3 * (max(t["fbar"], 1e-3) / e["n_tag"]) ** 0.5
        assert all(k in sy for k in ("stat", "pool", "extrapolation", "half_diff_methods", "total_rel")) and sy["total_rel"] < 0.5
        assert R["method_Bi"]["available"] and R["method_Bii"]["n_acc"] > 0
        assert R["inputs"]["n_random_trigger_RD0"] == 0 and "bez D" in R["inputs"]["b1_definition"]
        assert set(R["species"]) == {"Dplus", "D0", "Dstar"} and len(R["species"]["Dplus"]["rows"]) == len(RD.PT_EDGES)
        assert {g["gate"] for g in R["gates"]} >= set(M.VERDICT_GATES)


def test_matches_module_procedure(nominal):
    """N_acc z run_data (parquet, klasy trygera z trigPassed_*) == afp_mix.mix na tych samych tabelach."""
    d, res = nominal
    s = res[0]["sample"]
    ev, ca, meta = RD.load_sample(d, s)
    E = RD.standardize_events(ev)
    E["hasD"] = E["event"].isin(ca["event_number"]).to_numpy()
    T = E[E["hasD"]].reset_index(drop=True)
    mx = M.mix_region(T, E, M.MixConfig(tag_def="inwindow"), "AorC")[1]
    e = res[0]["estimate"]
    assert e["n_tag"] == float(mx.own.sum())
    assert abs(e["n_acc"] - mx.total()) < 1e-6 * mx.total()


def test_species_sum_consistent(nominal):
    """Suma po binach pT == bin inkluzywny (N_tag i N_acc — estymator liniowy w wagach biorców)."""
    _, res = nominal
    rows = res[0]["species"]["Dplus"]["rows"]
    assert abs(sum(r["peak"]["n_tag"] for r in rows[1:]) - rows[0]["peak"]["n_tag"]) < 1e-9
    assert abs(sum(r["peak"]["n_acc"] for r in rows[1:]) - rows[0]["peak"]["n_acc"]) < 1e-6


def test_afp_trigger_events_excluded(tmp_path):
    """Zdarzenia wyłącznie z trygera L1AFP (toy: 1 % zdarzeń z D) są wykluczane (zakres W2) → bramki domykają się."""
    d, res = _run(tmp_path, "afp_trigger")
    for R in res:
        assert R["inputs"]["n_afp_only"] > 0 and R["inputs"]["afp_only_tag_rate"] > 0.5
        g2 = next(g for g in R["gates"] if g["gate"] == "G2_disjoint_closure")
        assert g2["status"] == "PASS", g2
        assert abs(R["toy_truth"]["z_n_acc"]) < 3


def test_spoil_nside_ge_inwindow():
    """Definicja n_side == 1: każdy dodatkowy proton po otagowanej stronie psuje tag → f̄(n_side) ≥ f̄(1 w oknie)."""
    n = 400
    D = pd.DataFrame({"run": 1, "event": np.arange(n), "lb": 1, "mu": 0.02, "np_A": 0, "np_C": 0, "xi1_A": -1.0, "xi1_C": -1.0,
                      "xi2_A": -1.0, "xi2_C": -1.0, "xical_A": 0.5, "xical_C": 0.5})
    D.loc[:39, ["np_A", "xi1_A"]] = [1, 0.15]          # 10 % dawców: proton A poza oknem (większe ξ)
    D.loc[40:59, ["np_A", "xi1_A"]] = [1, 0.05]        # 5 %: proton A w oknie
    T = D.iloc[[399]].copy()
    T[["np_A", "xi1_A", "xical_A"]] = [1, 0.06, 0.01]
    T["event"] = 10_000
    T = T.reset_index(drop=True)
    kw = dict(use_tag_column=False, n_min_pool=10, donor_mode="all")
    q_ns = M.spoil_probability(T, D, M.MixConfig(tag_def="nside", **kw), "A")[0]
    q_in = M.spoil_probability(T, D, M.MixConfig(tag_def="inwindow", **kw), "A")[0]
    assert abs(q_ns - 60 / n) < 1e-12 and abs(q_in - 20 / n) < 1e-12
