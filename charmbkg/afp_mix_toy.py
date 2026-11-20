"""AFP-MIX v2: toy MC dla metody mieszania protonów (``charmbkg.afp_mix``) z wstrzykniętą dyfrakcją i protonami przypadkowymi.

Zdarzenie = oddziaływanie z trygera (ND/SD_A/SD_C/DD/CD) ⊕ pile-up Poisson(μ) oddziaływań tej samej mieszanki ⊕ halo
(Bernoulli h_LB po stronie, ξ płaskie). μ maleje w LB (fill), h_LB losowe per LB (warunki wiązki) → test warstw LB.
Protony: SD_s/CD dają proton po stronie s z ξ log-płaskim [0,003; 0,3]; akceptancja AFP 0,015–0,20 (ε = 0,95), σ_ξ = 0,002.
ξ_cal: SD_s → ξ·U(0,3; 1) + |N(0; 0,003)| po stronie s (przechodzi cięcie), reszta lognormal (ND mediana 0,25,
zdarzenia z D ×1,4). Pile-up nie zmienia ξ_cal (weto wierzchołków).
v2: nSelBL skorelowane z ξ_cal (λ = 4 + 40·(ξ_cal,A·ξ_cal,C)^¼, ×1,2 z D) — pula dawców (duże ξ_cal) ma wyższą aktywność
niż biorcy z małym ξ_cal, jak w danych; klasa trygera ``trig_class`` (1 = „dżetowy”, częstszy przy dużym nSelBL);
flaga ``has_D``; próbka B(i) ``b1_sample`` = klasa 0 bez D.
**Scenariusze naruszające założenia** (v1 spełniał je z konstrukcji; przegląd W7): częstość halo zależna od ξ_cal po stronie
(skok poniżej progu / potęgowo), od nSelBL, od klasy trygera, od obecności D; inne widmo ξ halo przy małym ξ_cal;
rozmycie ξ_cal (przeciek własnych protonów ponad cięcie). Zob. ``SCENARIOS`` i ``power_table``.
Kandydaci (tylko zdarzenia z D): m = pik D⁺ (8 %) + tło exp; η* > 0 z P = 0,6 w SD, symetrycznie w pozostałych.
Wyjście w NAZWACH ZEWNĘTRZNYCH (jak h5 v2.0) — test mapowania ``afp_mix.COLMAP``.
Prawda: ``own_region`` (tag tylko z własnych protonów), ``obs_region``, ``acc_region`` = obs ∧ ¬own (przypadkowy),
``acc_tag_nocut`` (przypadkowy bez cięcia ξ_cal), ``lost`` = own ∧ ¬obs (strata przez „dokładnie 1”).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .afp_mix import COLMAP, region_pass

M_DP = 1.8696
TYPES = ("ND", "SDA", "SDC", "DD", "CD")


@dataclass
class ToyConfig:
    n_runs: int = 2
    n_lb: int = 100
    n_per_lb: int = 4000
    mu_hi: float = 0.05
    mu_lo: float = 0.005
    p_type: tuple = (0.78, 0.07, 0.07, 0.07, 0.01)
    p_hasD: tuple = (0.30, 0.15, 0.15, 0.20, 0.15)
    halo0: float = 0.02            # P(proton halo po stronie, ξ ∈ [0; 0,25]) na BC, mnożone lognormal(0; 0,4) per LB
    halo_scale_D: float = 1.0      # halo ×K w zdarzeniach z D (niezależność złamana; G2)
    halo_scale_noD: float = 1.0    # halo ×K w zdarzeniach bez D (skład puli, W1; G2)
    pile_scale: float = 1.0        # skala liczby oddziaływań pile-up
    afp_trig_frac: float = 0.0     # ułamek zdarzeń z D dodanych „z trygera AFP” (wymuszony proton) — poza zakresem metody (G0)
    # --- naruszenia założeń (v2)
    halo_xical_step: tuple = (0.0, 1.0)   # (próg, K): halo ×K po stronie s, gdy ξ_cal,s < próg (K1 przeglądu)
    halo_xical_pow: float = 0.0           # halo × clip((ξ_cal,s/0,25)^(−α); 0,4; 3)
    halo_nsel: float = 0.0                # halo × (1 + a·exp(−nSelBL/8))
    halo_trig: float = 1.0                # halo ×K w klasie trygera 1
    halo_xi_low: tuple = (0.0, 0.0)       # (próg, t): przy ξ_cal,s < próg widmo ξ halo ∝ 1 − t·ξ/0,25 (zamiast płaskiego)
    xical_smear: float = 0.0              # lognormalne rozmycie ξ_cal własnej strony SD/CD (przeciek ponad cięcie)
    xi_lo: float = 0.03
    xi_hi: float = 0.10
    delta: float = 0.02
    seed: int = 1


def _loguni(rng, n, a, b):
    return np.exp(rng.uniform(np.log(a), np.log(b), n))


def _lognorm(rng, n, med, s):
    return med * np.exp(rng.normal(0, s, n))


def _own_protons(rng, typ):
    """ξ prawdziwe po stronach A, C (NaN = brak) dla oddziaływań o typach ``typ``."""
    n = len(typ)
    xi = {s: np.full(n, np.nan) for s in "AC"}
    for s, codes in (("A", (1, 4)), ("C", (2, 4))):
        m = np.isin(typ, codes)
        xi[s][m] = _loguni(rng, m.sum(), 0.003, 0.3)
    return xi


def _seen(rng, xi):
    ok = np.isfinite(xi) & (xi >= 0.015) & (xi <= 0.20) & (rng.random(len(xi)) < 0.95)
    out = np.where(ok, xi + rng.normal(0, 0.002, len(xi)), np.nan)
    return out


def _summ(P, lo, hi):
    """P: (n, k) ξ z NaN → n, ξ1, ξ2 (−1 gdy brak), tag (dokładnie 1 w oknie)."""
    Ps = np.sort(P, 1)
    n = np.isfinite(P).sum(1)
    xi1 = np.where(n >= 1, Ps[:, 0], -1.0)
    xi2 = np.where(n >= 2, Ps[:, 1] if P.shape[1] > 1 else -1.0, -1.0)
    nin = ((P >= lo) & (P <= hi)).sum(1)
    return n, xi1, xi2, nin == 1, nin


def generate(tc: ToyConfig | None = None):
    """→ (cands, events, truth): kandydaci i zdarzenia w nazwach h5 v2.0, prawda per zdarzenie (index = events)."""
    tc = tc or ToyConfig()
    rng = np.random.default_rng(tc.seed)
    nlb = tc.n_runs * tc.n_lb
    N = nlb * tc.n_per_lb
    run = np.repeat(np.arange(tc.n_runs) + 428770, tc.n_lb * tc.n_per_lb)
    lb_i = np.repeat(np.arange(nlb), tc.n_per_lb)
    lb = lb_i % tc.n_lb + 1
    frac = (lb - 1) / max(tc.n_lb - 1, 1)
    mu_lb = tc.mu_hi * (tc.mu_lo / tc.mu_hi) ** frac
    mu = mu_lb * (1 + rng.normal(0, 0.02, N))
    halo_lb = tc.halo0 * np.exp(rng.normal(0, 0.4, nlb))
    halo = halo_lb[lb_i]

    typ = rng.choice(len(TYPES), N, p=tc.p_type)
    hasD = rng.random(N) < np.asarray(tc.p_hasD)[typ]
    xi_true = _own_protons(rng, typ)

    # ξ_cal
    xical = {}
    for s, sd in (("A", 1), ("C", 2)):
        x = _lognorm(rng, N, 0.25, 0.9)
        x = np.where(typ == 3, _lognorm(rng, N, 0.15, 0.9), x)
        x = np.where((typ == 1) | (typ == 2), _lognorm(rng, N, 0.30, 0.9), x)
        x = np.where(hasD & (typ == 0), x * 1.4, x)
        own = (typ == sd) | (typ == 4)
        xo = xi_true[s] * rng.uniform(0.3, 1.0, N) + np.abs(rng.normal(0, 0.003, N))
        if tc.xical_smear > 0:
            xo = xo * np.exp(rng.normal(0, tc.xical_smear, N))
        x = np.where(own, xo, x)
        xical[s] = x

    # aktywność i klasa trygera (v2)
    lam = 4 + 40 * (np.clip(xical["A"], 0, 1) * np.clip(xical["C"], 0, 1)) ** 0.25
    nsel = rng.poisson(lam * np.where(hasD, 1.2, 1.0))
    trig = (rng.random(N) < 0.6 / (1 + np.exp(-(nsel - 28) / 6.0))).astype(np.int8)

    # protony: własny + pile-up + halo
    npu = rng.poisson(mu * tc.pile_scale)
    K = 4
    acc = {s: np.full((N, K), np.nan) for s in "AC"}
    own_seen = {s: _seen(rng, xi_true[s]) for s in "AC"}
    for j in range(K - 1):
        m = npu > j
        t_pu = rng.choice(len(TYPES), m.sum(), p=tc.p_type)
        xp = _own_protons(rng, t_pu)
        for s in "AC":
            acc[s][m, j] = _seen(rng, xp[s])
    hscale = np.where(hasD, tc.halo_scale_D, tc.halo_scale_noD) * (1 + tc.halo_nsel * np.exp(-nsel / 8.0)) * np.where(trig == 1, tc.halo_trig, 1.0)
    afp_trig = hasD & (rng.random(N) < tc.afp_trig_frac)
    for s in "AC":
        hs = hscale.copy()
        if tc.halo_xical_step[1] != 1.0:
            hs = hs * np.where(xical[s] < tc.halo_xical_step[0], tc.halo_xical_step[1], 1.0)
        if tc.halo_xical_pow:
            hs = hs * np.clip((xical[s] / 0.25) ** (-tc.halo_xical_pow), 0.4, 3.0)
        hm = (rng.random(N) < np.clip(halo * hs, 0, 1)) | (afp_trig & (s == "A"))
        xh = rng.uniform(0.0, 0.25, hm.sum())
        if tc.halo_xi_low[1]:
            thr, t = tc.halo_xi_low
            u = rng.random(hm.sum())           # gęstość ∝ 1 − t·x (x = ξ/0,25): odwrócenie dystrybuanty
            x_alt = (1 - np.sqrt(1 - u * t * (2 - t))) / t
            xh = np.where(xical[s][hm] < thr, 0.25 * x_alt, xh)
        acc[s][hm, K - 1] = _seen(rng, xh)
        if s == "A":   # wymuszony proton AFP-trygera w oknie
            acc[s][afp_trig, K - 1] = rng.uniform(tc.xi_lo, tc.xi_hi, afp_trig.sum())

    ev = {"run_number": run.astype(np.uint32), "event_number": np.arange(N, dtype=np.uint64) + 1, "lumiBlock": lb,
          "actual_mu": mu, "nPileupVertices": np.zeros(N, np.int32), "nSelBL": nsel.astype(np.int32),
          COLMAP["hasD"]: hasD, COLMAP["trig"]: trig, COLMAP["b1"]: (trig == 0) & ~hasD}
    tr = {"type": np.asarray(TYPES)[typ], "hasD": hasD, "npu": npu, "afp_trig": afp_trig}
    T, Town, Cut, Cown = {}, {}, {}, {}
    for s in "AC":
        P = np.concatenate([own_seen[s][:, None], acc[s]], 1)
        n, xi1, xi2, tag, nin = _summ(P, tc.xi_lo, tc.xi_hi)
        no, xo1, _, tago, _ = _summ(own_seen[s][:, None], tc.xi_lo, tc.xi_hi)
        ev.update({f"afp_nProtons_{s}": n.astype(np.int32), f"afp_xi1_{s}": xi1, f"afp_xi2_{s}": xi2,
                   f"afp_tag_{s}": tag, f"xi_cal_{s}": xical[s]})
        T[s], Town[s] = tag, tago
        Cut[s] = xical[s] <= xi1 + tc.delta
        Cown[s] = xical[s] <= xo1 + tc.delta
        tr[f"own_tag_{s}"] = tago
        tr[f"any_acc_inwin_{s}"] = ((acc[s] >= tc.xi_lo) & (acc[s] <= tc.xi_hi)).any(1)
    ones = np.ones(N, bool)
    obs = region_pass(T["A"], T["C"], Cut["A"], Cut["C"], "AorC")
    own = region_pass(Town["A"], Town["C"], Cown["A"], Cown["C"], "AorC")
    obs_nc = region_pass(T["A"], T["C"], ones, ones, "AorC")
    own_nc = region_pass(Town["A"], Town["C"], ones, ones, "AorC")
    tr.update({"obs_region": obs, "own_region": own, "acc_region": obs & ~own, "lost": own & ~obs,
               "acc_tag_nocut": obs_nc & ~own_nc})
    for s in "AC":
        o = T[s] & Cut[s]
        w = Town[s] & Cown[s]
        tr[f"acc_region_{s}"] = o & ~w
    events = pd.DataFrame(ev)
    truth = pd.DataFrame(tr)

    # kandydaci
    idx = np.flatnonzero(hasD)
    nc = 1 + rng.poisson(1.5, len(idx))
    ci = np.repeat(idx, nc)
    nC = len(ci)
    sig = rng.random(nC) < 0.08
    u = rng.random(nC)
    lam = 3.0
    mb = 1.70 - np.log(1 - u * (1 - np.exp(-lam * 0.5))) / lam
    m = np.where(sig, rng.normal(M_DP, 0.015, nC), mb)
    aeta = rng.uniform(0, 2.2, nC)
    t_c = typ[ci]
    sgn = np.where(rng.random(nC) < 0.5, 1.0, -1.0)
    sd = (t_c == 1) | (t_c == 2)
    sgn = np.where(sd, np.where(rng.random(nC) < 0.6, 1.0, -1.0), sgn)   # η* > 0 dla 60 %
    eta = np.where(t_c == 1, -sgn * aeta, sgn * aeta)                     # tag A: η = −η*; C: η = +η*
    cands = events.iloc[ci].reset_index(drop=True)
    cands[COLMAP["m"]] = m
    cands[COLMAP["eta"]] = eta
    cands[COLMAP["pt"]] = _lognorm(rng, nC, 4.0, 0.5)
    cands["is_signal"] = sig
    return cands, events, truth


SCENARIOS = {
    "nominal": {},
    # --- naruszenia, które estymator v2 usuwa z konstrukcji (warstwy aktywności / trygera)
    "rate_vs_nsel": {"halo_nsel": 3.0},                    # częstość zależna od aktywności (nSelBL)
    "rate_vs_trig": {"halo_trig": 1.8},                    # częstość zależna od klasy trygera
    "pool_noD_x2": {"halo_scale_noD": 2.0},                # dawcy bez D mają inną częstość (W1)
    # --- naruszenia w ξ_cal (K1): pula nie widzi obszaru użycia
    "xical_step_0.12_x1.5": {"halo_xical_step": (0.12, 1.5)},
    "xical_step_0.12_x2": {"halo_xical_step": (0.12, 2.0)},
    "xical_step_0.07_x2": {"halo_xical_step": (0.07, 2.0)},   # skok poniżej zasięgu obszaru kontrolnego
    "xical_pow_0.3": {"halo_xical_pow": 0.3},              # gładka zależność (widoczna także w puli)
    "halo_xi_tilt_low": {"halo_xi_low": (0.14, 0.9)},      # inne widmo ξ halo przy małym ξ_cal
    "xical_smear_0.5": {"xical_smear": 0.5},               # przeciek własnych protonów ponad cięcie
    # --- z v1
    "independence_broken": {"halo_scale_D": 1.6},          # halo w zdarzeniach z D ×1,6
    "afp_trigger": {"afp_trig_frac": 0.01},                # 1 % zdarzeń z D z trygera AFP (poza zakresem; G0 na danych je wyklucza)
}


def prepare(tc: ToyConfig):
    from . import afp_mix as M
    c, e, tr = generate(tc)
    C = M.standardize(c, required=M.EVENT_REQUIRED + ("m", "eta", "nsel", "hasD", "trig"))
    E = M.standardize(e, required=M.EVENT_REQUIRED + ("nsel", "hasD", "trig", "b1"))
    return C, E, tr


def run_scenario(name="nominal", donor_mode="xical_veto", seed=1, gc=None, cfg=None, **kw):
    """Toy → procedura v2 → (raport, prawda)."""
    from . import afp_mix as M
    C, E, tr = prepare(ToyConfig(seed=seed, **{**SCENARIOS[name], **kw}))
    cfg = cfg or M.MixConfig(donor_mode=donor_mode)
    rep = M.run_procedure(C, E, cfg=cfg, gc=gc)
    hd = tr[tr.hasD]
    truth = {"n_acc": int(hd.acc_region.sum()), "eps_xical_acc": float(hd.acc_region.sum() / max(hd.acc_tag_nocut.sum(), 1)),
             "fbar": float(hd.lost.sum() / max(hd.own_region.sum(), 1)), "n_diff": int(hd.own_region.sum())}
    return rep, truth


def power_table(scenarios=None, seeds=(1, 2, 3, 4, 5, 6), n_per_lb=2500, n_lb=60, gc=None, log=print):
    """Tabela mocy: per scenariusz bias N_acc/prawda (A v2, bez warstw aktywności ≈ v1, B(i), B(ii)), średnia deklarowana
    niepewność, pokrycie (|A − prawda| < 2σ_tot) i częstość FAIL każdej bramki. → (wiersze, markdown)."""
    from . import afp_mix as M
    gc = gc or M.GateConfig(n_boot=100, n_pseudo=500)
    rows = []
    for sc in (scenarios or list(SCENARIOS)):
        acc = {k: [] for k in ("A", "v1", "Bi", "Bii", "tot", "cover")}
        fails = {g: 0 for g in M.VERDICT_GATES}
        nver = 0
        for sd in seeds:
            rep, tr = run_scenario(sc, seed=sd, gc=gc, n_per_lb=n_per_lb, n_lb=n_lb)
            t = max(tr["n_acc"], 1)
            nA = rep["estimate"]["n_acc"]
            acc["A"].append(nA / t)
            acc["v1"].append(rep["info_variants"]["v1_like_LB_only"]["n_acc"] / t)
            acc["Bi"].append(rep["method_Bi"].get("n_acc", np.nan) / t)
            acc["Bii"].append(rep["method_Bii"].get("n_acc", np.nan) / t)
            tot = rep["syst_rel"]["total_rel"]
            acc["tot"].append(tot)
            acc["cover"].append(abs(nA - t) < 2 * np.hypot(tot * nA, np.sqrt(t)))
            for g in rep["gates"]:
                if g["gate"] in fails and g["status"] == "FAIL":
                    fails[g["gate"]] += 1
            nver += rep["verdict"] == "FAIL"
        n = len(seeds)
        row = {"scenario": sc, **{k: float(np.nanmean(v)) for k, v in acc.items()}, "A_sd": float(np.nanstd(acc["A"], ddof=1) / np.sqrt(n)) if n > 1 else 0.0,
               **{g: fails[g] / n for g in fails}, "verdict_fail": nver / n, "n_seeds": n}
        rows.append(row)
        if log:
            log(f"{sc}: A/prawda {row['A']:.2f} ± {row['A_sd']:.2f}, FAIL {nver}/{n}")
    short = {"G0_inputs": "G0", "G1_pool_homogeneity": "G1", "G2_disjoint_closure": "G2", "G5_cr_shape": "G5", "G7_methods_A_vs_Bii": "G7",
             "G8_methods_A_vs_Bi": "G8", "G9_opposite_side": "G9"}
    L = ["| scenariusz | A/prawda | tylko LB (≈ v1)/prawda | B(i)/prawda | B(ii)/prawda | ⟨σ_tot⟩ dekl. | pokrycie 2σ | "
         + " | ".join(short[g] for g in M.VERDICT_GATES) + " | werdykt FAIL |", "|---|---|---|---|---|---|---|" + "---|" * (len(M.VERDICT_GATES) + 1)]
    for r in rows:
        L.append(f"| {r['scenario']} | {r['A']:.2f} ± {r['A_sd']:.2f} | {r['v1']:.2f} | {r['Bi']:.2f} | {r['Bii']:.2f} | {r['tot']:.2f} | {r['cover']:.2f} | "
                 + " | ".join(f"{r[g]:.2f}" for g in M.VERDICT_GATES) + f" | {r['verdict_fail']:.2f} |")
    return rows, "\n".join(L)


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(description="AFP-MIX v2 toy: domknięcie metody mieszania protonów, scenariusze naruszające założenia")
    ap.add_argument("--scenario", default="nominal", choices=list(SCENARIOS))
    ap.add_argument("--donor-mode", default="xical_veto", choices=["xical_veto", "all"])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-per-lb", type=int, default=4000)
    ap.add_argument("--power", action="store_true", help="tabela mocy bramek (wszystkie scenariusze × --seeds)")
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--out", default="", help="--power: zapisz markdown/json pod tym prefiksem")
    a = ap.parse_args(argv)
    if a.power:
        rows, md = power_table(seeds=tuple(range(1, a.seeds + 1)))
        print(md)
        if a.out:
            open(a.out + ".md", "w").write(md + "\n")
            json.dump(rows, open(a.out + ".json", "w"), indent=1)
        return
    rep, truth = run_scenario(a.scenario, a.donor_mode, a.seed, n_per_lb=a.n_per_lb)
    est, sy = rep["estimate"], rep["syst_rel"]
    print(f"# AFP-MIX v2 toy: scenariusz={a.scenario} pula={a.donor_mode} seed={a.seed}")
    print("| wielkość | mieszanie | prawda |\n|---|---|---|")
    print(f"| N_acc A (A∨C, po cięciu ξ_cal) | {est['n_acc']:.1f} ± {est['n_acc_err_stat']:.1f} (stat) | {truth['n_acc']} |")
    for k, lab in (("method_Bi", "N_acc B(i)"), ("method_Bii", "N_acc B(ii)")):
        b = rep[k]
        if "n_acc" in b:
            print(f"| {lab} | {b['n_acc']:.1f} ± {b['n_acc_err']:.1f} | {truth['n_acc']} |")
    print(f"| N_acc bez warstw aktywności i klasy D (≈ v1) | {rep['info_variants']['v1_like_LB_only']['n_acc']:.1f} | {truth['n_acc']} |")
    print(f"| ε_ξcal(przypadkowe) | {est['eps_xical_acc']:.4f} | {truth['eps_xical_acc']:.4f} |")
    print(f"| f̄ (strata prawdziwego tagu) | {est['fbar']:.4f} | {truth['fbar']:.4f} |")
    print(f"| N_diff | {est['n_diff']:.0f} ± {est['n_diff_err']:.0f} | {truth['n_diff']} |")
    print(f"| f_acc = N_acc/N_tag | {est['f_acc']:.4f} | {truth['n_acc'] / est['n_tag']:.4f} |")
    print(f"\nbłąd względny N_acc: stat {sy['stat']:.3f}, pula {sy['pool']:.3f}, ekstrapolacja {sy['extrapolation']:.3f} (κ = {sy['kappa']:.3f} ± {sy['kappa_err']:.3f}), "
          f"połowa różnicy metod {sy['half_diff_methods']:.3f} → razem {sy['total_rel']:.3f}")
    print("warianty puli: " + ", ".join(f"{k} ×{v['ratio']:.3f}±{v['delta_err']:.3f}" for k, v in rep["variants"].items()) + "\n")
    print("| bramka | status | kluczowe |\n|---|---|---|")
    keys = ("coverage", "p", "lr", "ndf", "n_obs", "n_pred", "ratio", "ratio_err", "z", "z_slope", "f_acc_eta", "f_acc_mix")
    for g in rep["gates"]:
        kv = ", ".join(f"{k}={g[k]:.3g}" for k in keys if k in g and isinstance(g[k], (int, float, np.floating)))
        print(f"| {g['gate']} | {g['status']} | {kv} |")
    print(f"\n**Werdykt: {rep['verdict']}**" + (f" (FAIL: {', '.join(rep['verdict_failed'])})" if rep["verdict_failed"] else ""))


if __name__ == "__main__":
    main()
