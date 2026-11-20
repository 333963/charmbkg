"""FIT_scorebins — Asimov (optymalizacja granic, DP), pseudo-eksperymenty, odporność (SF, rzeźbienie).

Użycie: fsb_run.py [NTOY=600] [NPROC=4]   →  out/results_<gatunek>_<ib>.json, out/toys_<gatunek>_<ib>.npz
"""
import os, sys, json, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fsb_common import EFF_EDGES, NCELL, WPS, PTB, merge, Model, build

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
KS = [2, 3, 5, 8]
SF_SCEN = {"A_428770_5-8": [(0.5, 0.96), (0.8, 0.87), (0.95, 0.84)],
           "B_428770_8-30": [(0.5, 0.94), (0.8, 0.97), (0.95, 1.01)],
           "C_427929_5-8": [(0.5, 0.99), (0.8, 1.03), (0.95, 1.05)]}
SCULPT = {"R089": 0.11, "R097": 0.03}


def load(sp, ib):
    z = np.load(f"{OUT}/inputs_{sp}.npz")
    return dict(f=z[f"f_{ib}"], Us=z[f"Us_{ib}"], Ur=z[f"Ur_{ib}"], B=z[f"B_{ib}"], NS=float(z[f"NS_{ib}"]), X=z["X"],
                edges=z["edges"])


def interval_var(inp, NS):
    """V[s, e]: wariancja Asimova amplitudy sygnału w kategorii = komórki s..e-1 (tło exp(pol2) swobodne)."""
    T = inp["f"][:, None] * (inp["Us"] + inp["Ur"])
    V = np.full((NCELL + 1, NCELL + 1), np.inf)
    F = np.zeros_like(V)
    for s in range(NCELL):
        for e in range(s + 1, NCELL + 1):
            fc = inp["f"][s:e].sum()
            U = T[s:e].sum(0) / fc
            b = inp["B"][s:e].sum(0)
            mu = NS * fc * U + b
            J = np.c_[U, b[:, None] * inp["X"]]
            H = (J / mu[:, None]).T @ J
            V[s, e] = np.linalg.inv(H)[0, 0]
            F[s, e] = fc
    return V, F


def dp_best(V, F, K, s0=0):
    """Podział [s0, NCELL) na K przedziałów maksymalizujący sumę f_c^2 / V_c (informacja o N_S przy ułamkach ustalonych)."""
    I = np.where(np.isfinite(V), F ** 2 / V, -np.inf)
    best = {(0, s0): (0.0, [])}
    for k in range(1, K + 1):
        for e in range(s0 + k, NCELL + 1):
            cand = [(best[(k - 1, s)][0] + I[s, e], best[(k - 1, s)][1] + [s]) for s in range(s0 + k - 1, e) if (k - 1, s) in best]
            if cand:
                best[(k, e)] = max(cand, key=lambda t: t[0])
    info, starts = best[(K, NCELL)]
    return starts, 1 / np.sqrt(info)


def sf_fractions(f, pts):
    """f' komórek, gdy skumulowana wydajność E'(E) = E * SF(E); SF liniowo między punktami, stałe poniżej, → 1 przy E = 1."""
    xs = [p[0] for p in pts] + [1.0]
    ys = [p[1] for p in pts] + [1.0]
    E = EFF_EDGES.copy()                       # E[j] = wydajność MC powyżej dolnej granicy komórki j; E[NCELL] = 0
    Ep = E * np.interp(E, xs, ys)
    return Ep[:-1] - Ep[1:]


def configs(inp):
    V, F = interval_var(inp, inp["NS"])
    cfg = {}
    cfg["cut100"] = dict(bounds=[0], mode="fixed")
    for w in WPS:
        k = int(np.argmin(np.abs(EFF_EDGES - w)))
        cfg[f"cut{int(round(w * 100))}"] = dict(bounds=[k], mode="fixed")
    asim = {n: float(np.sqrt(V[c["bounds"][0], NCELL]) / F[c["bounds"][0], NCELL]) for n, c in cfg.items()}
    for K in KS:
        st, sg = dp_best(V, F, K)
        cfg[f"fix{K}"] = dict(bounds=st, mode="fixed"); asim[f"fix{K}"] = float(sg)
    cfg["fix17"] = dict(bounds=list(range(NCELL)), mode="fixed")
    asim["fix17"] = float(1 / np.sqrt(sum(F[s, s + 1] ** 2 / V[s, s + 1] for s in range(NCELL))))
    st, sg = dp_best(V, F, 3, s0=1)
    cfg["fix3_t95"] = dict(bounds=st, mode="fixed"); asim["fix3_t95"] = float(sg)
    for K in (2, 3, 5):
        st = cfg[f"fix{K}"]["bounds"]
        cfg[f"free{K}"] = dict(bounds=st, mode="free")
        asim[f"free{K}"] = float(np.sqrt(sum(V[a, b] for a, b in zip(st, st[1:] + [NCELL]))))
        for sg_, tag in ((0.05, "c05"), (0.10, "c10")):
            cfg[f"{tag}_{K}"] = dict(bounds=st, mode="constr", sig=sg_)
    return cfg, asim, V, F


def truth(inp, f=None, sculpt=0.0, NS=None):
    f = inp["f"] if f is None else f
    NS = inp["NS"] if NS is None else NS
    S = NS * f[:, None] * (inp["Us"] + inp["Ur"])
    B = inp["B"].copy()
    if sculpt:
        cen = 0.5 * (inp["edges"][1:] + inp["edges"][:-1])
        m0 = cen[np.argmax(inp["Us"].sum(0))]
        g = np.exp(-0.5 * ((cen - m0) / 0.055) ** 2)
        emid = 0.5 * (EFF_EDGES[:-1] + EFF_EDGES[1:])
        amp = sculpt * np.minimum(1.0, (1 - emid) / 0.3)          # pełna amplituda powyżej WP70, zero przy score → 0
        B = B * (1 - amp[:, None] * g[None, :])
    return S, B


def fit_cfg(inp, c, n_cells, B_cells):
    m, fc = build(inp, c["bounds"], c["mode"], sig=c.get("sig"))
    n = merge(n_cells, c["bounds"])
    th0 = m.theta0(inp["NS"], merge(B_cells, c["bounds"]))
    th, N, e, ok, _ = m.fit(n, th0)
    if c["bounds"][0] > 0 and c["mode"] == "fixed" and len(c["bounds"]) == 1:
        pass
    return N, e, ok


def run(args):
    sp, ib, ntoy, seed = args
    t0 = time.time()
    inp = load(sp, ib)
    cfg, asim, V, F = configs(inp)
    names = list(cfg)
    S, B = truth(inp)
    res = {"sp": sp, "pt": PTB[ib], "NS": inp["NS"], "bounds": {n: [float(EFF_EDGES[b]) for b in c["bounds"]] for n, c in cfg.items()},
           "asimov_rel": {n: v / inp["NS"] for n, v in asim.items()}}
    # --- Asimov: nominał i scenariusze (obciążenie bez fluktuacji) ---
    def asimov_set(S_, B_, which):
        o = {}
        for n in which:
            N, e, ok = fit_cfg(inp, cfg[n], S_ + B_, B)
            o[n] = [float(N / inp["NS"] - 1), float(e / inp["NS"]), bool(ok)]
        return o
    res["asimov_fit"] = asimov_set(S, B, names)
    best_wp = min((n for n in names if n.startswith("cut")), key=lambda n: res["asimov_fit"][n][1])
    res["best_wp"] = best_wp
    rob = {}
    for sc, pts in SF_SCEN.items():
        fp = sf_fractions(inp["f"], pts)
        S2, B2 = truth(inp, f=fp)
        rob["SF_" + sc] = asimov_set(S2, B2, names)
        # ułamki skorygowane o SF (fit zna prawdziwe f): kontrola = brak obciążenia
    for sc, a in SCULPT.items():
        S2, B2 = truth(inp, sculpt=a)
        rob["sculpt_" + sc] = asimov_set(S2, B2, names)
    res["robust"] = rob
    # --- skan normalizacji sygnału (Asimov, DP) ---
    scan = {}
    for k in (0.5, 1.0, 2.0):
        V2, F2 = interval_var(inp, inp["NS"] * k)
        cuts = {w: np.sqrt(V2[int(np.argmin(np.abs(EFF_EDGES - w))), NCELL]) / w for w in WPS}
        bw = min(cuts, key=cuts.get)
        scan[str(k)] = {"best_wp": bw, "cut": float(cuts[bw] / (inp["NS"] * k)),
                        **{f"fix{K}": float(dp_best(V2, F2, K)[1] / (inp["NS"] * k)) for K in KS}}
    res["scan_NS"] = scan
    # --- pseudo-eksperymenty (wspólne toye dla wszystkich konfiguracji → sparowane porównanie) ---
    rng = np.random.default_rng(seed)
    Nh = np.full((ntoy, len(names)), np.nan); Eh = np.full_like(Nh, np.nan); OK = np.zeros_like(Nh, bool)
    lam = S + B
    for t in range(ntoy):
        n_cells = rng.poisson(lam).astype(float)
        for i, n in enumerate(names):
            Nh[t, i], Eh[t, i], OK[t, i] = fit_cfg(inp, cfg[n], n_cells, B)
    np.savez_compressed(f"{OUT}/toys_{sp}_{ib}.npz", N=Nh, E=Eh, OK=OK, names=np.array(names), NS=inp["NS"])
    res["sec"] = time.time() - t0
    json.dump(res, open(f"{OUT}/results_{sp}_{ib}.json", "w"), indent=1)
    return sp, ib, res["sec"], best_wp


if __name__ == "__main__":
    ntoy = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else ["D0", "Dplus"]
    tasks = [(sp, ib, ntoy, 1000 + 10 * i + ib) for i, sp in enumerate(["D0", "Dplus"]) if sp in only for ib in range(3)]
    with Pool(nproc) as p:
        for r in p.imap_unordered(run, tasks):
            print(r, flush=True)
