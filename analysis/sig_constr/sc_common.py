"""SIG_constr — wspólne: wczytanie par D0←D* (bt_*.npz + wybrane kolumny tj_*.npz), klasy triggera, foldy po zdarzeniach,
wagi sygnału z danych (SBS z ISO_origin oraz sWeights z fitu Δm), estymatory ε. Projekt B. Tylko odczyt wejść."""
import os, sys, json, zlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "B_trig", "ISO_origin"))
sys.path.insert(0, os.path.join(CW, "B_trig"))
import iso_common as IC  # noqa: E402  (sbs_w, sbs, fit_pf, classes, C, DR)
C, DR = IC.C, IC.DR
OUT = os.path.join(HERE, "out")
F18 = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0", "pointing", "abs_cos_theta_star",
       "abs_impact_sig", "d_min_abs_d0", "d_max_abs_z0", "d_max_abs_eta", "dedx_K", "log_p_K", "npts_K", "dedx_pi_max", "log_p_pi_min"]
ISO = ["ptcone40", "ptcone100", "rel40", "rel100"]
EVT = ["e_nSelBL", "e_lead_pt", "e_ntrk500"]
TJ = ["o_ptmax08", "h_s_away", "o_drn5", "a_jD_pt", "z"]          # z TRKJET (O = bez torów stożka; a_jD_pt = [koło])
BASE = ["m", "pt", "eta", "phi", "score", "dm", "truth", "truth_ds", "event", "file_idx", "run", "pis_pt", "e_grl", "e_actual_mu"]
PTB = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
MWIN = 0.040
E_DM = C.E_DM
_cc = 0.5 * (E_DM[1:] + E_DM[:-1])
SIDE_HI = (148.5, 165.0)


def load(name, tj=True):
    z = np.load(os.path.join(CW, "B_trig", "runs", "lxplus", f"bt_{name}.npz"))
    d = {k: z[k] for k in BASE + F18 + ["ptcone40", "ptcone100"] + EVT}
    for k in z.files:
        if k.startswith("e_tp_"):
            d[k[2:]] = z[k]
    d["rel40"] = d["ptcone40"] / d["pt"]; d["rel100"] = d["ptcone100"] / d["pt"]
    d["ann"] = d["ptcone100"] - d["ptcone40"]
    if tj:
        p = os.path.join(CW, "B_trig", "TRKJET", "runs", "lxplus", f"tj_{name}.npz")
        t = np.load(p)
        assert len(t["m"]) == len(d["m"]) and np.array_equal(t["event"], d["event"]) and np.allclose(t["m"], d["m"]), "tj != bt"
        for k in TJ:
            d["tj_" + k if k == "z" else k] = t[k]
    # fold po zdarzeniu (deterministyczny hasz run, plik, event)
    key = (d["run"].astype(np.uint64) * np.uint64(1000003) + d["file_idx"].astype(np.uint64)) * np.uint64(2654435761) + d["event"].astype(np.uint64)
    key ^= key >> np.uint64(33); key *= np.uint64(0xff51afd7ed558ccd); key ^= key >> np.uint64(33)
    d["evkey"] = key
    d["kf"] = (key % np.uint64(5)).astype(np.int8)
    return d


def data_classes(D, run):
    cl = IC.classes(D, run)
    n = len(D["m"])
    tp = lambda c: D.get("tp_" + c, np.zeros(n, np.int8)) == 1  # noqa: E731
    if run == 427929:
        mb = cl["MB427929"]
        return {"9:MB": mb, "9:MB&j20&!J12": mb & tp("HLT_j20_L1MBTS_2") & ~tp("HLT_noalg_L1J12"),
                "9:MB&!j20&!J12": mb & ~tp("HLT_j20_L1MBTS_2") & ~tp("HLT_noalg_L1J12")}
    return {"8:all": cl["all"], "8:J12only": cl["J12only"], "8:j20&!J12": cl["j20_notJ12"], "8:MBTSbased": cl["MBTSbased"], "8:other": cl["other"]}


def base_sel(D, sthr=None, lo=3.0, hi=30.0, mwin=MWIN):
    s = (np.abs(D["m"] - C.M0_D0) < mwin) & (D["pt"] >= lo) & (D["pt"] < hi) & (D["dm"] * 1e3 >= E_DM[0]) & (D["dm"] * 1e3 < E_DM[-1])
    if sthr is not None:
        s &= D["score"] > sthr
    return s


def mc_sig(M, sthr=None, lo=3.0, hi=30.0, mwin=MWIN):
    """Prawdziwe D0 z prawdziwego D* (fd_nd), okno Δm 143,5–147,5 jak w ISO_origin."""
    s = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < mwin) & (M["pt"] >= lo) & (M["pt"] < hi)
    s &= (M["dm"] > 0.1435) & (M["dm"] < 0.1475)
    if sthr is not None:
        s &= M["score"] > sthr
    return s


def side_hi(D):
    dm = D["dm"] * 1e3
    return (dm > SIDE_HI[0]) & (dm < SIDE_HI[1])


def neff(w):
    w = np.asarray(w, float)
    s2 = (w ** 2).sum()
    return float(w.sum() ** 2 / s2) if s2 > 0 else 0.0


# ------------------------------------------------------------------------------------------------------------ sWeights
def _sig_pdf_bins(mu, sc, k, f, shift, scale):
    from scipy.stats import norm
    m = mu + shift
    s1, s2 = sc * scale, sc * k * scale
    cdf = f * norm.cdf(E_DM, m, s1) + (1 - f) * norm.cdf(E_DM, m, s2)
    return np.diff(cdf)


def sweights(dm_gev, shape, lo=140.6, free_shape=True):
    """sWeights (sPlot) z binowanego fitu ML Δm: N_s·S + N_b·u^a·e^{-bu}; S = 2 gaussy z MC (wolne przesunięcie i skala).
    Zwraca (w per kandydat [0 poza zakresem fitu], dict z N_s, N_b, błędem N_s, n_eff). Wagi liczone po ustaleniu kształtów (klasyczny sPlot)."""
    from scipy.optimize import minimize
    dm = np.asarray(dm_gev) * 1e3
    h = np.histogram(dm, E_DM)[0].astype(float)
    use = _cc > lo
    u = _cc - IC.MPI

    def comps(p):
        ns, nb, a, b, sh, scl = p
        S = _sig_pdf_bins(shape["mu"], shape["sc"], shape["k"], shape["f"], sh, scl) * use
        S = S / max(S.sum(), 1e-12)
        B = np.where(use, u ** a * np.exp(-b * u), 0.0)
        B = B / max(B.sum(), 1e-12)
        return S, B

    def nll(ns, nb, a, b, sh, scl):
        S, B = comps((ns, nb, a, b, sh, scl))
        lam = np.maximum((ns * S + nb * B)[use], 1e-9)
        return float(2 * np.sum(lam - h[use] * np.log(lam)))
    ntot = h[use].sum()
    if ntot < 30:
        return np.zeros(len(dm)), {"Ns": 0.0, "Nb": float(ntot), "Ns_err": np.nan, "neff": 0.0, "ok": False}
    from iminuit import Minuit
    win = use & (_cc > 144.4) & (_cc < 146.6)
    side = use & (((_cc > 140.6) & (_cc < 143.4)) | ((_cc > 148.5) & (_cc < 165.0)))
    B0, _, N0 = IC._bkg_yield(h, E_DM)
    ns0 = max((N0 - B0) / 0.82, 0.02 * N0, 5.0)
    best = None
    for a0, b0 in ((0.5, 0.02),):
        m = Minuit(nll, ns=ns0, nb=max(ntot - ns0, 1.0), a=a0, b=b0, sh=-0.1, scl=1.1)
        m.errordef = 1.0
        m.limits["ns"] = (0, ntot); m.limits["nb"] = (0, 1.5 * ntot); m.limits["a"] = (0.0, 4.0); m.limits["b"] = (-0.3, 1.0)
        m.limits["sh"] = (-0.5, 0.3); m.limits["scl"] = (0.8, 1.6)
        if not free_shape:
            m.fixed["scl"] = True; m.values["scl"] = 1.0
        m.migrad(ncall=20000)
        if best is None or m.fval < best.fval:
            best = m
    p = np.array([best.values[k] for k in ("ns", "nb", "a", "b", "sh", "scl")])
    S, B = comps(p)
    ns, nb = p[0], p[1]
    den = np.maximum(ns * S + nb * B, 1e-12)
    # macierz sPlot na danych binowanych
    Wss = (h * S * S / den ** 2)[use].sum(); Wsb = (h * S * B / den ** 2)[use].sum(); Wbb = (h * B * B / den ** 2)[use].sum()
    V = np.linalg.inv(np.array([[Wss, Wsb], [Wsb, Wbb]]))
    wbin = np.where(use, (V[0, 0] * S + V[0, 1] * B) / den, 0.0)
    ib = np.clip(np.digitize(dm, E_DM) - 1, 0, len(_cc) - 1)
    inr = (dm >= E_DM[0]) & (dm < E_DM[-1])
    w = np.where(inr, wbin[ib], 0.0)
    return w, {"Ns": float(ns), "Nb": float(nb), "Ns_err": float(np.sqrt(max(V[0, 0], 0))), "neff": neff(w), "shift": float(p[4]), "scale": float(p[5]), "valid": bool(best.valid),
               "sumw": float(w.sum()), "ok": bool(ns > 0)}


def mc_shapes(M):
    sh = {}
    f = np.isin(M["truth"], [0, 1]) & np.isin(M["truth_ds"], [0, 1]) & (np.abs(M["m"] - C.M0_D0) < MWIN)
    for b, (lo, hi) in list(PTB.items()) + [("3-30", (3.0, 30.0))]:
        sh[b] = DR.shape_from_mc(M["dm"][f & (M["pt"] >= lo) & (M["pt"] < hi)] * 1e3)
    return sh


def sw_per_ptbin(D, sel, shapes):
    """sWeights liczone osobno w binach pT (3–5, 5–8, 8–30); zwraca w (dla sel) i słownik fitów."""
    w = np.zeros(int(sel.sum())); info = {}
    pt = D["pt"][sel]; dm = D["dm"][sel]
    for b, (lo, hi) in PTB.items():
        g = (pt >= lo) & (pt < hi)
        if g.sum() < 30:
            info[b] = {"Ns": 0.0, "neff": 0.0, "ok": False}; continue
        w[g], info[b] = sweights(dm[g], shapes[b])
    return w, info


def sbs_weights(D, sel, pas=None):
    pas = np.ones(int(sel.sum()), bool) if pas is None else pas
    return IC.sbs_w(D["dm"][sel], D["pt"][sel], pas)


def weff(w, pas):
    return IC.eff_w(w, pas.astype(float))


def boot_ev(evkey, nb=200, seed=1):
    """Wagi bootstrapu po zdarzeniach: generator Poissona(1) per zdarzenie → per kandydat."""
    rng = np.random.default_rng(seed)
    u, inv = np.unique(evkey, return_inverse=True)
    for _ in range(nb):
        yield rng.poisson(1.0, len(u))[inv].astype(float)


def jdump(obj, name):
    def cv(o):
        if isinstance(o, dict):
            return {str(k): cv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [cv(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return o
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(cv(obj), f, indent=1)
