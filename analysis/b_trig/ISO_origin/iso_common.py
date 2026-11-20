"""ISO_origin — wspólne: wczytanie, klasy, estymatory ε_iso (fit Δm pass/fail i odejmowanie tła z pasm Δm)."""
import os, sys, time, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
BT = os.path.dirname(HERE)
sys.path.insert(0, BT)
import btrig_common as C  # noqa
C.add_paths()
from charmtag import dstar_refit as DR  # noqa
R = os.path.join(BT, "runs", "lxplus")
E_DM = C.E_DM
PAIR_COLS = ["m", "pt", "score", "dm", "ptcone40", "ptcone100", "pt_ds", "pis_pt", "eta", "truth", "truth_ds", "event",
             "e_nSelBL", "e_lead_pt", "e_actual_mu", "e_grl", "e_npu", "e_lb", "e_ntrk500"]


def load_pairs(name):
    z = np.load(os.path.join(R, name + ".npz"), allow_pickle=False)
    d = {k: z[k] for k in PAIR_COLS if k in z.files}
    for k in z.files:
        if k.startswith("e_tp_"):
            d[k[2:]] = z[k]
    d["rel40"] = d["ptcone40"] / d["pt"]
    d["z"] = d["pt"] / np.maximum(d["e_lead_pt"], 1e-3)
    d["ann"] = d["ptcone100"] - d["ptcone40"]
    return d


def classes(D, run):
    """Klasy rozłączne i diagnostyczne (maski na parach)."""
    tp = lambda c: D.get("tp_" + c, np.zeros(len(D["m"]), np.int8)) == 1  # noqa
    if run == 427929:
        g = D["e_grl"] == 1
        lowmu = D["e_actual_mu"] < 0.02
        return {"MB427929": tp("HLT_mb_mbts_L1MBTS_1") & g, "MB9_J12": tp("HLT_mb_mbts_L1MBTS_1") & g & tp("HLT_noalg_L1J12"),
                "MB9_notJ12": tp("HLT_mb_mbts_L1MBTS_1") & g & ~tp("HLT_noalg_L1J12"),
                "MB9_j20": tp("HLT_mb_mbts_L1MBTS_1") & g & tp("HLT_j20_L1MBTS_2"),
                "MB9_sptrk4": tp("HLT_mb_mbts_L1MBTS_1") & g & tp("HLT_mb_sptrk_pt4_L1MBTS_2"),
                "MB9_sptrk4_notJ12": tp("HLT_mb_mbts_L1MBTS_1") & g & tp("HLT_mb_sptrk_pt4_L1MBTS_2") & ~tp("HLT_noalg_L1J12"),
                "nonGRL": ~g, "_lowmu": lowmu}
    cl = C.classify({k: v for k, v in D.items() if k.startswith("tp_")}, 428770)
    j12 = tp("HLT_noalg_L1J12")
    unb = tp("HLT_noalg_mb_L1MBTS_1") | tp("HLT_noalg_mb_L1MBTS_2") | tp("HLT_mb_mbts_L1MBTS_1")
    hard = tp("HLT_j20_L1MBTS_2") | tp("HLT_mb_sptrk_pt4_L1MBTS_2")
    out = {"all": cl["all"], "J12only": cl["J12only"], "MBTSbased": cl["MBTSbased"], "other": cl["other"],
           "J12any": j12, "notJ12": ~j12, "MBTS_notJ12": cl["MBTSbased"] & ~j12, "MBTS_andJ12": cl["MBTSbased"] & j12,
           "sptrk4": cl["sptrk4"], "sptrk4_notJ12": cl["sptrk4"] & ~j12, "j20MBTS": cl["j20MBTS"], "j20_notJ12": cl["j20MBTS"] & ~j12,
           "unbMBTS": unb, "unbMBTS_or_pt2": unb | tp("HLT_mb_sptrk_pt2_L1MBTS_2"),
           "MBTSsoft_only": (unb | tp("HLT_mb_sptrk_pt2_L1MBTS_2")) & ~hard & ~j12,
           "AFPsptrk4": tp("HLT_mb_sptrk_pt4_L1AFP_A_OR_C") & ~j12 & ~cl["MBTSbased"]}
    return out


def hist(dm):
    return np.histogram(np.asarray(dm) * 1e3, E_DM)[0]


def fit_pf(dm_pass, dm_fail, shape, interval=False):
    f = C.fit_eps(dm_pass, dm_fail, shape, interval=interval, model_check=False)
    err = C.eps_sym_err(f) if interval else f.get("eps_err", np.nan)
    return {"eps": f["eps"], "err": err, "N": f["N"], "N_err": f["N_err"], "det": bool(f["determined"]), "valid": bool(f["valid"]),
            "n_pass": f["n_pass"], "n_fail": f["n_fail"], "method": "fit"}


# --- odejmowanie tła: okno sygnału Δm 144,4–146,6; tło z fitu progowego do pasm (140,5–143,4 ∪ 148,5–165) ---------------
SIG = (144.4, 146.6)
MPI = 139.57


def _bkg_yield(h, e):
    """Tło w oknie SIG z fitu u^a·exp(-b·u)·N do pasm; zwraca (B, σ_B). Liniowy fit ważony w log (Poisson ML iteracyjnie)."""
    from scipy.optimize import minimize
    c = 0.5 * (e[1:] + e[:-1])
    u = c - MPI
    side = ((c > 140.6) & (c < 143.4)) | ((c > 148.5) & (c < 165.0))
    win = (c > SIG[0]) & (c < SIG[1])
    hs = h[side].astype(float)
    if hs.sum() < 5:
        return 0.0, 1.0, float(h[win].sum())
    X = np.vstack([np.ones(side.sum()), np.log(u[side]), u[side]]).T

    def nll(p):
        lam = np.exp(np.clip(X @ p, -30, 30))
        return float(np.sum(lam - hs * np.log(lam)))

    def grad(p):
        lam = np.exp(np.clip(X @ p, -30, 30))
        return X.T @ (lam - hs)
    p0 = np.array([np.log(max(hs.mean(), 0.1)), 0.5, 0.0])
    r = minimize(nll, p0, jac=grad, method="BFGS")
    lam_w = np.exp(np.vstack([np.ones(win.sum()), np.log(u[win]), u[win]]).T @ r.x)
    B = float(lam_w.sum())
    # błąd B: propagacja przez kowariancję (inv Hessian)
    lam = np.exp(np.clip(X @ r.x, -30, 30))
    H = X.T @ (X * lam[:, None])
    try:
        cov = np.linalg.inv(H)
        Xw = np.vstack([np.ones(win.sum()), np.log(u[win]), u[win]]).T
        g = Xw.T @ lam_w
        vB = float(g @ cov @ g)
    except np.linalg.LinAlgError:
        vB = B
    return B, np.sqrt(max(vB, 0)), float(h[win].sum())


def sbs(dm_pass, dm_fail, w_pass=None, w_fail=None):
    """ε = S_pass/(S_pass+S_fail), S = N(okno) − B(fit pasm). Błąd: Poisson okna ⊕ błąd B, propagowany."""
    hp, hf = hist(dm_pass), hist(dm_fail)
    Bp, eBp, Np = _bkg_yield(hp, E_DM)
    Bf, eBf, Nf = _bkg_yield(hf, E_DM)
    Sp, Sf = Np - Bp, Nf - Bf
    vp, vf = Np + eBp ** 2, Nf + eBf ** 2
    T = Sp + Sf
    if T <= 0:
        return {"eps": np.nan, "err": np.nan, "N": T, "N_err": np.sqrt(vp + vf), "det": False, "method": "sbs", "Sp": Sp, "Sf": Sf}
    eps = Sp / T
    err = np.sqrt(Sf ** 2 * vp + Sp ** 2 * vf) / T ** 2
    return {"eps": float(eps), "err": float(err), "N": float(T), "N_err": float(np.sqrt(vp + vf)), "det": bool(T > 3 * np.sqrt(vp + vf)),
            "method": "sbs", "Sp": float(Sp), "Sf": float(Sf), "vp": float(vp), "vf": float(vf)}


def sig_yield(dm):
    """Uzysk sygnału (okno − tło) i błąd — do rozkładów sygnału w danych."""
    h = hist(dm)
    B, eB, N = _bkg_yield(h, E_DM)
    return N - B, np.sqrt(N + eB ** 2)


def chi2_const(v, e):
    v, e = np.asarray(v, float), np.asarray(e, float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    return C.chi2_consistency(v[ok], e[ok])


# --- wagi odejmowania tła per kandydat (rozkłady sygnału / ε różniczkowe / przeważanie MC) ------------------------------
PT_PARENT = np.array([3.0, 4.0, 5.0, 8.0, 12.0, 1e9])
_cc = 0.5 * (E_DM[1:] + E_DM[:-1])
_SIDE = ((_cc > 140.6) & (_cc < 143.4)) | ((_cc > 148.5) & (_cc < 165.0))
_WIN = (_cc > SIG[0]) & (_cc < SIG[1])


def sbs_w(dm, pt, pas):
    """Waga per kandydat: +1 w oknie sygnału Δm, −r w pasmach, 0 poza; r = B_okno/N_pasm z fitu progowego w grupie
    (bin pT rodzica × pass/fail). Zakłada kształt tła Δm niezależny od pozostałych zmiennych wewnątrz grupy."""
    dm = np.asarray(dm) * 1e3
    ib = np.clip(np.digitize(dm, E_DM) - 1, 0, len(_cc) - 1)
    inr = (dm >= E_DM[0]) & (dm < E_DM[-1])
    win, side = _WIN[ib] & inr, _SIDE[ib] & inr
    w = np.zeros(len(dm))
    w[win] = 1.0
    ip = np.clip(np.digitize(pt, PT_PARENT) - 1, 0, len(PT_PARENT) - 2)
    for ps in (True, False):
        g0 = (pas == ps)
        h0 = np.histogram(dm[g0], E_DM)[0]
        B0, _, _ = _bkg_yield(h0, E_DM)
        r0 = B0 / max(h0[_SIDE].sum(), 1)
        for k in range(len(PT_PARENT) - 1):
            g = g0 & (ip == k)
            if not g.any():
                continue
            h = np.histogram(dm[g], E_DM)[0]
            if h[_SIDE].sum() >= 60:
                B, _, _ = _bkg_yield(h, E_DM)
                r = B / h[_SIDE].sum()
            else:
                r = r0
            w[g & side] = -r
    return w


def eff_w(w, pas):
    """ε = Σw·pass/Σw i błąd (metoda delta, wagi niezależne)."""
    T = w.sum()
    if T <= 0:
        return np.nan, np.nan, float(T), float(np.sqrt((w ** 2).sum()))
    e = (w * pas).sum() / T
    v = (w ** 2 * (pas - e) ** 2).sum() / T ** 2
    return float(e), float(np.sqrt(v)), float(T), float(np.sqrt((w ** 2).sum()))
