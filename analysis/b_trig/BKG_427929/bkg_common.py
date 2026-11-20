"""BKG_427929 — wspólne: wczytanie all_* (inkluzywne D0, 1,75–2,00 GeV), regiony, progi WP, estymatory ε_B pod pikiem.
Tylko odczyt z ../runs/lxplus/*.npz. Definicje zgodne z btrig_analyze.cmd_epsb / charmbkg.val_epsb.default_regions('D0')."""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "runs", "lxplus")
OUT = os.path.join(HERE, "out")
M0 = 1.86484
LO, HI = 1.75, 2.00
SR_HALF, SB_MIN, SB_MAX = 0.055, 0.075, 0.125
NB = 50                                        # biny masy 5 MeV
BINS = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
WPS = ["0.95", "0.90", "0.80", "0.70", "0.50"]
MBTS_CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
               "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
VARS = ["pt", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500"]     # dostępne w all_* (bez masy i score)
FILES = {"d427": "bt_data_427929.npz", "d428": "bt_data_428770.npz", "mb": "bt_mb_nd.npz", "fd": "bt_fd_nd.npz"}


def load_all(tag, extra=()):
    """Słownik kolumn all_* (bez prefiksu). Dane 427929: GRL & HLT_mb_mbts_L1MBTS_1 (klasa MB427929 z btrig_common.classify)."""
    z = np.load(os.path.join(RUNS, FILES[tag]), allow_pickle=True)
    cols = ["run", "event", "score", "pt", "m", "truth", "ptcone40", "ptcone100", "e_nSelBL", "e_lead_pt", "e_ntrk500",
            "e_actual_mu", "e_grl"] + list(extra)
    D = {c: z["all_" + c] for c in cols}
    n0 = len(D["m"])
    if tag in ("d427", "d428"):
        tp = {c: z["all_e_tp_" + c] == 1 for c in MBTS_CHAINS + ["HLT_noalg_L1J12"]}
        D["mbts_any"] = np.any([tp[c] for c in MBTS_CHAINS], axis=0)
        D["j12"] = tp["HLT_noalg_L1J12"]
        D["j20"] = tp["HLT_j20_L1MBTS_2"]
        D["mb_mbts"] = tp["HLT_mb_mbts_L1MBTS_1"]
    if tag == "d427":
        keep = (D["e_grl"] == 1) & D["mb_mbts"]
        D = {k: v[keep] for k, v in D.items()}
    D["_n0"] = n0
    return D


def thresholds():
    """Progi WP = kwantyle score prawdziwych D0 (truth 0/1, |m−m0|<55 MeV) w fd_nd — jak btrig_analyze.thresholds."""
    z = np.load(os.path.join(RUNS, FILES["fd"]), allow_pickle=True)
    sc, pt, tr, m = z["all_score"], z["all_pt"], z["all_truth"], z["all_m"]
    sig = np.isin(tr, [0, 1]) & (np.abs(m - M0) < SR_HALF)
    thr = {}
    for b, (lo, hi) in BINS.items():
        s = sc[sig & (pt >= lo) & (pt < hi)]
        thr[b] = [float(np.quantile(s, 1 - float(wp))) for wp in WPS]
    return thr


EDGES = np.linspace(LO, HI, NB + 1)
CEN = 0.5 * (EDGES[1:] + EDGES[:-1])
T = (CEN - M0) / 0.1
SBm = np.abs(CEN - M0) >= SB_MIN               # biny SB (jak btrig: środek binu ≥ 75 MeV od m0)
SRm = np.abs(CEN - M0) < SR_HALF               # 22 biny w pełni w SR (±55 MeV)
_XF = np.linspace(M0 - SR_HALF, M0 + SR_HALF, 401)
_TF = (_XF - M0) / 0.1
_BW = EDGES[1] - EDGES[0]


def cat_hist(m, score, thr_list, w=None):
    """Histogram [kategoria score (0..5 = liczba spełnionych progów WP), bin masy]. pass(WP j) = suma kategorii ≥ j+1."""
    cat = np.searchsorted(np.asarray(thr_list), score, side="right")
    im = np.clip(((m - LO) / _BW).astype(int), 0, NB - 1)
    H = np.zeros((len(thr_list) + 1, NB))
    np.add.at(H, (cat, im), 1.0 if w is None else w)
    return H


def glm_pois(h, order, mask=SBm):
    """Poisson ML: N(bin) = exp(Σ p_k t^k) na binach `mask`; zwraca p (IRLS/Newton)."""
    X = np.vander(T[mask], order + 1, increasing=True)
    y = h[mask]
    if y.sum() <= 0:
        return None
    p = np.zeros(order + 1)
    p[0] = np.log(max(y.mean(), 1e-9))
    for _ in range(100):
        mu = np.exp(np.clip(X @ p, -50, 50))
        g = X.T @ (y - mu)
        Hs = X.T @ (X * mu[:, None])
        try:
            dp = np.linalg.solve(Hs + 1e-10 * np.eye(order + 1), g)
        except np.linalg.LinAlgError:
            return None
        st = np.abs(dp).max()
        if st > 5:
            dp *= 5 / st
        p = p + dp
        if st < 1e-9:
            break
    return p


def pred_sr(p):
    """Całka dopasowanej gęstości po SR (liczba kandydatów)."""
    if p is None:
        return np.nan
    dens = np.exp(np.clip(np.vander(_TF, len(p), increasing=True) @ p, -50, 50)) / _BW
    return float(np.sum(0.5 * (dens[1:] + dens[:-1]) * np.diff(_XF)))


def eps_expk(h_all, h_pass, order):
    """ε_B(SR) = N_pass_pred / N_all_pred, oba z ML exp(pol-k) na SB (order=2: nominał „pol2” jak btrig eps_b_pol2_ml, kształt wolny)."""
    a, b = pred_sr(glm_pois(h_all, order)), pred_sr(glm_pois(h_pass, order))
    return b / a if a > 0 else np.nan


def eps_sym(h_all, h_pass):
    """SB „symetryczne” wg val_epsb.eps_b_sb_sym: |m−m0| ∈ [75, 125] MeV w oknie 1,75–2,00 (lewe pasmo ucięte do 40 MeV)."""
    s = SBm & (np.abs(CEN - M0) <= SB_MAX)
    return h_pass[s].sum() / max(h_all[s].sum(), 1e-12)


def eps_lr(h_all, h_pass):
    """Średnia ε lewego i prawego pasma interpolowana liniowo do m0 (wariant niezależny od kształtu widma)."""
    L, R = SBm & (CEN < M0), SBm & (CEN > M0)
    eL, eR = h_pass[L].sum() / max(h_all[L].sum(), 1e-12), h_pass[R].sum() / max(h_all[R].sum(), 1e-12)
    mL, mR = np.average(CEN[L], weights=h_all[L] + 1e-12), np.average(CEN[R], weights=h_all[R] + 1e-12)
    return eL + (eR - eL) * (M0 - mL) / (mR - mL)


def eps_true_sr(h_all, h_pass):
    return h_pass[SRm].sum() / max(h_all[SRm].sum(), 1e-12)


METHODS = {"pol2": lambda a, p: eps_expk(a, p, 2), "pol1": lambda a, p: eps_expk(a, p, 1), "pol3": lambda a, p: eps_expk(a, p, 3),
           "sym": eps_sym, "lr": eps_lr}


def all_eps(H, method):
    """H = cat_hist → wektor ε po WP."""
    ha = H.sum(0)
    f = METHODS[method] if method != "true" else eps_true_sr
    return np.array([f(ha, H[j + 1:].sum(0)) for j in range(H.shape[0] - 1)])


def boot(H, rng):
    return rng.poisson(H).astype(float)


def atlas_label(ax, extra="", x=0.04, y=0.95, fs=11):
    t = ax.text(x, y, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", fontsize=fs, va="top")
    ax.annotate("", xy=(1, 0), xycoords=t, xytext=(4, 0), textcoords="offset points", fontsize=fs, va="bottom")
    if extra:
        ax.text(x, y - 0.075, extra, transform=ax.transAxes, fontsize=fs - 3, va="top")
