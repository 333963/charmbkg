"""HIGHSCORE — analiza EKSPLORACYJNA PO FAKCIE (18.11.2026, projekt B): deficyt sygnału D0<-D* przy score H > WP25, pT 8–30 GeV.
Wspólne: wczytanie bt_*.npz (tylko odczyt), estymator wydajków zliczeniem (okno Δm − tło z fitu pasm bocznych, kształt tła z R3),
bootstrap po zdarzeniach. Nic nie zapisuje poza własnym katalogiem. Dane inkluzywne (bez tagu AFP)."""
import os, sys, json
import numpy as np
from iminuit import Minuit

HERE = os.path.dirname(os.path.abspath(__file__))
R3 = os.path.join(os.path.dirname(HERE), "R3_2026-11-18")
sys.path.insert(0, R3)
import r3_common as R  # noqa: E402  (silnik R3, tylko import)

OUT = os.path.join(HERE, "out")
MU = {b: R.official()["shapes"][b]["mu"] for b in R.BINS}
THR = R.official()["thresholds"]
FEATS = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0", "pointing", "abs_cos_theta_star",
         "abs_impact_sig", "d_min_abs_d0", "d_max_abs_z0", "d_max_abs_eta", "dedx_K", "log_p_K", "npts_K", "dedx_pi_max", "log_p_pi_min"]
EXTRA = ["m", "pis_pt", "e_nSelBL", "e_ntrk500", "ptcone40", "ptcone100", "e_lead_pt"]
MBTS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
        "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
BASE = ["run", "event", "m", "pt", "eta", "score", "dm", "truth", "truth_ds", "e_grl", "e_npu", "e_lb"] + FEATS + EXTRA
_cache = {}


def load(name):
    if name in _cache:
        return _cache[name]
    z = np.load(os.path.join(R.BT, f"bt_{name}.npz"), allow_pickle=False)
    d = {k: z[k] for k in dict.fromkeys(BASE)}
    d["dm"] = d["dm"].astype(float) * 1e3
    d["score"] = d["score"].astype(float)
    if name.startswith("data"):
        mb = np.zeros(len(d["m"]), bool)
        for c in MBTS:
            mb |= z[f"e_tp_{c}"] == 1
        j12 = z["e_tp_HLT_noalg_L1J12"] == 1
        d["cls_J12only"] = j12 & ~mb
        d["cls_MBTSbased"] = mb
        d["cls_other"] = ~j12 & ~mb
        _, d["inv"] = np.unique(d["event"], return_inverse=True)
        d["nev"] = int(d["inv"].max()) + 1
    else:
        d["sig"] = np.isin(d["truth"], [0, 1]) & np.isin(d["truth_ds"], [0, 1])
        d["swap"] = np.isin(d["truth"], [100, 101]) & np.isin(d["truth_ds"], [0, 1])
        d["np"] = d["truth"] == 1          # interpretacja: 1 = nieprompt (z b) — wniosek z rozkładów IP/Lxy, nie z dokumentacji
    _cache[name] = d
    return d


def base_mask(d, b, mwin=R.NOM_WIN, isdata=True, run=None):
    lo, hi = b if isinstance(b, tuple) else R.BINS[b]
    s = (np.abs(d["m"] - R.M0_D0) < mwin) & (d["pt"] >= lo) & (d["pt"] < hi)
    if isdata and run == 427929:
        s &= d["e_grl"] == 1
    if not isdata:
        s &= d["sig"]
    return s


def thr_for(b, wps):
    """Progi score dla skumulowanych wydajności wps (malejąco); oficjalne gdzie są, inaczej kwantyl par MC (jak R3.edges_for)."""
    M = load("fd_nd")
    s = base_mask(M, b, isdata=False)
    return np.array([THR[b][f"{w:.2f}"] if f"{w:.2f}" in THR[b] else float(np.quantile(M["score"][s], 1 - w)) for w in wps])


def fit_bkg(h, sb, nco, x0=None):
    """Fit tła u^a(1+c1 t+c2 t^2) do binów sb (Poisson). Zwraca (B, a, c1, c2)."""
    hs = h[sb]; pos = hs > 0

    def f(B, a, c1, c2):
        L = np.clip(B * R.bkg_shape(a, c1, c2)[sb], 1e-9, None)
        d = np.where(pos, L / np.where(pos, hs, 1.0) - 1.0, 0.0)
        return float(np.sum(np.where(pos, hs * (d - np.log1p(d)), L)))
    best = None
    starts = [x0] if x0 is not None else [(max(h.sum(), 1.0), a0, 0.0, 0.0) for a0 in (0.5, 1.5)]
    for st in starts:
        m = Minuit(f, B=max(st[0], 1.0), a=st[1], c1=st[2], c2=st[3])
        m.errordef = 0.5
        m.errors["B"] = max(np.sqrt(h.sum()), 1.0)
        m.limits["B"] = (0, None); m.limits["a"] = (0.01, 3); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40)
        m.fixed["c2"] = nco < 2; m.fixed["c1"] = nco < 1
        m.migrad(ncall=4000)
        if best is None or m.fval < best.fval:
            best = m
    return [float(best.values[k]) for k in ("B", "a", "c1", "c2")]


def count_yields(H, mu, W=3.0, side=3.0, nco=None, x0=None, emax=None):
    """H: (ncat, NB). S_c = N(|dm-mu|<W) − tło(fit pasm |dm-mu|>side [, dm<emax]) w oknie. Zwraca S, lista parametrów tła, nco."""
    win = np.abs(R.CEN - mu) < W
    sb = np.abs(R.CEN - mu) > side
    if emax is not None:
        sb &= R.CEN < emax
    if nco is None:
        nsb = H[:, sb].sum(1)
        nco = np.where(nsb >= R.NSB_FULL, 2, np.where(nsb >= R.NSB_LIN, 1, 0))
    S, P = [], []
    for c in range(H.shape[0]):
        p = fit_bkg(H[c], sb, int(nco[c]), None if x0 is None else x0[c])
        S.append(H[c][win].sum() - p[0] * R.bkg_shape(*p[1:])[win].sum()); P.append(p)
    return np.array(S), P, np.asarray(nco)


def dm_bin(dm):
    ib = np.floor((dm - R.E[0]) / 0.1 + 1e-9).astype(int)
    return ib, (ib >= 0) & (ib < R.NB)


def derived(S, fm):
    """Z wydajków danych S (ncat, kategorie rosnąco w score, ostatnia = > WP25 dla 6 kat.) i ułamków MC fm: R_c, R_top, rho."""
    f = S / S.sum()
    out = {"R": f / fm, "N": S.sum()}
    if len(S) == 6:
        out["rho50"] = (S[5] / (S[4] + S[5])) / (fm[5] / (fm[4] + fm[5]))          # N(>WP25)/N(>WP50), dane/MC
        out["rho70"] = (S[5] / S[3:].sum()) / (fm[5] / fm[3:].sum())                # N(>WP25)/N(>WP70)
        out["SF50"] = (S[4:].sum() / S.sum()) / fm[4:].sum()
        out["SF70"] = (S[3:].sum() / S.sum()) / fm[3:].sum()
    return out


PT_EDGES = {428770: [8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 23, 30], 427929: [8, 10, 12, 15, 30]}


def pt_weights(run, W=3.0):
    """Wagi MC w pT(D0), by widmo sygnału MC odpowiadało danym. Widmo danych: wydajek zliczeniem dla score >= WP70 w binach pT,
    podzielony przez wydajność MC cięcia WP70 w tym binie (założenie: SF(WP70) niezależne od pT; R3: SF70 = 0,99 ± 0,04).
    Zwraca (wagi dla wszystkich wierszy MC [0 poza 8-30], tabela)."""
    d = load(f"data_{run}"); M = load("fd_nd")
    mu = MU["8-30"]; t70 = THR["8-30"]["0.70"]
    ed = np.array(PT_EDGES[run], float)
    sd = base_mask(d, "8-30", isdata=True, run=run) & (d["score"] >= t70)
    sm = base_mask(M, "8-30", isdata=False)
    ib, ok = dm_bin(d["dm"][sd])
    pb = np.clip(np.searchsorted(ed, d["pt"][sd], side="right") - 1, 0, len(ed) - 2)
    H = np.bincount((pb * R.NB + ib)[ok], minlength=(len(ed) - 1) * R.NB).reshape(len(ed) - 1, R.NB).astype(float)
    S, _, _ = count_yields(H, mu, W, 3.0)
    pm = np.clip(np.searchsorted(ed, M["pt"], side="right") - 1, 0, len(ed) - 2)
    inw = sm & (np.abs(M["dm"] - mu) < W)
    n_all = np.bincount(pm[inw], minlength=len(ed) - 1).astype(float)
    n_p = np.bincount(pm[inw & (M["score"] >= t70)], minlength=len(ed) - 1).astype(float)
    Sd = np.clip(S, 0, None) / (n_p / n_all)
    wbin = (Sd / Sd.sum()) / (n_all / n_all.sum())
    w = np.where(sm, wbin[pm], 0.0)
    tab = {"edges": ed.tolist(), "S70": S.tolist(), "eff70_mc": (n_p / n_all).tolist(), "frac_data": (Sd / Sd.sum()).tolist(),
           "frac_mc": (n_all / n_all.sum()).tolist(), "w": wbin.tolist()}
    return w, tab
