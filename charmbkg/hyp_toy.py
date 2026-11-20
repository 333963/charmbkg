"""HYP: toy MC D0→Kπ z obiema hipotezami K/π na parę torów (schemat jak parquet z `extract_d0hyp.py`).

Sygnał: D0 (m0, rozpad izotropowy, p* = 0,861 GeV), tło: „rodzice" o płaskim rozkładzie masy z losowych par (π/K/p).
Pomiar: pT torów z rozdzielczością 2 % (σ_m ≈ 15–25 MeV), dE/dx = skala·BB(p, m_prawdziwa)·lognormal(0, 12 %), npts 0–5.
Obie hipotezy (tor a = K / tor b = K) → osobne wiersze, gdy m w [1,75; 2,00] GeV. Score „D0_H" niezależny od masy,
ale zależny od pulla K hipotezy wiersza (jak cechy dedx_K w H)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit

from .constants import M_D0, M_K, M_PI
from .bethe_bloch import dedx_expected
from .hyp import DEDX_SCALE

M_P = 0.93827
LO, HI = 1.75, 2.00


def _decay(rng, M, m1, m2, pt, eta, phi):
    n = len(M)
    ps = np.sqrt(np.clip((M ** 2 - (m1 + m2) ** 2) * (M ** 2 - (m1 - m2) ** 2), 0, None)) / (2 * M)
    ct = rng.uniform(-1, 1, n)
    st = np.sqrt(1 - ct ** 2)
    fz = rng.uniform(0, 2 * np.pi, n)
    p1 = np.stack([ps * st * np.cos(fz), ps * st * np.sin(fz), ps * ct], 1)
    P = np.stack([pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)], 1)
    E = np.sqrt((P ** 2).sum(1) + M ** 2)
    b = P / E[:, None]
    b2 = (b ** 2).sum(1)
    g = 1 / np.sqrt(1 - b2)

    def boost(p, m):
        e = np.sqrt((p ** 2).sum(1) + m ** 2)
        bp = (b * p).sum(1)
        k = np.where(b2 > 0, (g - 1) * bp / b2, 0)
        return p + (k + g * e)[:, None] * b

    return boost(p1, m1), boost(-p1, m2), ct


def _kin(p):
    pt = np.hypot(p[:, 0], p[:, 1])
    return pt, np.arcsinh(p[:, 2] / pt), np.arctan2(p[:, 1], p[:, 0]), np.linalg.norm(p, axis=1)


def _mass(pa, pb, ma, mb):
    ea = np.sqrt((pa ** 2).sum(1) + ma ** 2)
    eb = np.sqrt((pb ** 2).sum(1) + mb ** 2)
    s = pa + pb
    return np.sqrt(np.clip((ea + eb) ** 2 - (s ** 2).sum(1), 0, None))


def _cts(pk, ppi, mk, mpi):
    """cosθ* kaonu (hipoteza) w układzie spoczynkowym pary."""
    ek = np.sqrt((pk ** 2).sum(1) + mk ** 2)
    epi = np.sqrt((ppi ** 2).sum(1) + mpi ** 2)
    P = pk + ppi
    E = ek + epi
    M = np.sqrt(np.clip(E ** 2 - (P ** 2).sum(1), 1e-9, None))
    pn = np.linalg.norm(P, axis=1)
    u = P / pn[:, None]
    g, bg = E / M, pn / M
    pl = (pk * u).sum(1)
    pl_star = g * pl - bg * ek
    pstar = np.sqrt(np.clip((M ** 2 - (mk + mpi) ** 2) * (M ** 2 - (mk - mpi) ** 2), 0, None)) / (2 * M)
    return np.clip(pl_star / np.maximum(pstar, 1e-6), -1, 1)


def _pdg(mt):
    return np.where(mt == M_K, 321, np.where(mt == M_PI, 211, 2212)).astype(np.int64)


def _dedx(rng, p, mtrue):
    n = len(p)
    npts = np.clip(rng.poisson(2.6, n) + 1, 1, 6)
    npts[rng.random(n) < 0.05] = 0
    res = 0.12 * np.sqrt(3.0 / np.maximum(npts, 1))
    dx = DEDX_SCALE * dedx_expected(p, mtrue) * np.exp(rng.normal(0, res))
    dx[npts == 0] = -1.0
    return dx, npts


def make(sample: str, n_sig: int, n_bkg: int, seed: int = 1, row0: int = 0, frac_b: float = 0.2, res: float = 0.020,
         pt_lo: float = 3.0, d_pt_min: float = 0.5, eta_max: float = 2.4) -> pd.DataFrame:
    """Kandydaci D0 (wiersze = hipotezy) próbki `sample` (fd_nd | mb_nd | data_428770). res = rozdzielczość pT torów,
    pt_lo = min pT D, d_pt_min = cięcie pT córek (testy: res duże / pt_lo małe → reżim, w którym dE/dx wygrywa z masą)."""
    rng = np.random.default_rng(seed)
    parts = []
    for kind, N in (("sig", n_sig), ("bkg", n_bkg)):
        if N <= 0:
            continue
        pt = np.clip(pt_lo + rng.exponential(3.2 if kind == "sig" else 2.4, N), pt_lo, 29.9)
        eta = rng.uniform(-eta_max, eta_max, N)
        phi = rng.uniform(-np.pi, np.pi, N)
        if kind == "sig":
            M = np.full(N, M_D0)
            m1, m2 = np.full(N, M_K), np.full(N, M_PI)          # tor 1 = K, tor 2 = π
        else:
            M = rng.uniform(1.30, 2.60, N)
            u = rng.random(N)
            m1 = np.where(u < 0.6, M_PI, np.where(u < 0.85, M_K, M_P))
            m2 = np.where(rng.random(N) < 0.85, M_PI, M_K)
        p1, p2, _ = _decay(rng, M, m1, m2, pt, eta, phi)
        for p in (p1, p2):
            p *= (1 + rng.normal(0, res, N))[:, None]
        ok = (np.hypot(p1[:, 0], p1[:, 1]) > d_pt_min) & (np.hypot(p2[:, 0], p2[:, 1]) > d_pt_min)
        dx1, n1 = _dedx(rng, np.linalg.norm(p1, axis=1), m1)
        dx2, n2 = _dedx(rng, np.linalg.norm(p2, axis=1), m2)
        cid = np.arange(N)
        for hyp in (0, 1):                                     # 0: tor1 = K (poprawna dla sygnału), 1: tor2 = K
            pa, pb = (p1, p2) if hyp == 0 else (p2, p1)
            da, db = (dx1, dx2) if hyp == 0 else (dx2, dx1)
            na, nb = (n1, n2) if hyp == 0 else (n2, n1)
            ta, tb = (m1, m2) if hyp == 0 else (m2, m1)                # prawdziwa masa toru → d_tpdg (v2 pertrack)
            m = _mass(pa, pb, M_K, M_PI)
            sel = ok & (m >= LO) & (m < HI)
            a = _kin(pa[sel])
            b = _kin(pb[sel])
            P = pa[sel] + pb[sel]
            Dpt, Deta, Dphi, _ = _kin(P)
            parts.append(pd.DataFrame({
                "kind": kind, "cid": cid[sel], "hyp": hyp, "m": m[sel], "pt": Dpt, "eta": Deta, "phi": Dphi,
                "d_pt0": a[0], "d_eta0": a[1], "d_phi0": a[2], "d_pt1": b[0], "d_eta1": b[1], "d_phi1": b[2],
                "d_dedx0": da[sel], "d_dedx1": db[sel], "d_npts0": na[sel], "d_npts1": nb[sel],
                "cos_theta_star": _cts(pa[sel], pb[sel], M_K, M_PI),
                "d_tpdg0": _pdg(ta[sel]), "d_tpdg1": _pdg(tb[sel])}))
    L = pd.concat(parts, ignore_index=True)
    sig = (L["kind"] == "sig").to_numpy()
    uid = np.where(sig, 0, 1) * 10_000_000 + L["cid"].to_numpy()
    L["run"] = 900000 + {"fd_nd": 1, "mb_nd": 2}.get(sample, 3)
    L["event"] = (uid // 3).astype(np.int64) + 50_000_000 * (~sig)
    L["d_trk0"] = np.where(L["hyp"] == 0, 2 * uid, 2 * uid + 1).astype(np.int64)
    L["d_trk1"] = np.where(L["hyp"] == 0, 2 * uid + 1, 2 * uid).astype(np.int64)
    L = L.sort_values(["event", "cid", "kind", "hyp"]).reset_index(drop=True)
    sig = (L["kind"] == "sig").to_numpy()
    from .hyp import _mix64
    isb = (_mix64((np.where(sig, 0, 1) * 10_000_000 + L["cid"].to_numpy()).astype(np.uint64) + np.uint64(seed)) % np.uint64(1000)).astype(int) < 1000 * frac_b
    L["truth"] = np.where(sig, np.where(isb, 1.0, 0.0), 8.0)   # c/b per kandydat (obie hipotezy tak samo)
    L["hypo_swap"] = np.where(sig, L["hyp"].to_numpy(), -1).astype(np.int32)
    from .bethe_bloch import pull
    zk = np.nan_to_num(pull(np.where(L["d_dedx0"] > 0, L["d_dedx0"], np.nan), L["d_pt0"] * np.cosh(L["d_eta0"]), M_K, DEDX_SCALE))
    lat = 2.6 * sig - 1.3 + rng.normal(0, 1.0, len(L)) - 0.10 * np.minimum(zk ** 2, 16)
    L["score"] = expit(lat).astype(np.float32)
    L["fold"] = -1
    L["sample"] = sample
    L["row"] = row0 + np.arange(len(L), dtype=np.int64)
    if sample.startswith("data"):
        L["truth_hidden"] = L["truth"]
        L["hypo_swap_hidden"] = L["hypo_swap"]
        L["truth"] = np.nan
        L["hypo_swap"] = -1
    return L.drop(columns=["kind", "cid", "hyp"])


def make_all(scale: float = 1.0, seed: int = 7, **kw) -> dict:
    s = lambda x: int(x * scale)  # noqa: E731
    return {"fd_nd": make("fd_nd", s(60000), s(60000), seed, **kw),
            "mb_nd": make("mb_nd", s(1500), s(150000), seed + 1, **kw),
            "data_428770": make("data_428770", s(14000), s(450000), seed + 2, **kw)}


def write_h5(L: pd.DataFrame, path: str, n_other: int = 500, seed: int = 0):
    """Toy → h5 w układzie v1.2 (compound `events`, nazwy ntupli, pT córek w MeV) + `n_other` wierszy D⁺ (decayType 4)
    wplecionych między zdarzenia (test filtra i globalnego `row`). Zwraca ramkę z kolumną `h5row` (wiersz w pliku)."""
    import h5py
    rng = np.random.default_rng(seed)
    n = len(L)
    hidden = "truth_hidden" in L
    t = L["truth_hidden"] if hidden else L["truth"]
    hs = L["hypo_swap_hidden"] if hidden else L["hypo_swap"]
    d = {"run_number": L["run"].to_numpy(np.int32), "event_number": L["event"].to_numpy(np.int64),
         "event_weight": np.ones(n, np.float32), "Dmeson_decayType": np.full(n, 10, np.int32),
         "Dmeson_mass": L["m"].to_numpy(np.float64), "Dmeson_pt": L["pt"].to_numpy(np.float32), "Dmeson_eta": L["eta"].to_numpy(np.float32),
         "Dmeson_phi": L["phi"].to_numpy(np.float32), "Dmeson_charge": np.zeros(n, np.int32), "Dmeson_delta_mass": np.full(n, -0.999),
         "Dmeson_Lxy": rng.exponential(1.0, n) + 0.25, "Dmeson_chi2": rng.chisquare(3, n), "Dmeson_impact": rng.normal(0, 0.03, n),
         "Dmeson_impactSignificance": rng.normal(0, 1.5, n), "Dmeson_z0SinTheta": rng.normal(0, 0.05, n),
         "Dmeson_cosThetaStar": L["cos_theta_star"].to_numpy(np.float32), "Dmeson_cosThetaP": 1 - rng.exponential(0.002, n),
         "Dmeson_ptcone40": rng.exponential(1.0, n), "Dmeson_truth_category": np.nan_to_num(t.to_numpy(float) if not hidden else
                                                                                           np.full(n, -1.0), nan=-1).astype(np.int32),
         "Dmeson_mKpi1": np.full(n, -0.999, np.float32), "Dmeson_mKpi2": np.full(n, -0.999, np.float32),
         "Dmeson_hypo_swap": (hs.to_numpy() if not hidden else np.full(n, -1)).astype(np.int32)}
    for i in range(3):
        on = i < 2
        g = lambda c, dflt=np.nan: L[f"{c}{i}"].to_numpy(float) if on else np.full(n, dflt)  # noqa: E731
        d[f"Dmeson_daughter{i}_pt"] = g("d_pt", 0.0) * 1000.0
        d[f"Dmeson_daughter{i}_eta"] = g("d_eta")
        d[f"Dmeson_daughter{i}_phi"] = g("d_phi")
        tid = L[f"d_trk{i}"].to_numpy(np.int64) if on else np.full(n, -999)
        d[f"Dmeson_daughter{i}_d0"] = ((tid % 9973) * 1e-4).astype(np.float32) if on else np.zeros(n, np.float32)
        d[f"Dmeson_daughter{i}_z0SinThetaPV"] = ((tid // 9973 % 9973) * 1e-4).astype(np.float32) if on else np.zeros(n, np.float32)
        d[f"Dmeson_daughter{i}_npts_dEdx"] = (g("d_npts", 0) if on else np.zeros(n)).astype(np.int32)
        d[f"Dmeson_daughter{i}_dEdx"] = g("d_dedx", -1.0) if on else np.full(n, -1.0)
        d[f"Dmeson_daughter{i}_trackId"] = (tid % 2_000_000_000).astype(np.int32)
        kp = np.where(L["hypo_swap"].to_numpy() if not hidden else hs.to_numpy(), 1, 0)
        d[f"Dmeson_daughter{i}_truthPdgId"] = (np.where((i == 0) ^ (kp == 1), 321, 211) if on else np.full(n, -999)).astype(np.int32)
    # wplecione D⁺
    k = rng.integers(0, n + 1, n_other)
    ins = np.sort(k)
    order = np.argsort(np.concatenate([np.arange(n) + np.searchsorted(ins, np.arange(n), side="right"),
                                       ins + np.arange(n_other)]), kind="stable")
    dt = np.dtype([(c, v.dtype) for c, v in d.items()])
    arr = np.zeros(n + n_other, dt)
    tot = np.concatenate([np.arange(n), np.full(n_other, -1)])[order]
    for c, v in d.items():
        full = np.concatenate([v, np.zeros(n_other, v.dtype)])[order]
        arr[c] = full
    other = tot < 0
    arr["Dmeson_decayType"][other] = 4
    arr["Dmeson_mass"][other] = 1.87
    with h5py.File(path, "w") as f:
        f.create_dataset("events", data=arr)
    out = L.copy()
    out["h5row"] = np.flatnonzero(~other)[np.argsort(tot[~other])]
    return out
