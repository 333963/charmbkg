"""DLZ_d0: nowe cechy D0→Kπ dopisane do ramek FEAT_charmnet (h5 v1.2; FEAT bez edycji) → DLZ_d0/frames/<próbka>_D0.pkl.
Córka 0 = K, 1 = π (hipoteza kandydata). Jednostki: pT GeV, d0/z0/Lxy mm.

kinematyczne (mogą rekonstruować masę → test rzeźbienia):
  ctsK_nf      cosθ* ze znakiem (H ma tylko |cosθ*|)
  ptAsym_nf    (pT_K − pT_π)/(pT_K + pT_π)
  dMinPt_nf, dMaxPt_nf   min/max pT córek [GeV];  dMinPtRatio_cn, dMaxPtRatio_cn (FEAT)
  dR_nf        ΔR(K, π)
parametr zderzenia / wierzchołek (nie zawierają pędu ani masy poza σ(pT)):
  d0prod_nf    d0_K · d0_π  (ze znakiem; sygnał: przeciwne znaki)
  d0sigK_nf, d0sigPi_nf   |d0|/σ_param, σ_param = 0,012 mm ⊕ 0,075 mm·GeV/(pT √sinθ)  — PRZYBLIŻONA parametryzacja (h5 nie ma błędów d0)
  d0sigMin_nf, d0sigMax_nf, d0chi2_nf = Σ (d0/σ)²
  dz0_nf       |z0sinθ_K − z0sinθ_π|  (zastępnik DCA torów wzdłuż z; DCA 3D nie ma w h5)
  absImpact_cn, d0AbsZ0_cn, d1AbsZ0_cn (FEAT)
czas:
  ctau_nf      L_xy · m0(PDG) / pT  [mm]  (m0 stałe — bez masy kandydata)
PID (dE/dx pikselowe, surowe — bez przesunięcia MC, spójnie z cechami H):
  bbK_K_nf, bbK_pi_nf    log10(dE/dx_K-tor / BB(p, m_K)), … / BB(p, m_π)
  bbPi_pi_nf, bbPi_K_nf  to samo dla toru π;  NaN gdy npts ≤ 2 lub dE/dx poza (0,1; 8)
  pPi_nf = log p toru π, nptsPi_nf, dedxPi_nf
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

FEAT = os.environ.get("FEAT_DIR", "/home/k/Documents/charm/charm_work/FEAT_charmnet")
FRIN = os.environ.get("D0_FRAMES_IN", os.path.join(FEAT, "frames"))
OUT = os.path.join(os.environ.get("D0_BASE") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "frames")
_bb = importlib.util.spec_from_file_location("bb", os.path.join(os.environ.get("FEAT_SCRIPTS", os.path.join(FEAT, "scripts")), "ext", "bethe_bloch_charmpy.py"))
BB = importlib.util.module_from_spec(_bb)
_bb.loader.exec_module(BB)
M0 = 1.86484
A_D0, B_D0 = 0.012, 0.075


def add(L, is_data):
    pt = [L[f"d_pt{i}"].to_numpy(float) for i in (0, 1)]
    eta = [np.clip(L[f"d_eta{i}"].to_numpy(float), -5, 5) for i in (0, 1)]
    phi = [L[f"d_phi{i}"].to_numpy(float) for i in (0, 1)]
    d0 = [L[f"d_d0{i}"].to_numpy(float) for i in (0, 1)]
    z0 = [L[f"d_z0{i}"].to_numpy(float) for i in (0, 1)]
    L["ctsK_nf"] = L["cos_theta_star"].to_numpy(float)
    L["ptAsym_nf"] = (pt[0] - pt[1]) / (pt[0] + pt[1])
    L["dMinPt_nf"] = np.minimum(pt[0], pt[1])
    L["dMaxPt_nf"] = np.maximum(pt[0], pt[1])
    dphi = np.abs(phi[0] - phi[1])
    dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
    L["dR_nf"] = np.hypot(eta[0] - eta[1], dphi)
    L["d0prod_nf"] = d0[0] * d0[1]
    sig = [np.hypot(A_D0, B_D0 / (pt[i] * np.sqrt(1 / np.cosh(eta[i])))) for i in (0, 1)]
    s = [np.abs(d0[i]) / sig[i] for i in (0, 1)]
    L["d0sigK_nf"], L["d0sigPi_nf"] = s[0], s[1]
    L["d0sigMin_nf"], L["d0sigMax_nf"] = np.minimum(s[0], s[1]), np.maximum(s[0], s[1])
    L["d0chi2_nf"] = s[0] ** 2 + s[1] ** 2
    L["dz0_nf"] = np.abs(z0[0] - z0[1])
    L["ctau_nf"] = L["Lxy"].to_numpy(float) * M0 / L["pt"].to_numpy(float)
    sh = 0.0  # BEZ przesunięcia MC: cechy H (dedx_K, dedx_pi_max) są surowe; przesunięta kopia obok surowej = deterministyczny znacznik MC
    #           (test 17.11: C2ST tła 0,570 → 0,999 z przesunięciem −0,023). Kolumna kontrolna dedxPiSh_nf zachowuje przesunięcie (tylko do modelling).
    L["dedxPiSh_nf"] = (L["d_dedx1"].to_numpy(float) + (0.0 if is_data else -0.023)).astype(np.float32)
    p = [pt[i] * np.cosh(eta[i]) for i in (0, 1)]
    for i, nm in ((0, "K"), (1, "Pi")):
        dedx = L[f"d_dedx{i}"].to_numpy(float) + sh
        npts = L[f"d_npts{i}"].to_numpy(float)
        v = (npts > 2) & (dedx > 0.1) & (dedx < 8.0) & (p[i] > 0.1) & (p[i] < 50)
        for hyp, mass in (("K", BB.M_K), ("pi", BB.M_PI)):
            x = np.full(len(L), np.nan)
            x[v] = BB.distance(mass, p[i][v], dedx[v])
            L[f"bb{nm}_{hyp}_nf"] = x
    L["pPi_nf"] = np.log(p[1])
    L["nptsPi_nf"] = L["d_npts1"].to_numpy(float)
    L["dedxPi_nf"] = L["d_dedx1"].to_numpy(float) + sh
    for c in L.columns:
        if c.endswith("_nf"):
            L[c] = L[c].astype(np.float32)
    return L


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for s in ("fd_nd", "mb_nd", "data_428770"):
        L = pd.read_pickle(os.path.join(FRIN, f"{s}_D0.pkl"))
        L = add(L, s.startswith("data"))
        L.to_pickle(os.path.join(OUT, f"{s}_D0.pkl"))
        print(s, len(L), flush=True)
