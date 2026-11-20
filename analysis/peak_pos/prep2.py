"""PEAK_pos: wyciągi z większych próbek (30 % lumi danych, model analizy D0_H / Dplus_H):
 D0: HYP_d0/runs/lxplus/{data_428770,fd_nd}.parquet  (q = zapach hipotezy = −ładunek toru K = −d_charge0)
 D+: SWAP_d0/runs/lxplus/{data_428770_dp,mb_nd_dp}.parquet (score ≥ 0,59 wstępnie; weto D*: dstar_veto==0)
Dokłada m_swap (D0) i m_rec (kontrola) z kinematyki córek."""
import os
import numpy as np, pandas as pd, pyarrow.parquet as pq
W = "/home/k/Documents/charm/charm_work"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
MK, MPI = 0.493677, 0.13957039


def p4(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return np.sqrt(px**2 + py**2 + pz**2 + m**2), px, py, pz


def minv(*v):
    E = sum(a[0] for a in v); px = sum(a[1] for a in v); py = sum(a[2] for a in v); pz = sum(a[3] for a in v)
    return np.sqrt(np.maximum(E**2 - px**2 - py**2 - pz**2, 0))


def kin(L, nd):
    return [(L[f"d_pt{i}"].to_numpy(float), L[f"d_eta{i}"].to_numpy(float), L[f"d_phi{i}"].to_numpy(float)) for i in range(nd)]


D0C = ["m", "pt", "eta", "phi", "score", "weight", "truth", "hypo_swap", "d_charge0", "d_charge1", "evblk", "file"] + \
      [f"d_{k}{i}" for k in ("pt", "eta", "phi") for i in range(2)]
for s, tag in [("data_428770", "data"), ("fd_nd", "mc")]:
    L = pq.read_table(f"{W}/HYP_d0/runs/lxplus/{s}.parquet", columns=D0C).to_pandas()
    k = kin(L, 2)
    L["m_rec"] = minv(p4(*k[0], MK), p4(*k[1], MPI)).astype(np.float32)
    L["m_swap"] = minv(p4(*k[0], MPI), p4(*k[1], MK)).astype(np.float32)
    L["q"] = (-L["d_charge0"]).astype(np.int8)
    if tag == "mc":
        L = L[np.isin(L["truth"], [0, 1])].reset_index(drop=True)   # sygnał (oba przypisania; hypo_swap rozróżnia)
    L.to_parquet(f"{OUT}/D0_{tag}30.parquet")
    print("D0", s, len(L), "max|m_rec-m|", float(np.nanmax(np.abs(L.m_rec - L.m))), "q", L.q.value_counts().to_dict())
    del L
DPC = ["m", "pt", "eta", "phi", "score", "weight", "truth", "charge", "d_charge0", "dstar_veto", "evblk", "file"] + \
      [f"d_{k}{i}" for k in ("pt", "eta", "phi") for i in range(3)]
for s, tag in [("data_428770_dp", "data"), ("mb_nd_dp", "mbmc")]:
    L = pq.read_table(f"{W}/SWAP_d0/runs/lxplus/{s}.parquet", columns=DPC).to_pandas()
    L = L[L["dstar_veto"] == 0].reset_index(drop=True)
    if tag == "mbmc":
        L = L[np.isin(L["truth"], [0, 1])].reset_index(drop=True)
    k = kin(L, 3)
    L["m_rec"] = minv(p4(*k[0], MK), p4(*k[1], MPI), p4(*k[2], MPI)).astype(np.float32)
    L["q"] = L["charge"].astype(np.int8)
    L.to_parquet(f"{OUT}/Dplus_{tag}30.parquet")
    print("Dplus", s, len(L), "max|m_rec-m|", float(np.nanmax(np.abs(L.m_rec - L.m))))
    del L
