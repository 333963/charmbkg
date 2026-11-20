"""PEAK_lxy010: stosunek odbicie/sygnał D⁰ K↔π z tagu D* (metoda SWAP_d0/swap_dstar.py: okno Δm − sidebandy, hipoteza z ładunku π_s)
przy L_xy ≥ 0,25, ≥ 0,10 i w plasterku 0,10–0,25 mm. hypo_from_tag = 0: zgodne z tagiem, = 1: odbicie.
R = N_odbicie(1750–2000) / N_zgodne(okno 1840–1890). Kształt odbicia: χ² znormalizowanego histogramu względem MC fd (25 binów po 10 MeV)
oraz plaster 0,10–0,25 vs ≥ 0,25 w tych samych danych (próbki rozłączne).
Dane: 427929 = out/dstar_data_427929_v20.parquet (h5 v2.0, 226/458 plików, L_xy ≥ 0,10), 428770 = HYP_d0/runs/lxplus_dstar (v1.2, L_xy ≥ 0,25).
MC: fd_nd, mb_nd z HYP_d0/runs/lxplus_dstar (TYLKO L_xy ≥ 0,25).  Wynik: out/refl_rows.csv, out/refl_shape.csv"""
import os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
W = "/home/k/Documents/charm/charm_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
MPI = 0.13957
WIN, SB = (0.1435, 0.1475), (0.1500, 0.1650)
E = np.arange(1750, 2001, 10.0); C = 0.5 * (E[1:] + E[:-1])
REG = {"low": (1750, 1840), "win": (1840, 1890), "high": (1890, 2000)}
cols = ["m", "score", "dm", "hypo_from_tag", "pt", "Lxy"]
SRC = {"data_427929": f"{OUT}/dstar_data_427929_v20.parquet", "data_428770": f"{W}/HYP_d0/runs/lxplus_dstar/data_428770.parquet",
       "fd_nd": f"{W}/HYP_d0/runs/lxplus_dstar/fd_nd.parquet", "mb_nd": f"{W}/HYP_d0/runs/lxplus_dstar/mb_nd.parquet"}
CUTS = {"ge025": (0.25, 1e9), "ge010": (0.10, 1e9), "slice": (0.10, 0.25)}
PTB = [(3, 30), (3, 5), (5, 8), (8, 30)]


def measure(d):
    """jak swap_dstar.py: tło Δm a·(Δm−mπ)^b z sidebandów → k = I(okno)/I(SB); histogramy m w oknie − k·SB"""
    he, ee = np.histogram(d.dm, np.arange(0.1396, 0.1651, 0.0005)); ce = 0.5 * (ee[1:] + ee[:-1])
    o = (ce < 0.1420) | (ce > 0.1500)
    f = lambda x, a, b: a * (x - MPI) ** b  # noqa
    p, _ = curve_fit(f, ce[o], he[o], p0=[max(he[o].max(), 1) * 10, 0.5], sigma=np.sqrt(np.maximum(he[o], 1)), maxfev=20000)
    xs = np.linspace(MPI + 1e-6, 0.165, 20001)
    I = lambda a, b: np.trapezoid(f(xs[(xs > a) & (xs < b)], *p), xs[(xs > a) & (xs < b)])  # noqa
    kf = I(*WIN) / I(*SB)
    inw, insb = (d.dm > WIN[0]) & (d.dm < WIN[1]), (d.dm > SB[0]) & (d.dm < SB[1])
    r = {"k_sb": kf}
    for h in (0, 1):
        g = d.hypo_from_tag == h
        hw = np.histogram(d.m[inw & g] * 1000, E)[0].astype(float); hs = np.histogram(d.m[insb & g] * 1000, E)[0].astype(float)
        r[f"h{h}"] = hw - kf * hs; r[f"e{h}"] = np.sqrt(hw + kf ** 2 * hs)
    w = (C > REG["win"][0]) & (C < REG["win"][1]); lo = C < REG["win"][0]
    n0, e0 = r["h0"][w].sum(), np.sqrt((r["e0"][w] ** 2).sum())
    n1, e1 = r["h1"].sum(), np.sqrt((r["e1"] ** 2).sum())
    nl, el = r["h1"][lo].sum(), np.sqrt((r["e1"][lo] ** 2).sum())
    r.update(n0=n0, e_n0=e0, n1=n1, e_n1=e1, R=n1 / n0, e_R=abs(n1 / n0) * np.hypot(e1 / n1, e0 / n0) if n1 and n0 else np.nan,
             Rlow=nl / n0, e_Rlow=abs(nl / n0) * np.hypot(el / nl, e0 / n0) if nl and n0 else np.nan)
    return r


def chi2_shape(a, b):
    """χ² kształtu odbicia (unormowane do pola 1) a vs b; ndf = 24. Sensowne tylko przy N_odbicie ≳ 200 (kolumna n_refl)."""
    fl = lambda r: np.sqrt(r["e1"] ** 2 + 1.0)   # +1 zdarzenie: podłoga błędu dla pustych binów  # noqa
    ha, ea, hb, eb = a["h1"] / a["h1"].sum(), fl(a) / abs(a["h1"].sum()), b["h1"] / b["h1"].sum(), fl(b) / abs(b["h1"].sum())
    return float(np.sum((ha - hb) ** 2 / (ea ** 2 + eb ** 2))), len(C) - 1


L = {k: pd.read_parquet(v, columns=cols) for k, v in SRC.items()}
rows, shp, M = [], [], {}
for sc in (0.857, 0.387):
    for lo, hi in PTB:
        for k, d0 in L.items():
            d = d0[(d0.score > sc) & (d0.pt >= lo) & (d0.pt < hi)]
            for cut, (a, b) in CUTS.items():
                if d0.Lxy.min() > 0.2 and cut != "ge025":
                    continue
                x = d[(d.Lxy >= a) & (d.Lxy < b)]
                try:
                    r = measure(x)
                except Exception as ex:  # noqa
                    print("FAIL", sc, lo, hi, k, cut, ex); continue
                M[sc, lo, hi, k, cut] = r
                rows.append(dict(score=sc, bin=f"pT{lo}-{hi}", sample=k, cut=cut, n_tag_win=r["n0"], e_n_tag_win=r["e_n0"], n_refl=r["n1"],
                                 e_n_refl=r["e_n1"], R=r["R"], e_R=r["e_R"], Rlow=r["Rlow"], e_Rlow=r["e_Rlow"], k_sb=r["k_sb"]))
T = pd.DataFrame(rows)
# stosunki dane/MC (MC zawsze ≥ 0,25)
for mc in ("fd_nd", "mb_nd"):
    ref = T[(T["sample"] == mc)].set_index(["score", "bin"])[["R", "e_R"]]
    j = T.join(ref, on=["score", "bin"], rsuffix=f"_{mc}")
    T[f"R_over_{mc}"] = j["R"] / j[f"R_{mc}"]; T[f"e_R_over_{mc}"] = T[f"R_over_{mc}"] * np.hypot(j["e_R"] / j["R"], j[f"e_R_{mc}"] / j[f"R_{mc}"])
T.to_csv(f"{OUT}/refl_rows.csv", index=False)
for (sc, lo, hi, k, cut), r in M.items():
    if not k.startswith("data"):
        continue
    for ref, key in (("fd_nd@ge025", (sc, lo, hi, "fd_nd", "ge025")), ("mb_nd@ge025", (sc, lo, hi, "mb_nd", "ge025")),
                     ("same_data@ge025", (sc, lo, hi, k, "ge025"))):
        if ref.startswith("same") and cut != "slice":
            continue
        c2, nd = chi2_shape(r, M[key])
        shp.append(dict(score=sc, bin=f"pT{lo}-{hi}", sample=k, cut=cut, vs=ref, chi2=c2, ndf=nd, n_refl=r["n1"], n_refl_ref=M[key]["n1"]))
Sh = pd.DataFrame(shp); Sh.to_csv(f"{OUT}/refl_shape.csv", index=False)
pd.set_option("display.width", 250)
print(T.round(3).to_string()); print(Sh.round(1).to_string())
