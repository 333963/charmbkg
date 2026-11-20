"""PEAK_lxy010 v2 — odbicie K↔π z tagu D* (metoda ../refl_lxy.py = SWAP_d0/swap_dstar.py: okno Δm − sidebandy, hipoteza z ładunku π_s)
przy L_xy ≥ 0,25 / ≥ 0,10 / plaster 0,10–0,25 mm, dane OBU runów i MC (fd_nd, mb_nd) przy TYM SAMYM cięciu L_xy.
R = N_odbicie(1750–2000 MeV) / N_zgodne(1840–1890 MeV); Rlow = część niskomasowa odbicia (1750–1840) / N_zgodne.
Selekcje klasyfikatora: wp70 (progi D0_H per bin pT), nowp (bez), s857 i s387 (score > 0,857 / 0,387 — jak wersja wstępna i SWAP_d0).
Błędy: analityczne (Poisson po odjęciu sidebandów) ORAZ bootstrap po zdarzeniach (Poisson(1) na evblk; wiersz zgodny i odbity tej samej pary
torów są w tym samym zdarzeniu, więc błąd analityczny R jest zawyżony). Kształt: χ² unormowanego histogramu odbić (25 binów, ndf 24).
Użycie: python refl_v2.py [nboot]   → out/refl_rows.csv, out/refl_shape.csv, out/refl_boot.csv"""
import os, sys, time
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
W = "/home/k/Documents/charm/charm_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"; IN = f"{HERE}/../runs/lxplus"
sys.path.insert(0, f"{W}/PEAK_pos")
NBOOT = int(sys.argv[1]) if len(sys.argv) > 1 else 0
sys.argv = sys.argv[:1]
from common import THR, pass_wp   # noqa: E402
MPI = 0.13957
WIN, SB = (0.1435, 0.1475), (0.1500, 0.1650)
EDM = np.arange(0.1396, 0.1651, 0.0005); CDM = 0.5 * (EDM[1:] + EDM[:-1]); NDM = len(CDM)
E = np.arange(1750, 2001, 10.0); C = 0.5 * (E[1:] + E[:-1]); NM = len(C)
REG = {"low": (1750, 1840), "win": (1840, 1890), "high": (1890, 2000)}
CUTS = {"ge025": (0.25, 1e9), "ge010": (0.10, 1e9), "slice": (0.10, 0.25)}
PTB = [(3, 30), (3, 5), (5, 8), (8, 30)]
SELS = ("wp70", "nowp", "s857", "s387")
SAMPLES = ("data_428770", "data_427929", "fd_nd", "mb_nd")
XS = np.linspace(MPI + 1e-6, 0.165, 20001)
OSB = (CDM < 0.1420) | (CDM > 0.1500)
f_bkg = lambda x, a, b: a * (x - MPI) ** b  # noqa: E731


def prep(d):
    """indeksy binów i maski wiersza (liczone raz na próbkę × selekcję × bin pT)"""
    dm, m = d.dm.to_numpy(), d.m.to_numpy() * 1000
    idm = np.clip(np.searchsorted(EDM, dm, side="right") - 1, 0, NDM - 1); okdm = (dm >= EDM[0]) & (dm < EDM[-1])
    im = np.clip(np.searchsorted(E, m, side="right") - 1, 0, NM - 1); okm = (m >= E[0]) & (m < E[-1])
    code, u = pd.factorize(d.evblk.to_numpy())
    return dict(idm=idm, okdm=okdm, im=im, okm=okm, inw=(dm > WIN[0]) & (dm < WIN[1]), insb=(dm > SB[0]) & (dm < SB[1]),
                h=d.hypo_from_tag.to_numpy(), L=d.Lxy.to_numpy(), code=code, nev=len(u))


def measure(P, mask, w=None):
    """jak swap_dstar.py: tło Δm a·(Δm−mπ)^b z sidebandów → k = I(okno)/I(SB); histogramy m w oknie − k·SB. w = wagi wierszy (bootstrap)."""
    ww = np.ones(len(mask)) if w is None else w
    s = mask & P["okdm"]
    he = np.bincount(P["idm"][s], weights=ww[s], minlength=NDM)
    p, _ = curve_fit(f_bkg, CDM[OSB], he[OSB], p0=[max(he[OSB].max(), 1) * 10, 0.5], sigma=np.sqrt(np.maximum(he[OSB], 1)), maxfev=20000)
    I = lambda a, b: np.trapezoid(f_bkg(XS[(XS > a) & (XS < b)], *p), XS[(XS > a) & (XS < b)])  # noqa: E731
    kf = I(*WIN) / I(*SB)
    r = {"k_sb": kf}
    for h in (0, 1):
        g = mask & P["okm"] & (P["h"] == h)
        hw = np.bincount(P["im"][g & P["inw"]], weights=ww[g & P["inw"]], minlength=NM)
        hs = np.bincount(P["im"][g & P["insb"]], weights=ww[g & P["insb"]], minlength=NM)
        r[f"h{h}"] = hw - kf * hs; r[f"e{h}"] = np.sqrt(hw + kf ** 2 * hs)
    wn = (C > REG["win"][0]) & (C < REG["win"][1]); lo = C < REG["win"][0]
    n0, e0 = r["h0"][wn].sum(), np.sqrt((r["e0"][wn] ** 2).sum())
    n1, e1 = r["h1"].sum(), np.sqrt((r["e1"] ** 2).sum())
    nl, el = r["h1"][lo].sum(), np.sqrt((r["e1"][lo] ** 2).sum())
    r.update(n0=n0, e_n0=e0, n1=n1, e_n1=e1, R=n1 / n0, e_R=abs(n1 / n0) * np.hypot(e1 / n1, e0 / n0) if n1 and n0 else np.nan,
             Rlow=nl / n0, e_Rlow=abs(nl / n0) * np.hypot(el / nl, e0 / n0) if nl and n0 else np.nan)
    return r


def chi2_shape(a, b):
    fl = lambda r: np.sqrt(r["e1"] ** 2 + 1.0)   # +1 zdarzenie: podłoga błędu dla pustych binów  # noqa: E731
    ha, ea, hb, eb = a["h1"] / a["h1"].sum(), fl(a) / abs(a["h1"].sum()), b["h1"] / b["h1"].sum(), fl(b) / abs(b["h1"].sum())
    return float(np.sum((ha - hb) ** 2 / (ea ** 2 + eb ** 2))), NM - 1


cols = ["m", "score", "dm", "hypo_from_tag", "pt", "Lxy", "evblk"]
rows, shp, boot, M = [], [], [], {}
t0 = time.time()
for k in SAMPLES:
    d0 = pd.read_parquet(f"{IN}/{k}.parquet", columns=cols)
    wp = pass_wp(d0.score.to_numpy(), d0.pt.to_numpy(), THR["D0"], 0.7)
    SM = {"wp70": wp, "nowp": np.ones(len(d0), bool), "s857": (d0.score > 0.857).to_numpy(), "s387": (d0.score > 0.387).to_numpy()}
    for sel in SELS:
        for lo, hi in PTB:
            d = d0[SM[sel] & ((d0.pt >= lo) & (d0.pt < hi)).to_numpy()]
            P = prep(d)
            masks = {c: (P["L"] >= a) & (P["L"] < b) for c, (a, b) in CUTS.items()}
            for c, mk in masks.items():
                try:
                    r = measure(P, mk)
                except Exception as ex:  # noqa
                    print("FAIL", k, sel, lo, hi, c, ex); continue
                M[sel, lo, hi, k, c] = r
                rows.append(dict(sel=sel, bin=f"pT{lo}-{hi}", sample=k, cut=c, n_tag_win=r["n0"], e_n_tag_win=r["e_n0"], n_refl=r["n1"],
                                 e_n_refl=r["e_n1"], R=r["R"], e_R=r["e_R"], Rlow=r["Rlow"], e_Rlow=r["e_Rlow"], k_sb=r["k_sb"]))
            g = np.random.default_rng(777)
            for ib in range(NBOOT):
                w = g.poisson(1.0, P["nev"])[P["code"]].astype(float)
                o = dict(sel=sel, bin=f"pT{lo}-{hi}", sample=k, ib=ib)
                try:
                    for c, mk in masks.items():
                        r = measure(P, mk, w)
                        o[f"R_{c}"], o[f"Rlow_{c}"], o[f"n0_{c}"], o[f"n1_{c}"] = r["R"], r["Rlow"], r["n0"], r["n1"]
                except Exception:  # noqa
                    continue
                boot.append(o)
        print(k, sel, "%.0fs" % (time.time() - t0), flush=True)
    del d0
T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/refl_rows.csv", index=False)
if boot:
    pd.DataFrame(boot).to_csv(f"{OUT}/refl_boot.csv", index=False)
for (sel, lo, hi, k, c), r in M.items():
    if not k.startswith("data"):
        continue
    for ref, key in (("fd_nd", (sel, lo, hi, "fd_nd", c)), ("mb_nd", (sel, lo, hi, "mb_nd", c)), ("same_data@ge025", (sel, lo, hi, k, "ge025"))):
        if (ref.startswith("same") and c != "slice") or key not in M:
            continue
        c2, nd = chi2_shape(r, M[key])
        shp.append(dict(sel=sel, bin=f"pT{lo}-{hi}", sample=k, cut=c, vs=ref, chi2=c2, ndf=nd, n_refl=r["n1"], e_n_refl=r["e_n1"], n_refl_ref=M[key]["n1"],
                        e_n_refl_ref=M[key]["e_n1"]))
# MC fd vs mb (kształt) i plaster vs ≥ 0,25 w MC
for (sel, lo, hi, k, c), r in M.items():
    if k == "fd_nd" and (sel, lo, hi, "mb_nd", c) in M:
        c2, nd = chi2_shape(M[sel, lo, hi, "mb_nd", c], r)
        shp.append(dict(sel=sel, bin=f"pT{lo}-{hi}", sample="mb_nd", cut=c, vs="fd_nd", chi2=c2, ndf=nd, n_refl=M[sel, lo, hi, "mb_nd", c]["n1"],
                        e_n_refl=M[sel, lo, hi, "mb_nd", c]["e_n1"], n_refl_ref=r["n1"], e_n_refl_ref=r["e_n1"]))
    if k in ("fd_nd", "mb_nd") and c == "slice":
        c2, nd = chi2_shape(r, M[sel, lo, hi, k, "ge025"])
        shp.append(dict(sel=sel, bin=f"pT{lo}-{hi}", sample=k, cut=c, vs="same_mc@ge025", chi2=c2, ndf=nd, n_refl=r["n1"], e_n_refl=r["e_n1"],
                        n_refl_ref=M[sel, lo, hi, k, "ge025"]["n1"], e_n_refl_ref=M[sel, lo, hi, k, "ge025"]["e_n1"]))
pd.DataFrame(shp).to_csv(f"{OUT}/refl_shape.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(T.round(3).to_string())
