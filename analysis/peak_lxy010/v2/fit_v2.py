"""PEAK_lxy010 v2 — pik Δm(D*) i m(Kπ) z tagu D* przy L_xy ≥ 0,25 / ≥ 0,10 / plaster 0,10–0,25 mm, z MC przy TYM SAMYM cięciu L_xy.
Wejścia: ../runs/lxplus/{fd_nd,mb_nd,data_428770,data_427929}.parquet (h5 v2.0, L_xy ≥ 0,10, score D0_H z modeli foldów).
Silnik fitu jak w wersji wstępnej (PEAK_pos: fitlib.fit_mc / fit_data, ../dmfit.fit_dm), selekcja: hypo_from_tag = 0, dstar_conflict = 0,
WP70 D0_H (progi per bin pT) albo bez klasyfikatora; bin po pT(D0).

Użycie: python fit_v2.py <próbka> [nboot]
  fd_nd        kształt sygnału MC (prawda: truth ∈ {0,1}, hypo_swap = 0; 2 Gaussy, wszystko wolne) per sel × bin × cięcie
               + bootstrap po zdarzeniach (evblk) różnicy μ, σ_eff między ≥ 0,10 a ≥ 0,25            → out/mc_fd.csv, out/mc_fd_boot.csv
  mb_nd        (a) sygnał z prawdy dopasowany kształtem fd z tego samego cięcia (μ_mb − μ_fd, k_mb/fd) → out/rows_mb_nd_truth.csv
               (b) wszystkie kandydaty jak dane (pseudo-dane z tłem)                                    → out/rows_mb_nd.csv (+ boot)
  data_*       dane: kształt fd z tego samego cięcia (nominalnie), fd z ≥ 0,25 (jak wersja wstępna)     → out/rows_<próbka>.csv (+ boot)
Bootstrap: Poisson(1) na zdarzenie (evblk), te same repliki dla ≥ 0,10 i ≥ 0,25 → δμ(≥0,10 − ≥0,25), k(≥0,10)/k(≥0,25).
Najpierw uruchomić fd_nd (pozostałe czytają out/mc_fd.csv)."""
import sys, os, time
import numpy as np, pandas as pd
W = "/home/k/Documents/charm/charm_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"; IN = f"{HERE}/../runs/lxplus"
sys.path.insert(0, f"{W}/PEAK_pos"); sys.path.insert(0, f"{HERE}/..")
SAMPLE = sys.argv[1]; NBOOT = int(sys.argv[2]) if len(sys.argv) > 2 else 0
sys.argv = sys.argv[:1]
from fitlib import fit_mc, fit_data                    # noqa: E402
from common import THR, M0, pass_wp                    # noqa: E402
from dmfit import fit_dm, DM0                          # noqa: E402
PTB = [(3, 5), (5, 8), (8, 30), (3, 30)]
CUTS = {"ge025": (0.25, 1e9), "ge010": (0.10, 1e9), "slice": (0.10, 0.25)}
SELS = ("wp70", "nowp")
COLS = ["m", "pt", "dm", "score", "hypo_from_tag", "dstar_conflict", "Lxy", "evblk", "truth", "hypo_swap"]


def load(s, truth=False):
    L = pd.read_parquet(f"{IN}/{s}.parquet", columns=COLS)
    L = L[(L.hypo_from_tag == 0) & (L.dstar_conflict == 0)]
    if truth:
        L = L[np.isin(L.truth, [0, 1]) & (L.hypo_swap == 0)]
    # tylko wiersze używane przez któryś z fitów
    L = L[((np.abs(L.m - M0["D0"]) < 0.04) & (L.dm < 0.160)) | (np.abs(L.dm - DM0) < 0.0012)]
    L["wp70"] = pass_wp(L.score.to_numpy(), L.pt.to_numpy(), THR["D0"], 0.7)
    return L.reset_index(drop=True)


def xdm(d):
    return d.dm.to_numpy()[(np.abs(d.m - M0["D0"]) < 0.04).to_numpy()]


def xm(d):
    return d.m.to_numpy()[(np.abs(d.dm - DM0) < 0.0012).to_numpy()]


def cut(d, c):
    a, b = CUTS[c]
    return d[(d.Lxy >= a) & (d.Lxy < b)]


def binsel(d, sel, lo, hi):
    d = d[(d.pt >= lo) & (d.pt < hi)]
    return d[d.wp70] if sel == "wp70" else d


def mc_shapes(d):
    r1 = fit_dm(xdm(d)); r2 = fit_mc(xm(d), 1.76, 1.98)
    return r1, r2


def like_data(d, shd, shm):
    r1 = fit_dm(xdm(d), shape=shd)
    r2 = fit_data(xm(d), 1.76, 1.98, shm, bkg="cheb2", mu0=M0["D0"])
    return {"dm": r1, "mKpi": r2}


def boot_weights(d, g):
    """Poisson(1) na zdarzenie → indeksy wierszy repliki"""
    code, u = pd.factorize(d.evblk.to_numpy())
    return np.repeat(np.arange(len(d)), g.poisson(1.0, len(u))[code])


t0 = time.time()
if SAMPLE == "fd_nd":
    S = load("fd_nd", truth=True)
    rows, boot = [], []
    for sel in SELS:
        for lo, hi in PTB:
            sb = binsel(S, sel, lo, hi)
            for c in CUTS:
                x = cut(sb, c); r1, r2 = mc_shapes(x)
                for obs, r, n in (("dm", r1, len(xdm(x))), ("mKpi", r2, len(xm(x)))):
                    rows.append(dict(sel=sel, bin=f"pT{lo}-{hi}", cut=c, obs=obs, mu=r["mu"], e_mu=r["e_mu"], s1=r["s1"], s2=r["s2"], f=r["f"],
                                     seff=r["seff"], ok=r["ok"], n=n))
            g = np.random.default_rng(20261118)
            for ib in range(NBOOT):
                xb = sb.iloc[boot_weights(sb, g)]
                try:
                    a = mc_shapes(cut(xb, "ge010")); b = mc_shapes(cut(xb, "ge025"))
                except Exception:   # noqa
                    continue
                for obs, ra, rb in (("dm", a[0], b[0]), ("mKpi", a[1], b[1])):
                    boot.append(dict(sel=sel, bin=f"pT{lo}-{hi}", obs=obs, ib=ib, dmu=ra["mu"] - rb["mu"], rs=ra["seff"] / rb["seff"],
                                     mu10=ra["mu"], mu25=rb["mu"], s10=ra["seff"], s25=rb["seff"], ok=ra["ok"] and rb["ok"]))
            print(sel, lo, hi, "%.0fs" % (time.time() - t0), flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/mc_fd.csv", index=False)
    if boot:
        pd.DataFrame(boot).to_csv(f"{OUT}/mc_fd_boot.csv", index=False)
    sys.exit(0)

MC = pd.read_csv(f"{OUT}/mc_fd.csv").set_index(["sel", "bin", "cut", "obs"])


def shapes(sel, lo, hi, c):
    o = []
    for obs in ("dm", "mKpi"):
        r = MC.loc[(sel, f"pT{lo}-{hi}", c, obs)]
        o.append(dict(mu=r.mu, e_mu=r.e_mu, s1=r.s1, s2=r.s2, f=r.f, seff=r.seff))
    return o


def run(D, tag, nboot):
    rows, boot = [], []
    for sel in SELS:
        for lo, hi in PTB:
            db = binsel(D, sel, lo, hi)
            res = {}
            for c in CUTS:
                x = cut(db, c)
                for var, cs in (("fd_same", c), ("fd_ge025", "ge025")):
                    if var == "fd_ge025" and c == "ge025":
                        continue
                    shd, shm = shapes(sel, lo, hi, cs)
                    o = like_data(x, shd, shm)
                    for obs, r in o.items():
                        sh = shd if obs == "dm" else shm
                        rows.append(dict(sample=tag, sel=sel, bin=f"pT{lo}-{hi}", cut=c, obs=obs, shape=var, mu=r["mu"], e_mu=r["e_mu"], k=r["k"],
                                         e_k=r["e_k"], Ns=r["Ns"], e_Ns=r["e_Ns"], ok=r["ok"], n=r["n"], mu_mc=sh["mu"], e_mu_mc=sh["e_mu"],
                                         seff_mc=sh["seff"]))
            print(tag, sel, lo, hi, "%.0fs" % (time.time() - t0), flush=True)
            g = np.random.default_rng(12345)
            sh10, sh25 = shapes(sel, lo, hi, "ge010"), shapes(sel, lo, hi, "ge025")
            for ib in range(nboot):
                xb = db.iloc[boot_weights(db, g)]
                try:
                    a = like_data(cut(xb, "ge010"), *sh10); b = like_data(cut(xb, "ge025"), *sh25)
                except Exception:   # noqa
                    continue
                for obs in a:
                    boot.append(dict(sample=tag, sel=sel, bin=f"pT{lo}-{hi}", obs=obs, ib=ib, mu10=a[obs]["mu"], mu25=b[obs]["mu"], k10=a[obs]["k"],
                                     k25=b[obs]["k"], Ns10=a[obs]["Ns"], Ns25=b[obs]["Ns"], ok=a[obs]["ok"] and b[obs]["ok"]))
            if nboot:
                print("  boot", tag, sel, lo, hi, "%.0fs" % (time.time() - t0), flush=True)
    T = pd.DataFrame(rows)
    sc = np.where(T.obs == "dm", 1e6, 1e3)
    T["dmu"] = (T.mu - T.mu_mc) * sc; T["e_dmu"] = np.hypot(T.e_mu, T.e_mu_mc) * sc; T["Z"] = T.Ns / T.e_Ns
    T.to_csv(f"{OUT}/rows_{tag}.csv", index=False)
    if boot:
        pd.DataFrame(boot).to_csv(f"{OUT}/boot_{tag}.csv", index=False)
    return T


if SAMPLE == "mb_nd":
    T = run(load("mb_nd", truth=True), "mb_nd_truth", 0)
    T = run(load("mb_nd"), "mb_nd", NBOOT)
else:
    T = run(load(SAMPLE), SAMPLE, NBOOT)
pd.set_option("display.width", 250)
print(T[T["shape"] == "fd_same"][["sample", "sel", "bin", "cut", "obs", "Ns", "e_Ns", "Z", "dmu", "e_dmu", "k", "e_k", "ok"]].round(3).to_string())
