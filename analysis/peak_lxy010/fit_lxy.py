"""PEAK_lxy010: pik Δm(D*) i m(Kπ) z tagu D* przy L_xy ≥ 0,10 mm vs ≥ 0,25 mm (i sam plaster 0,10–0,25), osobno runy.
Silnik fitu = PEAK_pos (fitlib.fit_mc / fit_data importowane; fit_dm = dmfit.py, kopia 1:1 z PEAK_pos/fit_dstar.py, bo tamten plik
wykonuje się przy imporcie). Kształt sygnału (2 Gaussy) ZAWSZE z MC fd_nd przy L_xy ≥ 0,25 (MC z L_xy < 0,25 lokalnie nie ma).
Próbki danych:
  428770_wp70_v12  HYP_d0/runs/lxplus_dstar/data_428770.parquet (v1.2, pełna lumi, L_xy ≥ 0,25)  → KONTROLA liczby z PEAK_pos
  427929_wp70      out/dstar_data_427929_v20.parquet (h5 v2.0, 226/458 plików, L_xy ≥ 0,10; score D0_H z modeli foldów)
  427929_nowp      to samo bez cięcia klasyfikatora
  428770_nowp      AFP_mix/runs/lxplus/data_428770_cands.parquet, wiersze D* (h5 v2.0, 97/97 plików, bez score, bez ładunków;
                   pT = pT(D*), masa/Δm kwantowane 0,1 MeV → dither ±0,05 MeV)
Wynik: out/lxy_rows.csv (+ out/lxy_boot.csv z bootstrapu różnicy ≥0,10 − ≥0,25).
Użycie: python fit_lxy.py [nboot]"""
import sys, os, time
import numpy as np, pandas as pd, pyarrow.parquet as pq
W = "/home/k/Documents/charm/charm_work"
sys.path.insert(0, f"{W}/PEAK_pos")
from fitlib import fit_mc, fit_data                    # noqa: E402
from common import THR, M0, pass_wp                    # noqa: E402
from dmfit import fit_dm, MPI, DM0                     # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"
NBOOT = int(sys.argv[1]) if len(sys.argv) > 1 else 0
PTB = [(3, 5), (5, 8), (8, 30), (3, 30)]
CUTS = {"ge025": (0.25, 1e9), "ge010": (0.10, 1e9), "slice": (0.10, 0.25)}


def tagsel(L, wp):
    ok = (L.hypo_from_tag == 0) & (L.dstar_conflict == 0)
    if wp:
        ok &= pass_wp(L.score.to_numpy(), L.pt.to_numpy(), THR["D0"], wp)
    return ok


COLS = ["m", "pt", "eta", "dm", "score", "hypo_from_tag", "dstar_conflict", "Lxy", "truth", "hypo_swap"]
S0 = pd.read_parquet(f"{W}/HYP_d0/runs/lxplus_dstar/fd_nd.parquet", columns=COLS)
S0 = S0[np.isin(S0.truth, [0, 1]) & (S0.hypo_swap == 0)]
MC = {"wp70": S0[tagsel(S0, 0.7)].reset_index(drop=True), "nowp": S0[tagsel(S0, None)].reset_index(drop=True)}
D = {}
d = pd.read_parquet(f"{W}/HYP_d0/runs/lxplus_dstar/data_428770.parquet", columns=COLS[:8])
D["428770_wp70_v12"] = ("wp70", d[tagsel(d, 0.7)].reset_index(drop=True))
d = pd.read_parquet(f"{OUT}/dstar_data_427929_v20.parquet", columns=COLS[:8])
D["427929_wp70"] = ("wp70", d[tagsel(d, 0.7)].reset_index(drop=True))
D["427929_nowp"] = ("nowp", d[tagsel(d, None)].reset_index(drop=True))
t = pq.read_table(f"{W}/AFP_mix/runs/lxplus/data_428770_cands.parquet",
                  filters=[("Dmeson_species", "=", 1)]).to_pandas()
rng = np.random.default_rng(428770)
dm = t.Dmeson_delta_mass.to_numpy(float) + rng.uniform(-5e-5, 5e-5, len(t))
mk = t.Dmeson_mass.to_numpy(float) - t.Dmeson_delta_mass.to_numpy(float) + rng.uniform(-5e-5, 5e-5, len(t))
d = pd.DataFrame({"m": mk, "pt": t.Dmeson_pt.to_numpy(float), "eta": t.Dmeson_eta.to_numpy(float), "dm": dm, "Lxy": t.Dmeson_Lxy.to_numpy(float)})
D["428770_nowp"] = ("nowp", d[(d.pt < 30) & (d.m > 1.75) & (d.m < 2.0)].reset_index(drop=True))
del t, d
for k, (s, d) in D.items():
    print(k, len(d), "Lxy min %.3f" % d.Lxy.min(), "frac<0.25 %.3f" % (d.Lxy < 0.25).mean(), flush=True)

# kształty MC (L_xy ≥ 0,25) per selekcja × bin
SH = {}
for s, S in MC.items():
    for lo, hi in PTB:
        b = ((S.pt >= lo) & (S.pt < hi)).to_numpy()
        SH[s, lo, hi] = (fit_dm(S.dm.to_numpy()[b & (np.abs(S.m - M0["D0"]).to_numpy() < 0.04)]),
                         fit_mc(S.m.to_numpy()[b & (np.abs(S.dm - DM0).to_numpy() < 0.0012)], 1.76, 1.98),
                         int(b.sum()))


def fits(d, s, lo, hi):
    """zwraca dict wyników dla jednej próbki danych (już po cięciu L_xy) w binie pT"""
    b = ((d.pt >= lo) & (d.pt < hi)).to_numpy()
    shd, shm, nmc = SH[s, lo, hi]
    r1 = fit_dm(d.dm.to_numpy()[b & (np.abs(d.m - M0["D0"]).to_numpy() < 0.04)], shape=shd)
    r2 = fit_data(d.m.to_numpy()[b & (np.abs(d.dm - DM0).to_numpy() < 0.0012)], 1.76, 1.98, shm, bkg="cheb2", mu0=M0["D0"])
    o = {}
    for nm, r, sh in (("dm", r1, shd), ("mKpi", r2, shm)):
        o[nm] = dict(mu=r["mu"], e_mu=r["e_mu"], k=r["k"], e_k=r["e_k"], Ns=r["Ns"], e_Ns=r["e_Ns"], ok=r["ok"], n=r["n"],
                     mu_mc=sh["mu"], e_mu_mc=sh["e_mu"], seff_mc=sh["seff"], n_mc=nmc)
    return o


rows, boot = [], []
t0 = time.time()
for name, (s, d) in D.items():
    for lo, hi in PTB:
        for cut, (a, bb) in CUTS.items():
            if d.Lxy.min() > 0.2 and cut != "ge025":
                continue
            o = fits(d[(d.Lxy >= a) & (d.Lxy < bb)], s, lo, hi)
            for obs, r in o.items():
                rows.append(dict(sample=name, sel=s, bin=f"pT{lo}-{hi}", cut=cut, obs=obs, **r))
        print(name, lo, hi, "%.0fs" % (time.time() - t0), flush=True)
        if NBOOT and d.Lxy.min() < 0.2:
            g = np.random.default_rng(12345)
            dd = d[(d.pt >= lo) & (d.pt < hi)].reset_index(drop=True)
            for ib in range(NBOOT):
                x = dd.iloc[g.integers(0, len(dd), len(dd))]
                try:
                    oa, ob = fits(x, s, lo, hi), fits(x[x.Lxy >= 0.25], s, lo, hi)
                except Exception as ex:  # noqa
                    continue
                for obs in oa:
                    boot.append(dict(sample=name, bin=f"pT{lo}-{hi}", obs=obs, ib=ib, dmu=oa[obs]["mu"] - ob[obs]["mu"],
                                     rk=oa[obs]["k"] / ob[obs]["k"], mu10=oa[obs]["mu"], k10=oa[obs]["k"],
                                     ok=oa[obs]["ok"] and ob[obs]["ok"]))
            print("  boot", name, lo, hi, "%.0fs" % (time.time() - t0), flush=True)
T = pd.DataFrame(rows)
sc = np.where(T.obs == "dm", 1e6, 1e3)   # keV dla Δm, MeV dla m(Kπ)
T["dmu"] = (T.mu - T.mu_mc) * sc; T["e_dmu"] = np.hypot(T.e_mu, T.e_mu_mc) * sc; T["Z"] = T.Ns / T.e_Ns
T.to_csv(f"{OUT}/lxy_rows.csv", index=False)
if boot:
    pd.DataFrame(boot).to_csv(f"{OUT}/lxy_boot.csv", index=False)
pd.set_option("display.width", 250)
print(T[["sample", "bin", "cut", "obs", "Ns", "e_Ns", "Z", "dmu", "e_dmu", "k", "e_k", "ok", "n_mc"]].round(3).to_string())
