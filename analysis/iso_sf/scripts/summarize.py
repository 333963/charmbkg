"""ISO_sf: tabele SF z out/sf_fits.json + out/mc_ntrk.json → out/sf_table.md, out/sf_summary.json"""
import json, os
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.environ.get("ISO_OUT", os.path.join(BASE, "out"))   # 16.11 ISO-FULL; mc_ntrk.json (MC v9) zawsze z out/
F = json.load(open(os.path.join(OUT, "sf_fits.json"))); R, NTB = F["res"], F["NTB"]
MN = json.load(open(os.path.join(BASE, "out", "mc_ntrk.json")))
rng = np.random.default_rng(1)


def sd(r):
    if r.get("eps_lo") is not None and np.isfinite(r["eps_lo"]) and np.isfinite(r["eps_hi"]):
        return 0.5 * (r["eps_hi"] - r["eps_lo"])
    return r["eps_err"]


def fracs(ps):
    keys = [f"{ps}|3-30|frac_n_trkid[{a},{b})" for a, b in zip(NTB[:-1], NTB[1:])]
    fd = np.array([R[k]["eps"] for k in keys]); ed = np.array([sd(R[k]) for k in keys]); fm = np.array([R[k]["eps_mc"] for k in keys])
    return fd, ed, fm


def rw_eps(ps, epsbin, fd, ed):
    e = np.array([np.nan if x is None else x for x in epsbin]); ok = np.isfinite(e)
    def one(f):
        f = np.clip(f, 0, None)[ok]; return float((f * e[ok]).sum() / f.sum())
    toys = [one(fd + ed * rng.standard_normal(len(fd))) for _ in range(2000)]
    return one(fd), float(np.std(toys))


rows, summ = [], {}
for ps in ("A_score0.6", "B_none"):
    fd, ed, fm = fracs(ps)
    summ[f"{ps}|ntrkid_fracs"] = {"edges": NTB, "data": fd.tolist(), "data_err": ed.tolist(), "mc": fm.tolist(),
                                  "mean_bin_data": float(np.dot(fd / fd.sum(), np.arange(len(fd)))), "mean_bin_mc": float(np.dot(fm, np.arange(len(fm))))}
    for k, r in R.items():
        p = k.split("|")
        if p[0] != ps or "frac_" in k or "ntrkid_" in k:
            continue
        pb, cut = p[1], p[2]
        ed_, em = sd(r), r["eps_mc"]; dem = np.sqrt(em * (1 - em) / r["n_mc"])
        sf = r["eps"] / em; dsf = sf * np.hypot(ed_ / max(r["eps"], 1e-3), dem / em)
        erw, derw = rw_eps(ps, r["eps_mc_by_ntrkid"], fd, ed)
        sf_rw = r["eps"] / erw; dsf_rw = sf_rw * np.hypot(ed_ / max(r["eps"], 1e-3), derw / erw)
        v9 = MN["dstar"]["cuts"].get(f"{pb}|{cut}")
        ratio_w1 = v9["w1"][0] / v9["unw"][0] if v9 and v9["unw"][0] > 0 else np.nan
        e_w1 = em * ratio_w1
        out = dict(ps=ps, pt=pb, cut=cut, eps_data=r["eps"], err_data=ed_, eps_expo=r["eps_expo"], N=r["N_tot"], dN=r["N_tot_err"],
                   det=r["determined"], eps_mc=em, sf=sf, dsf=dsf, eps_mc_rw=erw, derw=derw, sf_rw=sf_rw, dsf_rw=dsf_rw,
                   eps_mc_w1=e_w1, sf_w1=r["eps"] / e_w1 if e_w1 > 0 else np.nan,
                   dsf_w1=(r["eps"] / e_w1) * np.hypot(ed_ / max(r["eps"], 1e-3), v9["w1"][1] / v9["w1"][0]) if e_w1 > 0 else np.nan)
        rows.append(out); summ[k] = out


def order(r):
    v, c = r["cut"].split("<"); return (["ptcone40", "ptcone100", "rel40", "rel100"].index(v), float(c))


L = []
for ps, pbs in (("A_score0.6", ("3-30", "3-5", "5-8", "8-30")), ("B_none", ("3-30",))):
    for pb in pbs:
        rr = sorted([r for r in rows if r["ps"] == ps and r["pt"] == pb], key=order)
        if not rr:
            continue
        L.append(f"\n### {ps}, pT {pb} GeV (N_sygn = {rr[0]['N']:.0f} ± {rr[0]['dN']:.0f})\n")
        L.append("| cięcie | ε dane | ε MC | **SF** | ε MC rw(n_trkid) | SF rw | ε MC w1(n_trk) | SF w1 | det |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for r in rr:
            L.append(f"| {r['cut']} | {r['eps_data']:.3f} ± {r['err_data']:.3f} | {r['eps_mc']:.3f} | **{r['sf']:.2f} ± {r['dsf']:.2f}** | "
                     f"{r['eps_mc_rw']:.3f} ± {r['derw']:.3f} | {r['sf_rw']:.2f} ± {r['dsf_rw']:.2f} | {r['eps_mc_w1']:.3f} | {r['sf_w1']:.2f} ± {r['dsf_w1']:.2f} | {'T' if r['det'] else 'N'} |")
# krotność
L.append("\n### Rozkład n_trkid sygnału D* (frakcje z fitu pass/fail „w binie”)\n")
L.append("| presel | bin n_trkid | f dane | f MC |"); L.append("|---|---|---|---|")
for ps in ("A_score0.6", "B_none"):
    s = summ[f"{ps}|ntrkid_fracs"]
    for i, (a, b) in enumerate(zip(NTB[:-1], NTB[1:])):
        L.append(f"| {ps} | [{a},{b}) | {s['data'][i]:.3f} ± {s['data_err'][i]:.3f} | {s['mc'][i]:.3f} |")
L.append("\n### SF w podziale na krotność (n_trkid < 25 / ≥ 25), pT 3–30\n")
L.append("| presel | cięcie | n_trkid | ε dane | ε MC | SF | N_sygn |"); L.append("|---|---|---|---|---|---|---|")
for k, r in sorted(R.items()):
    if "ntrkid_" not in k:
        continue
    ps, pb, tg, cut = k.split("|")
    em = r["eps_mc"]; e = r["eps"]; s_ = sd(r)
    dsf = (e / em) * np.hypot(s_ / max(e, 1e-3), np.sqrt(em * (1 - em) / r["n_mc"]) / em)
    L.append(f"| {ps} | {cut} | {tg} | {e:.3f} ± {s_:.3f} | {em:.3f} | {e / em:.2f} ± {dsf:.2f} | {r['N_tot']:.0f} ± {r['N_tot_err']:.0f} |")
    summ[k] = dict(r, sf=e / em, dsf=dsf)
open(os.path.join(OUT, "sf_table.md"), "w").write("\n".join(L) + "\n")
json.dump(summ, open(os.path.join(OUT, "sf_summary.json"), "w"), indent=1, default=float)
print("\n".join(L))
