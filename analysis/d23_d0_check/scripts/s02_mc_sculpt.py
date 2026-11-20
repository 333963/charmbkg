"""A: ksztalty MC. B (test 1a): tlo MC z prawda (mb_nd, 3-5 GeV — w 2-3 GeV MC minbias nie istnieje) po ciasnym cieciu. E (test 2): blizniaki K<->pi."""
import numpy as np, json
from vlib import *
out = {}
def seff(m, w=None):
    lo, hi = np.quantile(m, [0.1587, 0.8413]); return 0.5 * (hi - lo)
F23, F35, MB = load("fd_d23"), load("fd_nd_d150"), load("mb_nd_d150")
print("== A. ksztalt sygnalu MC (fd_nd), sigma_eff = polowa przedzialu 68 % [MeV], mediana - m0 [MeV]; odbicie: udzial w zakresie fitu 1.75-2.00 i w oknie ±55 wzgl. sygnalu")
for nm, F in (("2-3", F23), ("3-5", F35)):
    for j in (0, 1):
        for wp in (1.0, 0.7, 0.3, 0.2, 0.1):
            th = -99 if wp == 1.0 else mc_wp(F, j, wp)
            s = (F["cat"] == 0) & (F["lg"][:, j] >= th); r = (F["cat"] == 1) & (F["lg"][:, j] >= th)
            ms, mr = F["m"][s], F["m"][r]
            core = np.abs(ms - M0) < 0.08
            o = dict(n=int(s.sum()), eps=s.sum() / (F["cat"] == 0).sum(), seff=1e3 * seff(ms[core]), med=1e3 * (np.median(ms[core]) - M0), frefl=r.sum() / s.sum(),
                     refl_fit=((mr > 1.75) & (mr < 2.0)).sum() / s.sum(), refl_win=(np.abs(mr - M0) < WIN).sum() / s.sum(), refl_seff=1e3 * seff(mr), refl_med=1e3 * (np.median(mr) - M0))
            out[f"shape|{nm}|{NM[j]}|{wp}"] = o
            print(f"  {nm} {NM[j]} WP{int(wp*100):3d}: n {o['n']:6d} eps {o['eps']:.3f} sigma_eff {o['seff']:.1f} med {o['med']:+.1f} | odbicie/sygnal (1.70-2.05) {o['frefl']:.3f}, w zakresie fitu {o['refl_fit']:.3f}, w oknie {o['refl_win']:.3f}, sigma_eff odb. {o['refl_seff']:.0f}, med {o['refl_med']:+.0f}")

print("== E. blizniaki (test 2), MC 2-3 GeV: sygnal (poprawna hipoteza) po cieciu — czy blizniak (zamieniona hipoteza) jest w tabeli (m_swap w 1.70-2.05) i czy tez przechodzi")
for j in (0, 1):
    for wp in (0.5, 0.2, 0.1):
        th = mc_wp(F23, j, wp)
        s = (F23["cat"] == 0) & (F23["lg"][:, j] >= th)
        has = F23["has"][s]; twp = F23["tw"][s, j] >= th
        o = dict(n=int(s.sum()), twin_in_table=has.mean(), twin_pass=(has & twp).mean())
        r = (F23["cat"] == 1) & (F23["lg"][:, j] >= th)
        o["swap_pass_with_correct_pass"] = float(((F23["tw"][r, j] >= th) & F23["has"][r]).mean())
        out[f"twinMC|{NM[j]}|{wp}"] = o
        print(f"  {NM[j]} WP{int(wp*100)}: sygnal {o['n']}, blizniak w tabeli {o['twin_in_table']:.3f}, blizniak tez przechodzi {o['twin_pass']:.3f}; zamienione przechodzace: poprawny tez przechodzi {o['swap_pass_with_correct_pass']:.3f}")
print("   dane 2-3 GeV: wybrane kandydaty w zakresie fitu 1.75-2.00")
for run in ("427929", "428770"):
    D = load(f"data_{run}_d23")
    for j in (0, 1):
        for wp in (0.2, 0.1):
            th = mc_wp(F23, j, wp)
            s = (D["lg"][:, j] >= th) & (D["m"] > 1.75) & (D["m"] < 2.0)
            u, c = np.unique(D["uid"][s], return_counts=True)
            twp = D["has"][s] & (D["tw"][s, j] >= th)
            twfit = twp & (D["ms"][s] > 1.75) & (D["ms"][s] < 2.0)
            inw = np.abs(D["m"][s] - M0) < WIN
            o = dict(n=int(s.sum()), nev=len(u), mult=c.mean(), frac_ev_multi=(c > 1).mean(), twin_pass=twp.mean(), twin_pass_infit=twfit.mean(),
                     win_twin_pass=float(twp[inw].mean()), win_twin_in_win=float((twp & (np.abs(D["ms"][s] - M0) < WIN))[inw].mean()))
            out[f"twinData|{run}|{NM[j]}|{wp}"] = o
            print(f"  {run} {NM[j]} WP{int(wp*100)}: {o['n']} kand. / {o['nev']} zdarzen = {o['mult']:.3f}; zdarzen z >1: {o['frac_ev_multi']:.3f}; blizniak tez przechodzi {o['twin_pass']:.3f} (w zakresie fitu {o['twin_pass_infit']:.3f}); w oknie: {o['win_twin_pass']:.3f}, oba w oknie {o['win_twin_in_win']:.3f}")

print("== B. test 1a: tlo MC minbias (mb_nd, pT 3-5 GeV; wagi w) po cieciu; progi z fd_nd 3-5; fit sygnal+odbicie+tlo na SAMYM tle")
print("   wagi mb_nd:", np.unique(MB["w"])[:6], "kategorie:", {int(c): int((MB['cat'] == c).sum()) for c in range(4)})
D35 = {r: load(f"data_{r}_d150") for r in ("427929", "428770")}
for j in (0, 1):
    for wp in (0.7, 0.3, 0.2, 0.1):
        th = mc_wp(F35, j, wp); q = templ(F35, j, th)
        if th < (1.0, 1.4)[j]:
            print(f"  {NM[j]} WP{int(wp*100)}: prog {th:.2f} ponizej luznego ciecia redukcji — pomijam"); continue
        for lab, cats in (("truth>=8", (3,)), ("truth 2-7", (2,)), ("cale tlo", (2, 3))):
            b = np.isin(MB["cat"], cats) & (MB["lg"][:, j] >= th)
            H, e = hist(MB["m"][b], w=MB["w"][b].astype(float))
            H1, e1 = hist(MB["m"][b])
            f = fit2(H1, e1, q)        # fit na liczbach nieważonych (wagi ~1?) 
            ctr = 0.5 * (e[1:] + e[:-1]); wn = np.abs(ctr - M0) < WIN
            nwin, nsb = H1[wn].sum(), H1[~wn].sum()
            # R_sculpt: okno / interpolacja liniowa z pasm (fit pol1 exp do pasm)
            o = dict(n=int(b.sum()), nwin=float(nwin), Nfake=f["N"], err=f["err"], Bwin=f["B"], frac=f["N"] * q["fwin"] / max(f["B"], 1e-9), frac_err=f["err"] * q["fwin"] / max(f["B"], 1e-9), chi2=f["chi2"], ndf=f["ndf"])
            out[f"mbbkg|{NM[j]}|{wp}|{lab}"] = o
            print(f"  {NM[j]} WP{int(wp*100)} {lab:9s}: n(1.75-2.00) {int(H1.sum()):6d}, w oknie {int(nwin):5d}; N_falszywe = {f['N']:+.1f} ± {f['err']:.1f} (chi2 {f['chi2']:.0f}/{f['ndf']}); N_f*fwin/B_okno = {o['frac']:+.3f} ± {o['frac_err']:.3f}")
        # domkniecie z pelnym tlem (wstrzykniecie): mb_nd wszystko
        allm = MB["lg"][:, j] >= th
        H1, e1 = hist(MB["m"][allm]); f = fit2(H1, e1, q)
        ntrue = ((MB["cat"] == 0) & allm).sum()
        out[f"mbclos|{NM[j]}|{wp}"] = dict(N=f["N"], err=f["err"], ntrue=int(ntrue), ratio=f["N"] / ntrue, ratio_err=f["err"] / ntrue)
        print(f"     domkniecie mb_nd (sygnal+tlo): N_fit {f['N']:.0f} ± {f['err']:.0f} vs prawda {ntrue} -> {f['N']/ntrue:.3f} ± {f['err']/ntrue:.3f}; S/B MC {ntrue*q['fwin']/f['B']:.2f}")
json.dump(out, open(os.path.join(BASE, "out", "s02.json"), "w"), indent=1, default=float)
