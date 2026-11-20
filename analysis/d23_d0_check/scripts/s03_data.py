"""Testy na danych 2-3 GeV: C (eps_B(m)), D (kategorie rozlaczne drzewo/siec), F (ksztalt wolny), G (wycinki score, warianty fitu, LR), I (pT pod-biny, stosunek runow)."""
import numpy as np, json, math, time
from scipy.stats import chi2 as CH, norm
from vlib import *
out = {}
F = load("fd_d23")
D = {r: load(f"data_{r}_d23") for r in ("427929", "428770")}
THR = {(j, wp): mc_wp(F, j, wp) for j in (0, 1) for wp in (0.5, 0.4, 0.3, 0.2, 0.1)}
HEAD = {("427929", 0): 0.1, ("427929", 1): 0.2, ("428770", 0): 0.1, ("428770", 1): 0.1}     # punkty z werdyktu etapu 4
t0 = time.time()

print("== G1. warianty fitu przy punktach z werdyktu (fit2 = niezalezna implementacja); Z_LR = sqrt(Δ(-2lnL)) wzgl. N=0")
for (run, j), wp in HEAD.items():
    d = D[run]; th = THR[(j, wp)]; q = templ(F, j, th); s = d["lg"][:, j] >= th
    key = f"{run}|{NM[j]}|{wp}"; O = out.setdefault("var|" + key, {})
    print(f" {run} {NM[j]} WP{int(wp*100)} (eps {q['eps']:.4f}, frefl MC {q['frefl']:.3f}, fwin {q['fwin']:.3f})")
    for lab, rng, nb, kind in (("nom exp2 1.75-2.00", (1.75, 2.0), 2, "exp"), ("exp1", (1.75, 2.0), 1, "exp"), ("exp3", (1.75, 2.0), 3, "exp"), ("exp4", (1.75, 2.0), 4, "exp"), ("pol2", (1.75, 2.0), 2, "pol"),
                               ("pol3", (1.75, 2.0), 3, "pol"), ("exp2 1.70-2.05", (1.70, 2.05), 2, "exp"), ("exp3 1.70-2.05", (1.70, 2.05), 3, "exp"), ("exp4 1.70-2.05", (1.70, 2.05), 4, "exp"),
                               ("exp2 1.78-1.96", (1.78, 1.96), 2, "exp"), ("exp1 1.78-1.96", (1.78, 1.96), 1, "exp")):
        H, e = hist(d["m"][s], rng)
        f = fit2(H, e, q, nbkg=nb, kind=kind); f0 = fit2(H, e, q, nbkg=nb, kind=kind, fixN=0.0)
        zlr = math.sqrt(max(f0["fval"] - f["fval"], 0))
        O[lab] = dict(N=f["N"], err=f["err"], Z=f["N"] / f["err"], ZLR=zlr, chi2=f["chi2"], ndf=f["ndf"], chi2_0=f0["chi2"], B=f["B"], valid=f["valid"])
        print(f"   {lab:20s} N {f['N']:7.1f} ± {f['err']:5.1f}  N/σ {f['N']/f['err']:5.2f}  Z_LR {zlr:5.2f}  chi2 {f['chi2']:5.1f}/{f['ndf']} (bez sygnalu {f0['chi2']:.1f})  B_okno {f['B']:.0f} {'' if f['valid'] else 'NIEWAZNY'}")
    H, e = hist(d["m"][s])
    for lab, kw in (("frefl=0", dict(frefl=0.0)), ("frefl x0.5", dict(frefl=0.5 * q["frefl"])), ("frefl x1.5", dict(frefl=1.5 * q["frefl"])), ("frefl wolne", dict(frefl="free")),
                    ("frefl wolne 1.70-2.05 exp3", dict(frefl="free", nbkg=3, rng=(1.70, 2.05)))):
        rng = kw.pop("rng", (1.75, 2.0)); H_, e_ = hist(d["m"][s], rng)
        f = fit2(H_, e_, q, **kw)
        O[lab] = dict(N=f["N"], err=f["err"], fr=f["fr"], fr_err=f["fr_err"], chi2=f["chi2"], ndf=f["ndf"], valid=f["valid"])
        print(f"   {lab:26s} N {f['N']:7.1f} ± {f['err']:5.1f}  frefl {f['fr']:.2f} ± {f['fr_err']:.2f}  chi2 {f['chi2']:.1f}/{f['ndf']}")
    # ksztalt wolny: szablon (dmu, sc) i Gauss; bledy z bootstrapu po zdarzeniach
    sr = s & (d["m"] > 1.75) & (d["m"] < 2.0)
    ft = fit2(H, e, q, free_shape=True)
    bs, nok = boot2(d["m"][sr], d["uid"][sr], q, nb=100, keys=("N", "a", "b"), free_shape=True)
    b0, nok0 = boot2(d["m"][sr], d["uid"][sr], q, nb=100, keys=("N",))
    # Gauss: odniesienie = ten sam fit Gaussa na MC (sygnal+odbicie MC, bez tla -> tlo ~0)
    sm = (F["cat"] == 0) & (F["lg"][:, j] >= th)
    mu_mc, sg_mc = norm.fit(F["m"][sm][np.abs(F["m"][sm] - M0) < 0.045])
    Hm, em = hist(F["m"][(F["cat"] <= 1) & (F["lg"][:, j] >= th)])
    fm = fit2(Hm, em, q, gauss=(M0, 0.02), nbkg=1)
    fg = fit2(H, e, q, gauss=(M0, 0.02))
    bg, nokg = boot2(d["m"][sr], d["uid"][sr], q, nb=100, keys=("N", "a", "b"), gauss=(M0, 0.02))
    O["free"] = dict(N=ft["N"], N_bs=bs[0], dmu=ft["a"], dmu_hesse=ft["a_err"], dmu_bs=bs[1], sc=ft["b"], sc_hesse=ft["b_err"], sc_bs=bs[2], nok=nok, N_nom_bs=b0[0],
                     g_N=fg["N"], g_N_bs=bg[0], g_mu=fg["a"], g_mu_bs=bg[1], g_sig=fg["b"], g_sig_bs=bg[2], g_mu_mc=fm["a"], g_sig_mc=fm["b"], g_sig_mc_err=fm["b_err"], nokg=nokg, chi2_g=fg["chi2"], ndf_g=fg["ndf"])
    print(f"   bootstrap po zdarzeniach (fit nominalny): σ(N) = {b0[0]:.1f} ({nok0} replik)")
    print(f"   SZABLON wolny: N {ft['N']:.0f} ± {bs[0]:.0f}; dmu {1e3*ft['a']:+.2f} ± {1e3*bs[1]:.2f} MeV (HESSE {1e3*ft['a_err']:.2f}); sc {ft['b']:.3f} ± {bs[2]:.3f} (HESSE {ft['b_err']:.3f}); {nok} replik")
    print(f"   GAUSS wolny:   N {fg['N']:.0f} ± {bg[0]:.0f}; mu {1e3*(fg['a']-M0):+.2f} ± {1e3*bg[1]:.2f} MeV wzgl. PDG (MC: {1e3*(fm['a']-M0):+.2f}); sigma {1e3*fg['b']:.1f} ± {1e3*bg[2]:.1f} MeV (MC ten sam fit: {1e3*fm['b']:.1f} ± {1e3*fm['b_err']:.1f}) -> dane/MC {fg['b']/fm['b']:.2f} ± {bg[2]/fm['b']:.2f}; chi2 {fg['chi2']:.0f}/{fg['ndf']}")
print(f"[{time.time()-t0:.0f}s]")

print("== G2/test 4. wycinki rozlaczne score: y_i = N_i/Δeps_i (statystycznie niezalezne), chi2 stalej")
EDG = [0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
for run in D:
    d = D[run]
    for j in (0, 1):
        ys, es, txt = [], [], []
        for hi_wp, lo_wp in zip(EDG[:-1], EDG[1:]):
            lo = THR[(j, hi_wp)]; hi = THR[(j, lo_wp)] if lo_wp > 0 else np.inf
            q = templ(F, j, lo, hi); s = (d["lg"][:, j] >= lo) & (d["lg"][:, j] < hi)
            H, e = hist(d["m"][s]); f = fit2(H, e, q)
            sr = s & (d["m"] > 1.75) & (d["m"] < 2.0)
            bs, _ = boot2(d["m"][sr], d["uid"][sr], q, nb=40)
            er = max(f["err"], bs[0])
            ys.append(f["N"] / q["eps"]); es.append(er / q["eps"])
            txt.append(f"{hi_wp:.1f}-{lo_wp:.1f}: N {f['N']:.0f} ± {er:.0f}, S/B {f['N']*q['fwin']/f['B']:.3f}, y {ys[-1]:.0f} ± {es[-1]:.0f}")
            out[f"slice|{run}|{NM[j]}|{hi_wp}-{lo_wp}"] = dict(N=f["N"], err=er, eps=q["eps"], y=ys[-1], yerr=es[-1], B=f["B"], frefl=q["frefl"], chi2=f["chi2"])
        ys, es = np.array(ys), np.array(es)
        for nm, sl in (("5 wycinkow (eps<0.5)", slice(0, 5)), ("3 ciasne (eps<0.3)", slice(2, 5))):
            w = 1 / es[sl] ** 2; mean = (ys[sl] * w).sum() / w.sum(); c2 = float((((ys[sl] - mean) / es[sl]) ** 2).sum()); nd = len(ys[sl]) - 1
            out[f"slicefit|{run}|{NM[j]}|{nm}"] = dict(mean=mean, err=1 / math.sqrt(w.sum()), chi2=c2, ndf=nd, p=float(CH.sf(c2, nd)))
            print(f"  {run} {NM[j]} {nm}: <y> = {mean:.0f} ± {1/math.sqrt(w.sum()):.0f}, chi2/ndf = {c2:.1f}/{nd}, p = {CH.sf(c2, nd):.2f}")
        for t in txt: print("      ", t)
print(f"[{time.time()-t0:.0f}s]")

print("== D/test 1d. kategorie rozlaczne przy WP20 obu klasyfikatorow i przy punktach z werdyktu: b = oba, t = tylko drzewo, n = tylko siec")
for run in D:
    d = D[run]
    for wt, wn in ((0.2, 0.2), (0.1, 0.1), (HEAD[(run, 0)], HEAD[(run, 1)])):
        tt, tn = THR[(0, wt)], THR[(1, wn)]
        sig = F["cat"] == 0; rf = F["cat"] == 1; nsig = sig.sum()
        res = {}
        for c, (mf, md) in dict(b=((F["lg"][:, 0] >= tt) & (F["lg"][:, 1] >= tn), (d["lg"][:, 0] >= tt) & (d["lg"][:, 1] >= tn)),
                                t=((F["lg"][:, 0] >= tt) & (F["lg"][:, 1] < tn), (d["lg"][:, 0] >= tt) & (d["lg"][:, 1] < tn)),
                                n=((F["lg"][:, 0] < tt) & (F["lg"][:, 1] >= tn), (d["lg"][:, 0] < tt) & (d["lg"][:, 1] >= tn))).items():
            ps, pr = sig & mf, rf & mf
            q = dict(sigT=A.Tmpl(F["m"][ps], np.ones(ps.sum()), *MCR), refT=A.Tmpl(F["m"][pr], np.ones(pr.sum()), *MCR), frefl=pr.sum() / ps.sum(), eps=ps.sum() / nsig)
            H, e = hist(d["m"][md]); f = fit2(H, e, q)
            sr = md & (d["m"] > 1.75) & (d["m"] < 2.0); bs, _ = boot2(d["m"][sr], d["uid"][sr], q, nb=40)
            res[c] = (f["N"], max(f["err"], bs[0]), q["eps"], ps.sum())
        Nb, Nt, Nn = (res[c][0] for c in "btn"); sb_, st_, sn_ = (res[c][1] for c in "btn"); eb, et, en = (res[c][2] for c in "btn")
        R = ((Nb + Nn) / (eb + en)) / ((Nb + Nt) / (eb + et))
        dR = R * math.sqrt((sb_ * (1 / (Nb + Nn) - 1 / (Nb + Nt))) ** 2 + (sn_ / (Nb + Nn)) ** 2 + (st_ / (Nb + Nt)) ** 2)
        # roznica y miedzy kategoria 'tylko drzewo' a 'tylko siec' i 'oba'
        y = {c: (res[c][0] / res[c][2], res[c][1] / res[c][2]) for c in "btn"}
        out[f"cat|{run}|{wt}|{wn}"] = dict(R=R, R_err=dR, **{f"N_{c}": res[c][0] for c in "btn"}, **{f"e_{c}": res[c][2] for c in "btn"}, **{f"err_{c}": res[c][1] for c in "btn"}, **{f"y_{c}": y[c][0] for c in "btn"}, **{f"yerr_{c}": y[c][1] for c in "btn"})
        print(f"  {run} drzewo WP{int(wt*100)} / siec WP{int(wn*100)}: " + "; ".join(f"{c}: N {res[c][0]:.0f} ± {res[c][1]:.0f}, e {res[c][2]:.4f}, y {y[c][0]:.0f} ± {y[c][1]:.0f}" for c in "btn") + f" | R = y_siec/y_drzewo = {R:.3f} ± {dR:.3f}")
print(f"[{time.time()-t0:.0f}s]")

print("== I. pod-biny pT (2-2.5 / 2.5-3): y i pozycja piku; MC: udzial sygnalu po cieciu")
for (run, j), wp in HEAD.items():
    d = D[run]; th = THR[(j, wp)]
    for lo, hi in ((2.0, 2.5), (2.5, 3.0)):
        fm = (F["pt"] >= lo) & (F["pt"] < hi); Fs = {k: (v[fm] if hasattr(v, "shape") and v.shape[:1] == fm.shape else v) for k, v in F.items()}
        q = templ(Fs, j, th); nmc = ((Fs["cat"] == 0) & (Fs["lg"][:, j] >= th)).sum()
        s = (d["lg"][:, j] >= th) & (d["pt"] >= lo) & (d["pt"] < hi)
        H, e = hist(d["m"][s]); f = fit2(H, e, q); fg = fit2(H, e, q, free_shape=True)
        out[f"ptsub|{run}|{NM[j]}|{lo}"] = dict(N=f["N"], err=f["err"], nmc=int(nmc), dmu=fg["a"], dmu_err=fg["a_err"], sc=fg["b"], B=f["B"])
        print(f"  {run} {NM[j]} WP{int(wp*100)} pT {lo}-{hi}: N {f['N']:.0f} ± {f['err']:.0f} (Z {f['N']/f['err']:.1f}), S/B {f['N']*q['fwin']/f['B']:.2f}, sygnal MC po cieciu {nmc}, N/N_MC {f['N']/nmc:.4f} ± {f['err']/nmc:.4f}; dmu {1e3*fg['a']:+.1f} ± {1e3*fg['a_err']:.1f} MeV, sc {fg['b']:.2f}")
json.dump(out, open(os.path.join(BASE, "out", "s03.json"), "w"), indent=1, default=float)
print(f"[{time.time()-t0:.0f}s]")
