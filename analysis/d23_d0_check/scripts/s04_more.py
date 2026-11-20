"""C (test 1c): eps_B(m) w danych; 3-5 GeV: przenoszenie eps przy ciasnym cieciu; stosunek runow vs tlo (test 'pik ∝ tlo'); test 5 (wiarygodnosc); profil sc; liczba prob."""
import numpy as np, json, math
from scipy.stats import norm, chi2 as CH
from vlib import *
out = {}
F = load("fd_d23"); F35 = load("fd_nd_d150")
D = {r: load(f"data_{r}_d23") for r in ("427929", "428770")}
D35 = {r: load(f"data_{r}_d150") for r in ("427929", "428770")}
S3 = json.load(open(os.path.join(BASE, "out", "s03.json")))
HEAD = {("427929", 0): 0.1, ("427929", 1): 0.2, ("428770", 0): 0.1, ("428770", 1): 0.1}

print("== C. eps_B(m) = N_po_cieciu(m)/N_wszystkie(m), dane 2-3 GeV, biny 5 MeV (H: wszystkie kandydaty klasy); fit exp(pol2) do pasm |m-m0|>80 MeV w 1.70-2.05")
ME, LE = D["427929"]["ME"], D["427929"]["LE"]; ctr = 0.5 * (ME[1:] + ME[:-1])
for (run, j), wp in HEAD.items():
    d = D[run]; th = mc_wp(F, j, wp)
    allm = d["H"][0, j].sum(1)
    npass = np.histogram(d["m"][d["lg"][:, j] >= th], ME)[0].astype(float)
    r = npass / allm; er = np.sqrt(np.maximum(npass, 1)) / allm
    x = (ctr - M0) / 0.1
    sb = np.abs(ctr - M0) > 0.08
    def reg(mask, lab):
        return f"{lab} {1e4*npass[mask].sum()/allm[mask].sum():.3f}"
    cf, cov = np.polyfit(x[sb], np.log(r[sb]), 2, w=(r / er)[sb], cov="unscaled")
    pred = np.exp(np.polyval(cf, x)) * allm
    wn = np.abs(ctr - M0) < WIN
    exc = npass[wn].sum() - pred[wn].sum()
    # blad interpolacji: propagacja kowariancji
    J = np.array([(pred[wn] * x[wn] ** 2).sum(), (pred[wn] * x[wn]).sum(), pred[wn].sum()]); epred = math.sqrt(J @ cov @ J)
    q = templ(F, j, th); o23 = json.load(open(os.path.join(BASE, "out", "s02.json")))[f"shape|2-3|{NM[j]}|{wp}"]
    exp_exc = S3[f"var|{run}|{NM[j]}|{wp}"]["nom exp2 1.75-2.00"]["N"] * (q["fwin"] + o23["refl_win"])
    c2 = float((((npass - pred) / np.sqrt(np.maximum(pred, 1)))[sb] ** 2).sum())
    out[f"epsB|{run}|{NM[j]}"] = dict(slope=cf[1], slope_err=math.sqrt(cov[1, 1]), curv=cf[0], curv_err=math.sqrt(cov[0, 0]), excess=exc, excess_err=math.sqrt(npass[wn].sum() + epred ** 2), pred=float(pred[wn].sum()), exp_from_fit=exp_exc, chi2_sb=c2, nsb=int(sb.sum()) - 3)
    print(f"  {run} {NM[j]} WP{int(wp*100)}: eps_B x1e4: " + ", ".join(reg(m_, l_) for l_, m_ in (("1.70-1.75", ctr < 1.75), ("1.75-1.80", (ctr > 1.75) & (ctr < 1.80)), ("okno±55", wn), ("1.93-1.98", (ctr > 1.93) & (ctr < 1.98)), ("2.00-2.05", ctr > 2.0))))
    print(f"       ln eps_B: nachylenie {cf[1]:+.3f} ± {math.sqrt(cov[1,1]):.3f} /100 MeV, krzywizna {cf[0]:+.3f} ± {math.sqrt(cov[0,0]):.3f} /(100 MeV)^2, chi2 pasm {c2:.0f}/{sb.sum()-3};"
          f" okno: obserw. {npass[wn].sum():.0f}, interpolacja {pred[wn].sum():.0f} ± {epred:.0f} -> nadwyzka {exc:.0f} ± {math.sqrt(npass[wn].sum() + epred**2):.0f} (z fitu nominalnego: N*(fwin+odbicie w oknie) = {exp_exc:.0f})")

print("== 3-5 GeV, dane: y = N/eps drzewa i sieci przy WP50..10 (kontrola przenoszenia eps przy ciasnym cieciu tam, gdzie fit etapu 4 dziala od WP70)")
y35 = {}
for run, d in D35.items():
    for j in (0, 1):
        txt = []
        for wp in (0.5, 0.3, 0.2, 0.1):
            th = mc_wp(F35, j, wp)
            if th < (1.0, 1.4)[j]:
                txt.append(f"WP{int(wp*100)}: prog {th:.2f} < ciecia redukcji"); continue
            q = templ(F35, j, th); s = d["lg"][:, j] >= th
            H, e = hist(d["m"][s]); f = fit2(H, e, q); fs = fit2(H, e, q, free_shape=True)
            y35[(run, j, wp)] = (f["N"] / q["eps"], f["err"] / q["eps"])
            out[f"y35|{run}|{NM[j]}|{wp}"] = dict(N=f["N"], err=f["err"], eps=q["eps"], y=f["N"] / q["eps"], B=f["B"], dmu=fs["a"], sc=fs["b"])
            txt.append(f"WP{int(wp*100)}: N {f['N']:.0f} ± {f['err']:.0f}, y {f['N']/q['eps']:.0f} ± {f['err']/q['eps']:.0f}, S/B {f['N']*q['fwin']/f['B']:.2f}")
        print(f"  {run} {NM[j]}: " + "; ".join(txt))

print("== test 'pik ∝ tlo' (rzezbienie skaluje sie z tlem, sygnal nie): stosunek runow 427929/428770")
def sm(run, j, k="3 ciasne (eps<0.3)"):
    o = S3[f"slicefit|{run}|{NM[j]}|{k}"]; return o["mean"], o["err"]
rDp, eDp = 3441 / 3055, (3441 / 3055) * math.hypot(260 / 3441, 363 / 3055)     # D+ (D23 RAPORT §0)
print(f"  D+ (D23, H WP20): y427/y428 = {rDp:.2f} ± {eDp:.2f}")
rng = np.random.default_rng(5)
for j in (0, 1):
    (a, ea), (b, eb) = sm("427929", j), sm("428770", j)
    r = a / b; er = r * math.hypot(ea / a, eb / b)
    wp = 0.1
    Na, Nb_ = S3[f"var|427929|{NM[j]}|{HEAD[('427929', j)]}"]["nom exp2 1.75-2.00"], S3[f"var|428770|{NM[j]}|{HEAD[('428770', j)]}"]["nom exp2 1.75-2.00"]
    # przy tym samym WP10: N i B obu runow z wycinka 0.1-0.0
    sa, sb_ = S3[f"slice|427929|{NM[j]}|0.1-0.0"], S3[f"slice|428770|{NM[j]}|0.1-0.0"]
    rB = sa["B"] / sb_["B"]; rN = sa["N"] / sb_["N"]; erN = rN * math.hypot(sa["err"] / sa["N"], sb_["err"] / sb_["N"])
    # N_run = S_run + f*B_run, S427/S428 = rDp  ->  f = (N427 - rDp N428)/(B427 - rDp B428)
    smp = []
    for _ in range(20000):
        n1, n2, rr = rng.normal(sa["N"], sa["err"]), rng.normal(sb_["N"], sb_["err"]), rng.normal(rDp, eDp)
        smp.append((n1 - rr * n2) / (sa["B"] - rr * sb_["B"]))
    f0 = (sa["N"] - rDp * sb_["N"]) / (sa["B"] - rDp * sb_["B"]); lo, hi = np.quantile(smp, [0.16, 0.84])
    out[f"runratio|{NM[j]}"] = dict(r_y=r, r_y_err=er, rN10=rN, rN10_err=erN, rB10=rB, f=f0, f_lo=lo, f_hi=hi, fakeN427=f0 * sa["B"], N427=sa["N"])
    print(f"  D0 {NM[j]}: <y>(eps<0.3) 427929/428770 = {r:.2f} ± {er:.2f}; przy WP10: N427/N428 = {rN:.2f} ± {erN:.2f}, B427/B428 = {rB:.2f} -> hipoteza 'pik ∝ tlo' odrzucona na {(rN-rB)/erN:.1f}σ; "
          f"skladowa ∝ tlo f = {f0:+.3f} [{lo:+.3f}, {hi:+.3f}] (68 %) × B_okno -> N_falszywe(427929) = {f0*sa['B']:+.0f} [{lo*sa['B']:+.0f}, {hi*sa['B']:+.0f}] z N = {sa['N']:.0f}")

print("== test 5: wiarygodnosc uzysku")
nD0, nDp = 51920, 27709        # sygnal MC fd_nd 2-3 GeV po preselekcji (te same 15 plikow): D0 poprawna hipoteza (ten plik), D+ (dnet4_d23.json n_sig_mc)
BR0, BRp = 0.0395, 0.0938
for run, yDp, eDp_ in (("427929", 3441, 260), ("428770", 3055, 363)):
    for j in (0, 1):
        y, e = sm(run, j)
        rat = y / yDp; erat = rat * math.hypot(e / y, eDp_ / yDp)
        # skorygowany stosunek D+/D0 przy zalozeniu, ze fd_nd wymusza rozpad kazdego D (gen D0/D+ = 1/0.521 Pythia A3): eff0/eff+ = (nD0/nDp)*0.521
        effr = nD0 / nDp * 0.521
        corr = (yDp / y) * (BR0 / BRp) * effr
        out[f"plaus|{run}|{NM[j]}"] = dict(y=y, yerr=e, ratio=rat, ratio_err=erat, DpD0=corr, DpD0_err=corr * erat / rat, effr=effr)
        print(f"  {run} {NM[j]}: y(D0) = {y:.0f} ± {e:.0f}; y(D0)/y(D+) = {rat:.2f} ± {erat:.2f}; oczekiwane dla D+/D0 = 0.40-0.50: {(1/0.50)*(BR0/BRp)*effr:.2f}-{(1/0.40)*(BR0/BRp)*effr:.2f} (eff0/eff+ = {effr:.2f}, zalozenie);"
              f" odwrotnie: (D+/D0)_skoryg = {corr:.2f} ± {corr*erat/rat:.2f} (stat)")
        for wp in (0.3, 0.2):
            if (run, j, wp) in y35:
                y3, e3 = y35[(run, j, wp)]
                print(f"       y(2-3)/y(3-5, WP{int(wp*100)}) = {y/y3:.2f} ± {(y/y3)*math.hypot(e/y, e3/y3):.2f}   [D+: 0.79 ± 0.08 (427929), 0.52 ± 0.07 (428770)]")
print(f"  MC fd_nd: sygnal D0 po preselekcji 2-3: {nD0} (15 plikow), 3-5: {(F35['cat']==0).sum()} (24 pliki) — rozne produkcje, stosunku widma z MC nie da sie znormalizowac")

print("== profil -2lnL po sc (427929 drzewo WP10, dmu wolne) — czy szerokosc jest okreslona")
d = D["427929"]; th = mc_wp(F, 0, 0.1); q = templ(F, 0, th); H, e = hist(d["m"][d["lg"][:, 0] >= th])
from iminuit import Minuit
prof = []
for sc in np.arange(0.6, 1.61, 0.1):
    best = min(fit_fixed["fval"] for fit_fixed in [fit2(H, e, dict(q, sigT=type("T", (), {"bins": staticmethod(lambda ed, m0, a, b, _s=sc, _d=dm: q["sigT"].bins(ed, m0, _d, _s))})()), nbkg=2) for dm in (-0.004, -0.002, 0.0, 0.002)])
    prof.append((sc, best))
f0 = min(p[1] for p in prof)
print("  " + "; ".join(f"sc {s:.1f}: Δ {v-f0:.1f}" for s, v in prof))
out["sc_profile"] = [(float(s), float(v - f0)) for s, v in prof]

print("== liczba prob (D0 2-3): 2 runy × 2 klasyfikatory × (WP wybrany z pasm + sasiedzi) = 10 fitow obejrzanych; pozycja piku ustalona")
for (run, j), wp in HEAD.items():
    v = S3[f"var|{run}|{NM[j]}|{wp}"]
    for lab in ("nom exp2 1.75-2.00", "exp4", "exp2 1.78-1.96", "exp4 1.70-2.05"):
        z = v[lab]["ZLR"]; p = norm.sf(z); pg = 1 - (1 - p) ** 10
        out[f"trials|{run}|{NM[j]}|{lab}"] = dict(ZLR=z, Zglob=float(norm.isf(pg)))
        print(f"  {run} {NM[j]} {lab:20s}: Z_LR {z:.2f} -> po 10 probach {norm.isf(pg):.2f}")
json.dump(out, open(os.path.join(BASE, "out", "s04.json"), "w"), indent=1, default=float)
