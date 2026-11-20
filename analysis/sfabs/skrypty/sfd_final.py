#!/usr/bin/env python3
"""SFABS: B (tabela + bramki), C (SF_ref + przeniesienie), D (polaczenie, przekroj — wlasna kopia rachunku; XSEC tylko odczyt). -> tabele/{b_tab,c_tab,d_tab}.md, out/final.json"""
import json, math, numpy as np
from scipy.stats import chi2 as C2
D = '/home/k/Documents/charm/charm_work/SFABS_wynik_2026-11-20/'
W = '/home/k/Documents/charm/charm_work/'
bm, bd, bt = (json.load(open(D + f'out/b_{x}.json')) for x in ('mc', 'data', 'toys'))
gof = json.load(open(D + 'out/b_gof_fix.json'))
A = json.load(open(D + 'out/a_summary.json'))
sens = json.load(open(D + 'out/c_sens.json'))
dn = json.load(open(W + 'DNET/runs/lxplus4/dnet4_results.json'))['species']['D0']['sf']
xs = json.load(open(W + 'XSEC_2026-11-19/out/xsec_results.json'))
BINS = ['3-5', '5-8', '8-30']; RUNS = ['427929', '428770']; CLS = {'427929': 'MB427929', '428770': 'all'}
F = {}
# ---------------- B
L = ['### B1. Domknięcie na mb_nd (pełne tło, jedna realizacja; model nominalny z ε_R = ε) — G_B1\n', '| bin | ε prawda (N_S, N_R) | ε fit ± HESSE | rozrzut Poissona | pull | N fit / prawda | p(χ²) | wariant swobodne R | G_B1 |', '|---|---|---|---|---|---|---|---|---|']
pulls = []
for bn in BINS:
    r = bm[bn]; f = r['fit']; t = r['truth']; pl = (f['eps'] - t['eps']) / max(f['eps_err'], r['boot_rms']); pulls.append(pl)
    p = f['p'] if f['p'] == f['p'] and f['p'] is not None else gof.get(f'mb|{bn}', {}).get('p')
    L.append(f"| {bn} | {t['eps']:.3f} ({t['N']}, {t['NR']}) | {f['eps']:.3f} ± {f['eps_err']:.3f} | {r['boot_rms']:.3f} | {pl:+.2f} | {f['N']:.0f} ± {f['N_err']:.0f} / {t['N']} | {p:.2f} | {r['var']['freeR']['eps']:.3f} | {'TAK' if abs(pl) <= 2 else 'NIE'} |")
L.append(f"\nŚredni pull po binach: {np.mean(pulls):+.2f} (próg 1,0) → G_B1 {'spełnione' if abs(np.mean(pulls)) <= 1 and all(abs(x) <= 2 for x in pulls) else 'NIE'}.\n")
L += ['### B2. Zabawki przy czystości danych (sygnał fd_nd w tle z pasm masy danych + składowa R z mb_nd; 30 na punkt) — G_B2\n', '| run | bin | wstrz. ε = ε_MC: obc. ± błąd, pull śr./rms | wstrz. 0,85 ε_MC: obc., pull śr./rms | Δ(odzyskane) − Δ(prawda) | 0,5 σ_stat / 0,05 ε | G_B2 | mnożnik błędu |', '|---|---|---|---|---|---|---|---|']
B = {}
for run in RUNS:
    for bn in BINS:
        k = f'{run}|{bn}'; t = bt[k]; d = bd[k]; f = d['fit']; a, c = t['eps_MC'], t['0.85eps_MC']
        lim1 = 0.5 * f['eps_err']; lim2 = 0.05 * d['eps_mc']
        dd = (a['mean'] - c['mean']) - (a['eps_true'] - c['eps_true']); sdd = math.hypot(a['bias_err'], c['bias_err'])
        ok = all(abs(q['bias']) <= lim1 and abs(q['bias']) <= lim2 and q['pull_rms'] <= 1.6 for q in (a, c)) and abs(dd) <= 2 * sdd
        infl = max([q['pull_rms'] for q in (a, c) if 1.3 < q['pull_rms'] <= 1.6] + [1.0])
        p = f['p'] if (f['p'] is not None and f['p'] == f['p']) else gof[k]['p']
        g3 = bool(f['valid'] and not f['at_limit'] and p >= 0.01)
        stat = max(d['SF_hesse'] * infl, d['SF_boot'])
        tot = math.hypot(stat, d['SF_model'])
        B[k] = dict(SF=d['SF'], stat=stat, model=d['SF_model'], tot=tot, G2=ok, G3=g3, p=p, ok=ok and g3, N=f['N'], N_err=f['N_err'], eps=f['eps'], hesse=d['SF_hesse'], boot=d['SF_boot'], freeR=d['var']['freeR']['eps'] / d['eps_mc'])
        L.append(f"| {run} | {bn} | {a['bias']:+.3f} ± {a['bias_err']:.3f}, {a['pull_mean']:+.2f}/{a['pull_rms']:.2f} | {c['bias']:+.3f} ± {c['bias_err']:.3f}, {c['pull_mean']:+.2f}/{c['pull_rms']:.2f} | {dd:+.3f} ± {sdd:.3f} | {lim1:.3f} / {lim2:.3f} | {'TAK' if ok else 'NIE'} | {infl:.2f} |")
L += ['\n### B3. Dane: SF(WP70) z fitu 2D; porównanie z dnet4 (fit kształtu Δm)\n', '| run | bin | N_D* | ε_dane | SF | ± stat (HESSE / bootstrap) | ± model | razem | p(χ²) | G_B2 / G_B3 | status | wariant swobodne R | dnet4 SF(WP70) |', '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for run in RUNS:
    for bn in BINS:
        k = f'{run}|{bn}'; b = B[k]; q = dn[f"D0.xgb.H.sb.n.s1|{run}|{CLS[run]}|{bn}|0.7"]
        st = ('pomiar' if b['tot'] <= 0.25 else 'ograniczenie') if b['ok'] else 'NIE przechodzi'
        L.append(f"| {run} | {bn} | {b['N']:.0f} ± {b['N_err']:.0f} | {b['eps']:.3f} | **{b['SF']:.3f}** | {b['stat']:.3f} ({b['hesse']:.3f} / {b['boot']:.3f}) | {b['model']:.3f} | {b['tot']:.3f} | {b['p']:.3f} | {'+' if b['G2'] else '−'} / {'+' if b['G3'] else '−'} | {st} | {b['freeR']:.3f} | {q['SF']:.3f} ± {q['SF_err']:.3f}{'' if q['determined'] else ' (nieokr.)'} |")
open(D + 'tabele/b_tab.md', 'w').write('\n'.join(L) + '\n')
F['B'] = B
# ---------------- C
ref = []
for run in RUNS:
    for bn in ('5-8', '8-30'):
        k = f'{run}|{bn}'; b = B[k]
        if b['ok']:
            ref.append(dict(k=k, bn=bn, SF=b['SF'], stat=b['stat'], model=b['model'], src='fit 2D'))
        else:
            q = dn[f"D0.xgb.H.sb.n.s1|{run}|{CLS[run]}|{bn}|0.7"]; ref.append(dict(k=k, bn=bn, SF=q['SF'], stat=q['SF_err'], model=0.0, src='dnet4 (fit 2D nie przeszedł G_B2)'))
w = np.array([1 / r['stat'] ** 2 for r in ref]); v = np.array([r['SF'] for r in ref])
sfref = float((w * v).sum() / w.sum()); sref = float(1 / math.sqrt(w.sum())); chi = float((w * (v - sfref) ** 2).sum()); pchi = float(C2.sf(chi, len(ref) - 1))
mref = float((w * np.array([r['model'] for r in ref])).sum() / w.sum())
if pchi < 0.05: sref *= math.sqrt(chi / (len(ref) - 1))
wb = {bn: float(sum(wi for wi, r in zip(w, ref) if r['bn'] == bn) / w.sum()) for bn in ('5-8', '8-30')}
# wariant: tylko fity 2D, ktore przeszly
r2 = [r for r in ref if r['src'] == 'fit 2D']; w2 = np.array([1 / r['stat'] ** 2 for r in r2]); alt = float((w2 * np.array([r['SF'] for r in r2])).sum() / w2.sum()); salt = float(1 / math.sqrt(w2.sum()))
# wariant: same dnet4
rd = [dn[f"D0.xgb.H.sb.n.s1|{run}|{CLS[run]}|{bn}|0.7"] for run in RUNS for bn in ('5-8', '8-30')]; wd = np.array([1 / q['SF_err'] ** 2 for q in rd]); dalt = float((wd * np.array([q['SF'] for q in rd])).sum() / wd.sum()); sdalt = float(1 / math.sqrt(wd.sum()))
sv = lambda sp, bn, var, i: float(np.mean([x[i] for x in sens[f'{sp}|{bn}|{var}']]))
seedn = lambda sp, bn, var: abs(sens[f'{sp}|{bn}|{var}'][0][0] - sens[f'{sp}|{bn}|{var}'][1][0])
RES = ['T1.05', 'T1.10', 'F1.05', 'F1.10']; DEDX = ['dEdx0.95', 'dEdx1.05', 'dEdxS5']
L = ['### C1. SF_ref: tag D*, D⁰ ≥ 5 GeV, WP70\n', '| komórka | SF | ± stat | ± model | źródło |', '|---|---|---|---|---|']
for r in ref: L.append(f"| {r['k']} | {r['SF']:.3f} | {r['stat']:.3f} | {r['model']:.3f} | {r['src']} |")
L.append(f"\n**SF_ref = {sfref:.3f} ± {sref:.3f} (stat) ± {mref:.3f} (model fitu, średnia ważona liniowo)**; χ²/ndf = {chi:.2f}/{len(ref) - 1}, p = {pchi:.2f}; wagi binów: 5–8: {wb['5-8']:.2f}, 8–30: {wb['8-30']:.2f}. "
         f"Warianty: tylko fity 2D po bramkach {alt:.3f} ± {salt:.3f}; same liczby dnet4 {dalt:.3f} ± {sdalt:.3f}.\n")
L += ['### C2. Czułość SF(WP70) na MC (fd_nd, sygnał z prawdą): z migracją przez próg L_xy / bez (D⁺: wchodzących nie da się odtworzyć — jedna liczba)\n', '| wariant | D⁰ 3–5 | D⁰ 5–8 | D⁰ 8–30 | D⁺ 3–5 | D⁺ 5–8 | D⁺ 8–30 |', '|---|---|---|---|---|---|---|']
for var in ['T1.02', 'F1.02'] + RES + DEDX:
    L.append(f"| {var} | " + ' | '.join(f"{sv('D0', b, var, 0):.3f} / {sv('D0', b, var, 1):.3f}" for b in BINS) + ' | ' + ' | '.join(f"{sv('Dplus', b, var, 0):.3f}" for b in BINS) + ' |')
L += ['\n### C3. Przeniesienie: δ_T = rozdzielczość ⊕ dE/dx ⊕ statystyka MC; SF_T = SF_ref × 1\n', '| cel | max \\|SF_T/SF_ref − 1\\| rozdzielczość (wariant) | dE/dx (wariant) | szum MC | δ_T | SF_T(s = 1,05) − 1 (T / F) | **SF (C)** ± stat ± model ± przeniesienie | razem |', '|---|---|---|---|---|---|---|---|']
Cc = {}
for sp, bn in (('Dplus', '3-5'), ('Dplus', '5-8'), ('Dplus', '8-30'), ('D0', '3-5')):
    best = (0, ''); bde = (0, '')
    for var in RES + DEDX:
        for i in (0, 1):
            rr = sum(wb[b] * sv('D0', b, var, i) for b in wb)
            x = abs(sv(sp, bn, var, i) / rr - 1)
            if var in RES and x > best[0]: best = (x, f"{var}, {'z migr.' if i == 0 else 'bez migr.'}")
            if var in DEDX and x > bde[0]: bde = (x, var)
    sn = max(seedn(sp, bn, v_) for v_ in RES)
    dT = math.sqrt(best[0] ** 2 + bde[0] ** 2 + sn ** 2)
    tot = math.sqrt(sref ** 2 + mref ** 2 + (dT * sfref) ** 2)
    Cc[f'{sp}|{bn}'] = dict(SF=sfref, stat=sref, model=mref, transfer=dT * sfref, tot=tot)
    L.append(f"| {sp} {bn} | {best[0]:.3f} ({best[1]}) | {bde[0]:.3f} ({bde[1]}) | {sn:.3f} | {dT:.3f} | {sv(sp, bn, 'T1.05', 0) - 1:+.3f} / {sv(sp, bn, 'F1.05', 0) - 1:+.3f} | **{sfref:.3f}** ± {sref:.3f} ± {mref:.3f} ± {dT * sfref:.3f} | {tot:.3f} |")
open(D + 'tabele/c_tab.md', 'w').write('\n'.join(L) + '\n')
F['C'] = Cc; F['SF_ref'] = dict(SF=sfref, stat=sref, model=mref, chi2=chi, p=pchi, cells=ref, alt_2D=[alt, salt], alt_dnet4=[dalt, sdalt])
# ---------------- D
def a_cell(sp, bn):
    a7, a8 = A[f'{sp}|{bn}|0.9|427929'], A[f'{sp}|{bn}|0.9|428770']; av = A[f'AVG|{sp}|{bn}']['0.9']
    if a7['ok'] and a8['ok'] and av['p'] >= 0.05: return dict(SF=av['SF70'], tot=av['tot'], src='A: oba runy')
    if a7['ok'] and a8['ok']: return dict(SF=a7['SF70'], tot=max(a7['tot'], abs(a7['SF70'] - a8['SF70'])), src='A: runy niezgodne, 427929, rozrzut')
    if a7['ok']: return dict(SF=a7['SF70'], tot=a7['tot'], src='A: tylko 427929')
    if a8['ok']: return dict(SF=a8['SF70'], tot=a8['tot'], src='A: tylko 428770 (kontrola)')
    return None
L = ['### D1. Połączenie per gatunek/bin (SF drzewa H przy WP70)\n', '| gatunek | bin | A drabina | B tag D* 2D | C przeniesienie | zgodność (max odchylenie par, χ²/ndf, p) | **SF łączny** | tylko metody bezpośrednie (A, B) | waga C |', '|---|---|---|---|---|---|---|---|---|']
FIN = {}
for sp, bn in (('Dplus', '3-5'), ('Dplus', '5-8'), ('Dplus', '8-30'), ('D0', '3-5')):
    ms = {}
    a = a_cell(sp, bn)
    if a: ms['A'] = (a['SF'], a['tot'])
    if sp == 'D0':
        oks = [B[f'{r}|{bn}'] for r in RUNS if B[f'{r}|{bn}']['ok']]
        if len(oks) == 1: ms['B'] = (oks[0]['SF'], oks[0]['tot'])
        elif len(oks) == 2:
            ww = [1 / o['stat'] ** 2 for o in oks]; m_ = sum(wi * o['SF'] for wi, o in zip(ww, oks)) / sum(ww); c_ = sum(wi * (o['SF'] - m_) ** 2 for wi, o in zip(ww, oks))
            ms['B'] = (m_, math.hypot(1 / math.sqrt(sum(ww)), max(o['model'] for o in oks))) if C2.sf(c_, 1) >= 0.05 else (m_, abs(oks[0]['SF'] - oks[1]['SF']))
    c = Cc[f'{sp}|{bn}']; ms['C'] = (c['SF'], c['tot'])
    names = list(ms); val = np.array([ms[n][0] for n in names]); err = np.array([ms[n][1] for n in names]); ww = 1 / err ** 2
    mean = float((ww * val).sum() / ww.sum()); smean = float(1 / math.sqrt(ww.sum())); chi = float((ww * (val - mean) ** 2).sum()); ndf = len(names) - 1
    maxdev = max(abs(val[i] - val[j]) / math.hypot(err[i], err[j]) for i in range(len(names)) for j in range(i))
    direct = [(ms[n][0], ms[n][1]) for n in names if n != 'C']
    if direct:
        wd_ = [1 / e ** 2 for _, e in direct]; dm_ = sum(wi * v_ for wi, (v_, _) in zip(wd_, direct)) / sum(wd_); ds_ = 1 / math.sqrt(sum(wd_)); dtxt = f"{dm_:.3f} ± {ds_:.3f}"
    else: dm_ = ds_ = None; dtxt = 'brak (żadna nie przechodzi bramek)'
    if maxdev > 2:
        bestd = min(direct, key=lambda x: x[1]) if direct else (val[-1], err[-1]); mean, smean = bestd[0], float(val.max() - val.min()); note = 'NIEZGODNE > 2σ: bez uśredniania'
    else: note = ''
    FIN[f'{sp}|{bn}'] = dict(SF=mean, err=smean, methods={n: ms[n] for n in names}, chi2=chi, ndf=ndf, maxdev=maxdev, wC=float(ww[names.index('C')] / ww.sum()), direct=[dm_, ds_], A_src=(a or {}).get('src'))
    fm = lambda n: f"{ms[n][0]:.3f} ± {ms[n][1]:.3f}" if n in ms else '— (bramki)'
    L.append(f"| {sp} | {bn} | {fm('A')} | {fm('B') if sp == 'D0' else 'n.d.'} | {fm('C')} | {maxdev:.2f}σ, {chi:.2f}/{ndf}, {C2.sf(chi, ndf):.2f} | **{mean:.3f} ± {smean:.3f}** {note} | {dtxt} | {FIN[f'{sp}|{bn}']['wC']:.2f} |")
# test gatunkowy: drabina D+ po binach vs C
da = [FIN[f'Dplus|{b}']['methods'].get('A') for b in BINS]; da = [x for x in da if x]
wa = [1 / e ** 2 for _, e in da]; am = sum(wi * v_ for wi, (v_, _) in zip(wa, da)) / sum(wa); asd = 1 / math.sqrt(sum(wa)); achi = sum(wi * (v_ - am) ** 2 for wi, (v_, _) in zip(wa, da))
L.append(f"\nTest założenia C (opisowo): drabina D⁺ uśredniona po 3 binach = **{am:.3f} ± {asd:.3f}** (χ² zgodności binów {achi:.2f}/{len(da) - 1}) wobec SF_ref = {sfref:.3f} ± {sref:.3f} → różnica {(am - sfref) / math.hypot(asd, sref):+.2f}σ; "
         f"bezpośredni test ogranicza efekt specyficzny dla D⁺ (niewidoczny dla C) do ±{asd:.2f} (1σ).\n")
F['A_Dplus_allbins'] = [am, asd, achi]
# ---------------- przekroj (wlasna kopia rachunku)
L += ['### D2. Przekrój 427929 z nowym SF (własna kopia rachunku; wejście: `XSEC_2026-11-19/out/xsec_results.json` — tylko odczyt, stan z 19.11 23:49, D⁰ z modelem M8)\n',
      '| gatunek | bin | σ przy SF = 1 [µb] | SF | **σ [µb] ± stat ± syst ± L** | syst wzgl. (z SF) | SF wzgl. | razem | było (scen. 1 ± 10 %): razem | Pythia A3 |', '|---|---|---|---|---|---|---|---|---|---|']
PY = {'D0|3-5': 1230, 'Dplus|3-5': 641, 'Dplus|5-8': 254, 'Dplus|8-30': 85}
X = {}
for sp, bn in (('D0', '3-5'), ('Dplus', '3-5'), ('Dplus', '5-8'), ('Dplus', '8-30')):
    c = xs['cells'][f'{sp}|427929|{bn}']; f = FIN[f'{sp}|{bn}']
    s0 = c['sigma_nb'] / 1e3; rel = dict(c['rel']); stat = rel.pop('stat'); rel['SF'] = f['err'] / f['SF']
    syst = math.sqrt(sum(v_ ** 2 for v_ in rel.values())); s1 = s0 / f['SF']; tot = math.sqrt(stat ** 2 + syst ** 2 + c['lumi_rel'] ** 2)
    old = math.sqrt(c['tot'] ** 2 + c['lumi_rel'] ** 2)
    X[f'{sp}|{bn}'] = dict(sigma=s1, stat=s1 * stat, syst=s1 * syst, lumi=s1 * c['lumi_rel'], rel_tot=tot, sigma_SF1=s0, rel_syst=syst, rel_SF=rel['SF'], rel_stat=stat)
    L.append(f"| {sp} | {bn} | {s0:.0f} | {f['SF']:.3f} ± {f['err']:.3f} | **{s1:.0f} ± {s1 * stat:.0f} ± {s1 * syst:.0f} ± {s1 * c['lumi_rel']:.0f}** | {100 * syst:.1f} % | {100 * rel['SF']:.1f} % | {100 * tot:.1f} % | {100 * old:.1f} % | {PY[f'{sp}|{bn}']} (dane/Pythia {s1 / PY[f'{sp}|{bn}']:.2f} ± {s1 / PY[f'{sp}|{bn}'] * tot:.2f}) |")
# D+/D0 3-5
r0 = xs['ratio_DpD0']['427929|3-5']; fp, f0 = FIN['Dplus|3-5'], FIN['D0|3-5']
def comb_var(f):   # wariancja czesci C i reszty w BLUE
    wC = f['wC']; c = Cc
    return wC
vref = sref ** 2 + mref ** 2
cov = fp['wC'] * f0['wC'] * vref
rsf = fp['SF'] / f0['SF']
rel_rsf = math.sqrt(max((fp['err'] / fp['SF']) ** 2 + (f0['err'] / f0['SF']) ** 2 - 2 * cov / (fp['SF'] * f0['SF']), 0))
R = r0['R_sf_equal'] / rsf; relR = math.sqrt(r0['tot_noSF'] ** 2 + rel_rsf ** 2)
L.append(f"\n**D⁺/D⁰ (3–5 GeV, 427929): R = {r0['R_sf_equal']:.3f} (SF⁺ = SF⁰) → {R:.3f} ± {R * r0['rel']['stat']:.3f} (stat) ± {R * math.sqrt(relR ** 2 - r0['rel']['stat'] ** 2):.3f} (syst)**; SF⁺/SF⁰ = {rsf:.3f} ± {rsf * rel_rsf:.3f} "
         f"(wspólny SF_ref skraca się: kowariancja {cov:.5f}); błąd względny razem {100 * relR:.1f} % (było {100 * r0['tot_sfratio10']:.1f} % przy scenariuszu SF⁺/SF⁰ = 1 ± 10 %; bez SF {100 * r0['tot_noSF']:.1f} %). Pythia A3: {r0['pythia_A3_mbnd']:.3f}.\n")
F['D'] = FIN; F['xsec'] = X; F['ratio_3-5'] = dict(R=R, rel=relR, R_old=r0['R_sf_equal'], sf_ratio=rsf, sf_ratio_rel=rel_rsf, stat_rel=r0['rel']['stat'])
open(D + 'tabele/d_tab.md', 'w').write('\n'.join(L) + '\n')
json.dump(F, open(D + 'out/final.json', 'w'), indent=1, default=float)
for t in ('b_tab', 'c_tab', 'd_tab'): print(open(D + f'tabele/{t}.md').read())
