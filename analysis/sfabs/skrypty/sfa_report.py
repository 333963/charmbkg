#!/usr/bin/env python3
"""SFABS A: odczyt drabiny -> tabele/a_tab.md, out/a_summary.json. Wzory: KRYTERIA par. A."""
import json, sys, math, numpy as np
sys.path.insert(0, 'skrypty')
from sfa_ladder import kappa, sf70_k
D = '/home/k/Documents/charm/charm_work/SFABS_wynik_2026-11-20/'
MC = {(r['sp'], r['bin'], r['L']): r for r in json.load(open(D + 'out/a_mc.json'))}
DA = {(r['sp'], r['bin'], r['L'], r['run']): r for r in json.load(open(D + 'out/a_data.json'))}
try: TO = {(r['sp'], r['bin'], r['run']): r for r in json.load(open(D + 'out/a_toys.json'))}
except Exception: TO = {}
BW = {0.9: (0.90, 1.05), 0.95: (0.95, 1.025), 0.8: (0.80, 1.10)}
out = {}; L1 = []; L2 = []; L3 = []
L1.append('| gat. | bin | L | b_L (mb_nd) | wstrz. skok: SF70 odz. / prawda | wstrz. κ=1,3: odz. / prawda | σ_testu | G_A1 |\n|---|---|---|---|---|---|---|---|')
L2.append('| gat. | bin | L | run | r_L = f_dane/f_MC | r_L^kor | κ̂ | SF(L) | SF70 | ± stat | ± σ(b) | ± szczebel | ± razem | p(χ²) | G_A2 (pull śr./rms; obc.) | bramki |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for sp in ('Dplus', 'D0'):
    for bn in ('3-5', '5-8', '8-30'):
        for L in (0.9, 0.95, 0.8):
            m = MC[(sp, bn, L)]; e7, eL = m['e7'], m['eL']; z = m['modes']['zero']
            rep0 = np.array(z['f_rep'], float)
            b = z['f'] / z['f_true']; sb = max(np.nanstd(rep0, ddof=1), 1.3 * z['f_err']) / z['f_true']
            g1 = True; cells = []; st = 0
            for md in ('step', 'kap'):
                o = m['modes'][md]
                rr = o['f'] / z['f']; rt = o['f_true'] / z['f_true']
                d = np.array(o['f_rep'], float) / rep0; sd = float(np.nanstd(d, ddof=1))
                rec, _ = sf70_k(rr, e7, eL)
                srec = abs(sf70_k(rr + sd, e7, eL)[0] - sf70_k(rr - sd, e7, eL)[0]) / 2
                g1 &= (abs(rec - o['SF70_true']) <= 2 * srec) and (abs(rec - o['SF70_true']) <= 0.05)
                cells.append(f"{rec:.3f} / {o['SF70_true']:.3f}"); st = max(st, srec)
            L1.append(f"| {sp} | {bn} | {L:.2f} | {b:.3f} ± {sb:.3f} | {cells[0]} | {cells[1]} | {st:.3f} | {'TAK' if g1 else 'NIE'} |")
            for run in ('427929', '428770'):
                d = DA[(sp, bn, L, run)]
                r = d['f'] / d['f_mc']
                rep = np.array(d['f_rep'], float)
                sstat_f = max(np.nanstd(rep, ddof=1), 1.3 * (d['f_err'] or 0)) / d['f_mc']
                rk = r / b
                sf, k = sf70_k(rk, e7, eL)
                dsf = lambda x: sf70_k(x, e7, eL)[0]
                amp = (dsf(rk + 0.01) - dsf(rk - 0.01)) / 0.02
                s_stat = abs(amp) * sstat_f / b; s_b = abs(amp) * rk * sb / b
                sfL = (1 - k * (1 - eL)) / eL
                s_sz = max(abs(sf - rk) / math.sqrt(3), (BW[L][1] - BW[L][0]) / math.sqrt(12) * rk)
                tot = math.sqrt(s_stat ** 2 + s_b ** 2 + s_sz ** 2)
                kraw = kappa(rk, e7, eL)
                g3 = bool(d['valid'] and not d['at_limit'] and d['p'] >= 0.01 and 0.5 < kraw < 2.0)
                t = TO.get((sp, bn, run)) if L == 0.9 else None
                g2 = None; g2s = '—'; infl = 1.0
                if t:
                    g2 = True; ss = []
                    for lab in ('k1.0', 'k1.3'):
                        q = t[lab]
                        g2 &= abs(q['pull_mean']) <= 0.35 and q['pull_rms'] <= 1.6 and abs(q['bias_rel']) <= 0.03
                        infl = max(infl, q['pull_rms'] / 1.3 if q['pull_rms'] > 1.3 else 1.0)
                        ss.append(f"{q['pull_mean']:+.2f}/{q['pull_rms']:.2f}; {100 * q['bias_rel']:+.1f} %")
                    g2s = ' ∥ '.join(ss)
                ok = g1 and g3 and (g2 is not False)
                gate = ('TAK' if ok else 'NIE') + (' (pomiar)' if ok and tot <= 0.15 else (' (ograniczenie: σ > 0,15)' if ok else f" [G1 {'+' if g1 else '−'} G2 {'+' if g2 else ('−' if g2 is False else '?')} G3 {'+' if g3 else '−'}]"))
                L2.append(f"| {sp} | {bn} | {L:.2f} | {run} | {r:.3f} ± {sstat_f:.3f} | {rk:.3f} | {kraw:.2f} | {sfL:.3f} | **{sf:.3f}** | {s_stat:.3f} | {s_b:.3f} | {s_sz:.3f} | {tot:.3f} | {d['p']:.3f} | {g2s} | {gate} |")
                out[f"{sp}|{bn}|{L}|{run}"] = dict(r=r, r_stat=sstat_f, b=b, sb=sb, rk=rk, kappa=kraw, SFL=sfL, SF70=sf, stat=s_stat, sbias=s_b, sszcz=s_sz, tot=tot, amp=amp, p=d['p'], G1=g1, G2=g2, G3=g3, ok=ok,
                                                   NA=d['NA'], NB=d['NB'], e7=e7, eL=eL)
# srednia runow (L = 0,90): b wspolne
L3.append('| gat. | bin | SF70 427929 | SF70 428770 | χ² zgodności (stat) / p | średnia ważona r^kor | SF70 śr. | ± stat | ± σ(b) | ± szczebel | ± razem | kontrola L=0,95 | kontrola L=0,80 | status |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
from scipy.stats import chi2 as C2
for sp in ('Dplus', 'D0'):
    for bn in ('3-5', '5-8', '8-30'):
        res = {}
        for L in (0.9, 0.95, 0.8):
            a, c = out[f"{sp}|{bn}|{L}|427929"], out[f"{sp}|{bn}|{L}|428770"]
            w1, w2 = 1 / a['r_stat'] ** 2, 1 / c['r_stat'] ** 2
            rm = (a['r'] * w1 + c['r'] * w2) / (w1 + w2); srm = 1 / math.sqrt(w1 + w2)
            chi = (a['r'] - c['r']) ** 2 / (a['r_stat'] ** 2 + c['r_stat'] ** 2)
            rk = rm / a['b']; sf, k = sf70_k(rk, a['e7'], a['eL']); amp = (sf70_k(rk + 0.01, a['e7'], a['eL'])[0] - sf70_k(rk - 0.01, a['e7'], a['eL'])[0]) / 0.02
            s_stat = abs(amp) * srm / a['b']; s_b = abs(amp) * rk * a['sb'] / a['b']
            s_sz = max(abs(sf - rk) / math.sqrt(3), (BW[L][1] - BW[L][0]) / math.sqrt(12) * rk)
            res[L] = dict(SF70=sf, stat=s_stat, sbias=s_b, sszcz=s_sz, tot=math.sqrt(s_stat ** 2 + s_b ** 2 + s_sz ** 2), chi=chi, p=float(C2.sf(chi, 1)), rk=rk, ok=a['ok'] and c['ok'], ok427=a['ok'], kappa=kappa(rk, a['e7'], a['eL']))
        n = res[0.9]; a, c = out[f"{sp}|{bn}|0.9|427929"], out[f"{sp}|{bn}|0.9|428770"]
        status = ('pomiar' if n['tot'] <= 0.15 else 'ograniczenie (σ > 0,15)') if n['ok'] else ('tylko 427929' if n['ok427'] else 'NIE przechodzi bramek')
        if n['p'] < 0.05: status += '; runy niezgodne (p < 0,05) → bez średniej'
        L3.append(f"| {sp} | {bn} | {a['SF70']:.3f} ± {a['tot']:.3f} | {c['SF70']:.3f} ± {c['tot']:.3f} | {n['chi']:.2f} / {n['p']:.2f} | {n['rk']:.3f} | **{n['SF70']:.3f}** | {n['stat']:.3f} | {n['sbias']:.3f} | {n['sszcz']:.3f} | **{n['tot']:.3f}** | "
                  f"{res[0.95]['SF70']:.3f} ± {res[0.95]['tot']:.3f} | {res[0.8]['SF70']:.3f} ± {res[0.8]['tot']:.3f} | {status} |")
        out[f"AVG|{sp}|{bn}"] = {str(k): v for k, v in res.items()} | {'status': status}
open(D + 'tabele/a_tab.md', 'w').write('### A1. Domknięcie na mb_nd (pełne tło): poprawka b_L i wstrzyknięcia (G_A1)\n\n' + '\n'.join(L1) + '\n\n### A2. Dane per run i szczebel\n\n' + '\n'.join(L2) + '\n\n### A3. Średnia runów (b wspólne), L = 0,90 nominał\n\n' + '\n'.join(L3) + '\n')
json.dump(out, open(D + 'out/a_summary.json', 'w'), indent=1)
print('\n'.join(L2)); print(); print('\n'.join(L3))
