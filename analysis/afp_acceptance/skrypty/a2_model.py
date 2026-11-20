#!/usr/bin/env python3
"""a2: czynnik modelowy z SD_var/SD_pred (poziom cząstek Pythii) × akceptancja fast-sim A(pT) (res_mc.json).
Per model dyfrakcji (proces SD, D0, pT 3–30 GeV): f_win = udział D z protonem w oknie 0,03–0,10 (→ ekstrapolacja do SD całkowitego = 1/f_win),
udział |t| < 1 GeV², udział rdzenia 0,04–0,09 w oknie, A_fold = Σ dN/dt × A(t) (akceptancja geometryczna zwinięta z rozkładem t modelu).
Zapis: tabele/model_factor.md, tabele/model_factor.json."""
import json, glob, os, numpy as np
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); CW = os.path.join(W, '..')
mc = json.load(open(W + '/runs/lxplus/res_mc.json'))['mc_mb_sd']
XIB = [0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.095, 0.10, 0.105, 0.11, 0.12]; PTB = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0])
def A_pt(lo, hi):
    k = np.zeros(7); n = np.zeros(7)
    for s in (0, 1):
        for i, row in enumerate(mc[f'A_xi_pt_s{s}']):
            if XIB[i] >= lo - 1e-9 and XIB[i + 1] <= hi + 1e-9:
                for j, e in enumerate(row):
                    if e[0] is not None: k[j] += e[0] * e[2]; n[j] += e[2]
    return k / np.maximum(n, 1)
A_win = A_pt(0.03, 0.10); A_core = A_pt(0.04, 0.09)
tc = np.arange(41) * 0.05 + 0.025; tc[-1] = 2.5; ptc = np.sqrt(tc); jb = np.clip(np.searchsorted(PTB, ptc, side='right') - 1, 0, 6)
rows = []; out = {'A_pt_window': A_win.tolist(), 'A_pt_core': A_core.tolist(), 'pt_edges': PTB.tolist(), 'models': {}}
for fn in sorted(glob.glob(CW + '/SD_var_2026-11-19/runs/lxplus/merged/dv_*.json') + glob.glob(CW + '/SD_pred/runs/lxplus/merged/dif_*.json')):
    d = json.load(open(fn)); h = d['hists']; nm = os.path.basename(fn)[:-5]
    if 'SD|D0|pt_tag' not in h or 'SD|D0|t_tag' not in h: continue
    g = lambda k: np.array(h[k]['s'], float)
    inc = g('SD|D0|pt_incl')[1:].sum(); tag = g('SD|D0|pt_tag')[1:].sum(); t1 = g('SD|D0|pt_tag_t1')[1:].sum(); tt = g('SD|D0|t_tag'); xw = g('SD|D0|xiw')
    if tag < 200: continue
    r = dict(n_tag=tag, f_win=tag / inc, f_win_err=np.sqrt(tag) / inc, f_t1=t1 / tag, core_frac=xw[4:9].sum() / xw[3:10].sum(), A_fold_win=float((tt * A_win[jb]).sum() / tt.sum()),
             A_fold_core=float((tt * A_core[jb]).sum() / tt.sum()), b_eff=float(1 / max((tt * tc).sum() / tt.sum(), 1e-9)), xi_slope=float(xw[3:6].sum() / xw[7:10].sum()),
             dd_share=float(g('DD|D0|pt_tag')[1:].sum() * 1.0) if 'DD|D0|pt_tag' in h else None)
    out['models'][nm] = r
    rows.append(f"| {nm} | {tag:.0f} | {r['f_win'] * 100:.2f} ± {r['f_win_err'] * 100:.2f} | {1 / r['f_win']:.2f} | {r['f_t1']:.3f} | {r['core_frac']:.3f} | {r['b_eff']:.2f} | {r['A_fold_win']:.4f} | {r['A_fold_core']:.4f} | {r['xi_slope']:.2f} |")
M = out['models']; v = lambda k: np.array([m[k] for m in M.values()])
env = {k: [float(v(k).min()), float(M['dv_A3off'][k]), float(v(k).max())] for k in ('f_win', 'f_t1', 'core_frac', 'A_fold_win', 'A_fold_core', 'b_eff')}
out['envelope_min_nom_max'] = env
L = ['# Czynnik modelowy: SD z D⁰ (pT 3–30 GeV, |η| < 2,2), Pythia 8.317, poziom cząstek × A(pT) z fast-sim AFP\n',
     'f_win = N(D, proton w 0,03–0,10 po dowolnej stronie)/N(D) w procesie SD; 1/f_win = ekstrapolacja „widzialne → SD całkowite”; f_t1 = udział |t| < 1 GeV² wśród tagowanych;',
     'rdzeń = udział ξ 0,04–0,09 w ξ 0,03–0,10 (D pT ≥ 2); b_eff = 1/⟨|t|⟩; A_fold = akceptancja geometryczna fast-sim zwinięta z dN/dt modelu (okno / rdzeń); ξ-nachylenie = N(0,03–0,06)/N(0,07–0,10).\n',
     f"A(pT) fast-sim w oknie 0,03–0,10 (biny pT {PTB.tolist()}): {np.round(A_win, 4).tolist()}; w rdzeniu 0,04–0,09: {np.round(A_core, 4).tolist()}\n",
     '| model | N_tag | f_win [%] | 1/f_win | f_t1 | rdzeń | b_eff [GeV⁻²] | A_fold okno | A_fold rdzeń | ξ-nachylenie |\n|---|---|---|---|---|---|---|---|---|---|'] + rows
L.append('\n**Obwiednia (min / nominał A3off / max) po wszystkich modelach z tabeli:**\n')
for k, e in env.items(): L.append(f"- {k}: {e[0]:.4f} / {e[1]:.4f} / {e[2]:.4f}  (× {e[0] / e[1]:.3f} … × {e[2] / e[1]:.3f})")
open(W + '/tabele/model_factor.md', 'w').write('\n'.join(L) + '\n'); json.dump(out, open(W + '/tabele/model_factor.json', 'w'), indent=1)
print('\n'.join(L))
