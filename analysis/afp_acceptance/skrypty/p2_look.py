#!/usr/bin/env python3
"""p2: pierwsze spojrzenie na tabele p1 (jednostki, krotności, akceptancja fast-sim vs ξ prawdziwe)."""
import sys, numpy as np
np.set_printoptions(linewidth=200, precision=4, suppress=True)
d = np.load(sys.argv[1]); n = list(d['names']); r = d['rows']; c = {k: r[:, i] for i, k in enumerate(n)}
print(dict(zip(d['cnt_keys'], d['cnt_vals'])), r.shape)
for k in ('extra', 'npu', 'n_prot', 'n_near', 'n_far', 'n_win') + (('nt', 'nt_win') if 'nt' in c else ()):
    u, cn = np.unique(c[k], return_counts=True); print(k, dict(zip(u[:12].tolist(), cn[:12].tolist())))
q = c['n_prot'] == 1
for k in ('xi1', 'px', 'py', 'chi2', 'x_near', 'y_near', 'x_far', 'y_far', 'ncl_near', 'ncl_far', 'xical_same'):
    print(k, np.percentile(c[k][q], [0.1, 1, 16, 50, 84, 99, 99.9]))
if 'nt' in c:
    m = c['nt'] >= 1
    for k in ('txi1', 'tpx', 'tpy', 'tpz'): print(k, np.percentile(c[k][m], [0.1, 1, 16, 50, 84, 99, 99.9]))
    m = c['nt'] == 1
    for lo, hi in [(0, .01), (.01, .02), (.02, .025), (.025, .03), (.03, .035), (.035, .04), (.04, .05), (.05, .06), (.06, .07), (.07, .08), (.08, .09), (.09, .095), (.095, .10), (.10, .11), (.11, .12), (.12, .15), (.15, .25)]:
        s_ = f"{lo:.3f}-{hi:.3f}"
        for s in (0, 1):
            q = m & (c['txi1'] >= lo) & (c['txi1'] < hi) & (c['side'] == s)
            s_ += f" | s{s} N={q.sum():7d} near={np.mean(c['n_near'][q] > 0):.3f} far={np.mean(c['n_far'][q] > 0):.3f} prot={np.mean(c['n_prot'][q] > 0):.3f} tag={np.mean((c['n_prot'][q] == 1) & (c['n_win'][q] == 1)):.3f}"
        print(s_)
    q = (c['n_prot'] == 1) & (c['nt'] == 1)
    print('xi_reco - xi_truth', np.percentile((c['xi1'] - c['txi1'])[q], [2.5, 16, 50, 84, 97.5]))
