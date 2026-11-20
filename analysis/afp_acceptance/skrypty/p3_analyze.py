#!/usr/bin/env python3
"""p3: analiza tabel p1 (na lxplus, wynik = mały JSON). Bez mas kandydatów D.
MC (fast-sim AFP, prawda protonu): akceptancja A(ξ, pT, φ), wydajność tagu fiducjalnego, migracje, wrażliwość na nachylenie t i kształt ξ.
Dane: tag-and-probe near/far, wymóg dokładnie jednego protonu, wydajność L1 AFP względem MBTS, rozkłady ξ/x/y/pT protonów.
Użycie: p3_analyze.py out.json plik1.npz [plik2.npz ...]"""
import sys, json, numpy as np
XIB = [0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.095, 0.10, 0.105, 0.11, 0.12]
PTB = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0]
TRIG = None
TB = [0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.095, 0.10]

def load(fn):
    d = np.load(fn); n = list(d['names']); r = d['rows']
    c = {k: r[:, i] for i, k in enumerate(n)}
    return c, dict(zip([str(x) for x in d['cnt_keys']], [int(x) for x in d['cnt_vals']])), [str(x) for x in d['trig']], len(d['files']), int(d['nfiles_all'])

def eff(k, n):
    k = float(k); n = float(n)
    if n <= 0: return [None, None, 0]
    p = k / n; return [round(p, 5), round(float(np.sqrt(max(p * (1 - p), 1e-12) / n)), 5), int(n)]

def hist(x, lo, hi, nb, w=None):
    h, _ = np.histogram(x, bins=nb, range=(lo, hi), weights=w); return [float(v) for v in h]

def mc(c):
    o = {}
    pt = np.hypot(c['tpx'], c['tpy']) / 1000.; phi = np.arctan2(c['tpy'], c['tpx']); xi = c['txi1']
    has = c['n_prot'] > 0; tag = (c['n_prot'] == 1) & (c['n_win'] == 1)
    for s in (0, 1):
        ms = c['side'] == s
        o[f'A_xi_pt_s{s}'] = [[eff(has[q].sum(), q.sum()) for q in [ms & (xi >= XIB[i]) & (xi < XIB[i + 1]) & (pt >= PTB[j]) & (pt < PTB[j + 1])
                                                                    for j in range(len(PTB) - 1)]] for i in range(len(XIB) - 1)]
        core = ms & (xi >= 0.035) & (xi < 0.09)
        o[f'A_phi_pt_core_s{s}'] = [[eff(has[q].sum(), q.sum()) for q in [core & (pt >= PTB[j]) & (pt < PTB[j + 1]) & (phi >= -np.pi + k * np.pi / 4) & (phi < -np.pi + (k + 1) * np.pi / 4)
                                                                          for k in range(8)]] for j in range(len(PTB) - 1)]
    tw = (xi >= 0.03) & (xi <= 0.10)
    o['pt2_mean_window'] = float(np.mean(pt[tw] ** 2)); o['b_slope_window_GeV-2'] = float(1 / np.mean(pt[tw] ** 2))
    o['truth_pt_hist_window_0_2GeV_40'] = hist(pt[tw], 0, 2, 40); o['truth_xi_hist_0_0.2_40'] = hist(xi, 0, 0.2, 40)
    o['frac_window_pt_gt'] = {str(v): float(np.mean(pt[tw] > v)) for v in (0.5, 0.8, 1.0, 1.5)}
    sel = {'all': np.ones(len(xi), bool), 'D_2_30': (c['extra'].astype(int) & 8) > 0, 'D0_3_30': (c['extra'].astype(int) & 1) > 0, 'Dp_3_30': (c['extra'].astype(int) & 2) > 0}

    def fid(m, lo, hi, ptmax, w=None):
        w = np.ones(len(xi)) if w is None else w
        t = m & (xi >= lo) & (xi <= hi) & (pt < ptmax)
        r_ = m & (c['n_prot'] == 1) & (c['xi1'] >= lo) & (c['xi1'] <= hi)
        nt = w[t].sum(); nr = w[r_].sum(); nrt = w[r_ & t].sum()
        st = float(np.sqrt((w[t] ** 2).sum())) / max(nt, 1e-9)
        return {'N_truth': float(nt), 'N_reco': float(nr), 'eff_reco_given_truth': float(nrt / max(nt, 1e-9)), 'f_nonfid_in_reco': float(1 - nrt / max(nr, 1e-9)),
                'C_reco_over_truth': float(nr / max(nt, 1e-9)), 'rel_stat': st}
    o['fid'] = {}
    for nm, m in sel.items():
        for (lo, hi) in [(0.03, 0.10), (0.035, 0.09), (0.035, 0.08), (0.04, 0.09)]:
            for ptmax in (99., 1.0):
                o['fid'][f'{nm}|{lo}-{hi}|pt<{ptmax}'] = fid(m, lo, hi, ptmax)
            for s in (0, 1):
                o['fid'][f'{nm}|{lo}-{hi}|side{s}'] = fid(m & (c['side'] == s), lo, hi, 99.)
    # wrażliwość modelowa C (okno nominalne 0,03–0,10 i wąskie 0,035–0,09): przeważanie nachylenia t i kształtu ξ
    o['reweight'] = {}
    for nm, w in [('b-3', np.exp(3 * pt ** 2)), ('b-2', np.exp(2 * pt ** 2)), ('b+2', np.exp(-2 * pt ** 2)), ('b+4', np.exp(-4 * pt ** 2)),
                  ('xi^-0.5', np.clip(xi, 1e-4, 1) ** -0.5), ('xi^+0.5', np.clip(xi, 1e-4, 1) ** 0.5), ('xi^-1', np.clip(xi, 1e-4, 1) ** -1.0), ('xi^+1', np.clip(xi, 1e-4, 1) ** 1.0)]:
        w = np.where(pt < 2.0, w, 0.)
        for (lo, hi) in [(0.03, 0.10), (0.035, 0.09)]:
            o['reweight'][f'{nm}|{lo}-{hi}'] = fid(sel['all'] & (pt < 2.0), lo, hi, 99., w)['C_reco_over_truth']
    for (lo, hi) in [(0.03, 0.10), (0.035, 0.09)]:
        o['reweight'][f'nom|{lo}-{hi}'] = fid(sel['all'] & (pt < 2.0), lo, hi, 99.)['C_reco_over_truth']
    q = (c['n_prot'] == 1)
    o['xi_res_pct_2.5_16_50_84_97.5'] = [float(v) for v in np.percentile((c['xi1'] - xi)[q], [2.5, 16, 50, 84, 97.5])]
    o['reco'] = dists(c, tag)
    base = np.ones(len(xi), bool); c = dict(c); c['xical_same'] = np.full(len(xi), 0.01, np.float32)
    tp(c, base, o); o['all1'] = alldists(c, base)
    return o

def dists(c, m):
    o = {'n': int(m.sum())}
    for s in (0, 1):
        q = m & (c['side'] == s)
        o[f's{s}'] = {'n': int(q.sum()), 'xi_0.02_0.12_20': hist(c['xi1'][q], 0.02, 0.12, 20), 'x_near_-20_0_40': hist(c['x_near'][q], -20, 0, 40), 'x_far_-20_0_40': hist(c['x_far'][q], -20, 0, 40),
                      'y_near_-10_20_30': hist(c['y_near'][q], -10, 20, 30), 'y_far_-10_20_30': hist(c['y_far'][q], -10, 20, 30),
                      'px_-2_2_40': hist(c['px'][q], -2, 2, 40), 'py_-2_2_40': hist(c['py'][q], -2, 2, 40), 'pt_0_2_40': hist(np.hypot(c['px'][q], c['py'][q]), 0, 2, 40),
                      'pt_mean': float(np.mean(np.hypot(c['px'][q], c['py'][q]))) if q.sum() else None,
                      'x_near_pct_1_50_99': [float(v) for v in np.percentile(c['x_near'][q], [1, 50, 99])] if q.sum() else None}
    return o

def alldists(c, base):
    o = {}
    for s in (0, 1):
        ms = base & (c['side'] == s)
        o[f's{s}'] = {'xi_1prot_0_0.2_80': hist(c['xi1'][ms & (c['n_prot'] == 1)], 0, 0.2, 80), 'x_near_1trk_-20_0_80': hist(c['x_near'][ms & (c['n_near'] == 1)], -20, 0, 80),
                      'x_far_1trk_-20_0_80': hist(c['x_far'][ms & (c['n_far'] == 1)], -20, 0, 80),
                      'x_near_1trk_far1_-20_0_80': hist(c['x_near'][ms & (c['n_near'] == 1) & (c['n_far'] == 1)], -20, 0, 80),
                      'x_far_1trk_near0_-20_0_80': hist(c['x_far'][ms & (c['n_far'] == 1) & (c['n_near'] == 0)], -20, 0, 80),
                      'x_far_1trk_near1_-20_0_80': hist(c['x_far'][ms & (c['n_far'] == 1) & (c['n_near'] == 1)], -20, 0, 80),
                      'y_near_1trk_-10_20_60': hist(c['y_near'][ms & (c['n_near'] == 1)], -10, 20, 60), 'y_far_1trk_-10_20_60': hist(c['y_far'][ms & (c['n_far'] == 1)], -10, 20, 60),
                      'xi_in_xnearbins': None}
    return o

def tp(c, base, o):
    one = base & (c['n_prot'] == 1)
    xie = {}
    for s in (0, 1):
        q = one & (c['side'] == s) & (c['n_near'] == 1) & (c['n_far'] == 1) & (c['xi1'] > 0.02) & (c['xi1'] < 0.12) & (c['x_near'] > -20) & (c['x_far'] > -20)   # tylko zakres okna, float64
        pn = np.polyfit(c['x_near'][q].astype(float), c['xi1'][q].astype(float), 2); pf = np.polyfit(c['x_far'][q].astype(float), c['xi1'][q].astype(float), 2)
        o[f'xi_of_x_poly2_s{s}'] = {'near': [float(v) for v in pn], 'far': [float(v) for v in pf],
                                    'resid_rms_near': float(np.std(c['xi1'][q] - np.polyval(pn, c['x_near'][q]))), 'resid_rms_far': float(np.std(c['xi1'][q] - np.polyval(pf, c['x_far'][q]))),
                                    'dx_far_minus_near_pct_1_16_50_84_99': [float(v) for v in np.percentile((c['x_far'] - c['x_near'])[q], [1, 16, 50, 84, 99])],
                                    'dy_pct': [float(v) for v in np.percentile((c['y_far'] - c['y_near'])[q], [1, 16, 50, 84, 99])]}
        xie[s] = (pn, pf)
    # krotności torów
    o['mult'] = {}
    for s in (0, 1):
        ms = base & (c['side'] == s)
        o['mult'][f's{s}'] = {'n_near': np.bincount(np.minimum(c['n_near'][ms], 6).astype(int), minlength=7).tolist(), 'n_far': np.bincount(np.minimum(c['n_far'][ms], 6).astype(int), minlength=7).tolist(),
                              'n_prot': np.bincount(np.minimum(c['n_prot'][ms], 6).astype(int), minlength=7).tolist(),
                              'nprot_given_nwin>=1': np.bincount(np.minimum(c['n_prot'][ms & (c['n_win'] >= 1)], 6).astype(int), minlength=7).tolist(),
                              'ncl_near_1trk': np.bincount(np.clip(c['ncl_near'][ms & (c['n_near'] == 1)], 0, 6).astype(int), minlength=7).tolist(),
                              'ncl_far_1trk': np.bincount(np.clip(c['ncl_far'][ms & (c['n_far'] == 1)], 0, 6).astype(int), minlength=7).tolist()}
    # tag-and-probe
    TB = [0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.095, 0.10]
    o['tp'] = {}
    for s in (0, 1):
        ms = base & (c['side'] == s); pn, pf = xie[s]
        for tagst, (ntag, xt, pp, nprobe, xp, yt, yp) in {'tagNEAR_probeFAR': ('n_near', 'x_near', pn, 'n_far', 'x_far', 'y_near', 'y_far'), 'tagFAR_probeNEAR': ('n_far', 'x_far', pf, 'n_near', 'x_near', 'y_far', 'y_near')}.items():
            xe = np.polyval(pp, c[xt])
            for var, extra in {'nom': np.ones(len(xe), bool), 'xical': c['xical_same'] <= xe + 0.02, 'xical_tight': (c['xical_same'] <= xe + 0.02) & (c['xical_same'] > 0.005),
                               'ncl4': (c['ncl_far' if ntag == 'n_far' else 'ncl_near'] >= 4) | (c['ncl_near'].max() <= 0)}.items():
                t = ms & (c[ntag] == 1) & extra
                def row(lo, hi):
                    q = t & (xe >= lo) & (xe < hi); n = q.sum()
                    return {'N': int(n), 'probe>=1': eff((c[nprobe][q] >= 1).sum(), n), 'probe==1': eff((c[nprobe][q] == 1).sum(), n),
                            'prot>=1': eff((c['n_prot'][q] >= 1).sum(), n), 'tag_exactly1_inwin': eff(((c['n_prot'][q] == 1) & (c['n_win'][q] == 1)).sum(), n),
                            'prot==1': eff((c['n_prot'][q] == 1).sum(), n),
                            'probe_n': np.bincount(np.minimum(c[nprobe][q], 5).astype(int), minlength=6).tolist()}
                o['tp'][f's{s}|{tagst}|{var}'] = {'bins': {f'{TB[i]}-{TB[i + 1]}': row(TB[i], TB[i + 1]) for i in range(len(TB) - 1)}, '0.035-0.09': row(0.035, 0.09), '0.03-0.10': row(0.03, 0.10), '0.04-0.08': row(0.04, 0.08)}
    return o

def data(c, cnt, trig):
    o = {'cnt': cnt}
    ex = c['extra'].astype(np.int64); bit = lambda nm: ((ex >> trig.index(nm)) & 1) == 1
    ref = 'trigPassed_HLT_noalg_mb_L1MBTS_1' if cnt.get('grl_trigPassed_HLT_noalg_mb_L1MBTS_1', 0) > 0.5 * cnt['n_grl'] else 'trigPassed_HLT_noalg_L1J12'   # 428770: odniesienie J12
    o['ref_trigger'] = ref
    base = (c['grl'] == 1) & (c['npu'] == 0) & bit(ref)
    o['n_base_rows'] = int(base.sum())
    one = base & (c['n_prot'] == 1)
    tag = one & (c['n_win'] == 1); t2 = tag & (c['xical_same'] <= c['xi1'] + 0.02)
    o['n_tag'] = int(tag.sum()); o['n_t2'] = int(t2.sum())
    o['reco_t2'] = dists(c, t2)
    tp(c, base, o)
    o['all1'] = alldists(c, base)
    # wydajność L1 AFP względem MBTS dla tagów offline (T2)
    o['l1'] = {}
    for s, (b, fs, ns) in {0: ('trigTBP_L1_AFP_A', 'trigTBP_L1_AFP_FSA_BGRP12', 'trigTBP_L1_AFP_NSA_BGRP12'), 1: ('trigTBP_L1_AFP_C', 'trigTBP_L1_AFP_FSC_BGRP12', 'trigTBP_L1_AFP_NSC_BGRP12')}.items():
        q0 = t2 & (c['side'] == s)
        def row(q):
            n = q.sum(); return {'N': int(n), 'L1_AFP_side': eff(bit(b)[q].sum(), n), 'L1_FS': eff(bit(fs)[q].sum(), n), 'L1_NS': eff(bit(ns)[q].sum(), n), 'L1_A_OR_C': eff(bit('trigTBP_L1_AFP_A_OR_C')[q].sum(), n),
                                 'HLT_noalg_L1AFP': eff(bit('trigPassed_HLT_noalg_L1AFP_A_OR_C')[q].sum(), n)}
        o['l1'][f's{s}'] = {'all': row(q0), 'bins': {f'{TB[i]}-{TB[i + 1]}': row(q0 & (c['xi1'] >= TB[i]) & (c['xi1'] < TB[i + 1])) for i in range(len(TB) - 1)},
                            'lb_quart': [row(q0 & (c['lb'] >= a) & (c['lb'] < b_)) for a, b_ in zip(*(lambda e: (e[:-1], e[1:]))(np.quantile(c['lb'][q0], [0, .25, .5, .75, 1]) + np.array([0, 0, 0, 0, 1])))],
                            'y_bins_near': {f'{a}-{b_}': row(q0 & (c['y_near'] >= a) & (c['y_near'] < b_)) for a, b_ in [(-10, 0), (0, 2), (2, 4), (4, 6), (6, 10), (10, 20)]}}
    # MBTS względem AFP: zdarzenia z HLT_noalg_L1AFP (preskala 1 w GRL 427929) i tagiem T2 — ułamek z bitem L1_MBTS_1
    aa = (c['grl'] == 1) & (c['npu'] == 0) & bit('trigPassed_HLT_noalg_L1AFP_A_OR_C') & (c['n_prot'] == 1) & (c['n_win'] == 1) & (c['xical_same'] <= c['xi1'] + 0.02)
    o['mbts_given_afp_tag'] = {'L1_MBTS_1': eff(bit('trigTBP_L1_MBTS_1')[aa].sum(), aa.sum()), 'L1_MBTS_2': eff(bit('trigTBP_L1_MBTS_2')[aa].sum(), aa.sum()), 'HLT_mb_mbts': eff(bit('trigPassed_HLT_mb_mbts_L1MBTS_1')[aa].sum(), aa.sum())}
    # odniesienie J12 (428770)
    jj = (c['grl'] == 1) & (c['npu'] == 0) & bit('trigPassed_HLT_noalg_L1J12') & (c['n_prot'] == 1) & (c['n_win'] == 1) & (c['xical_same'] <= c['xi1'] + 0.02)
    o['l1_ref_J12'] = {f's{s}': {'N': int((jj & (c['side'] == s)).sum()), 'L1_AFP_side': eff(bit('trigTBP_L1_AFP_A' if s == 0 else 'trigTBP_L1_AFP_C')[jj & (c['side'] == s)].sum(), (jj & (c['side'] == s)).sum()),
                                 'L1_A_OR_C': eff(bit('trigTBP_L1_AFP_A_OR_C')[jj & (c['side'] == s)].sum(), (jj & (c['side'] == s)).sum())} for s in (0, 1)}
    return o

if __name__ == "__main__":
    out = {}
    for fn in sys.argv[2:]:
        c, cnt, trig, nf, nfa = load(fn)
        key = fn.split('/')[-1].replace('.npz', '')
        r = mc(c) if 'nt' in c else data(c, cnt, trig)
        r['files'] = [nf, nfa]; r['cnt'] = cnt; out[key] = r
        print(key, 'ok', flush=True)
    json.dump(out, open(sys.argv[1], 'w'))
