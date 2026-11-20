#!/usr/bin/env python3
"""SFABS metoda C: czulosc SF(WP70) drzewa H na skale rozdzielczosci wejsc i dE/dx dla D0 i D+ per bin (sygnal fd_nd z prawda; modele rozmycia T i F z INCAL inc_lib — import bez zmian).
Wyjscie: out/c_sens.json: SF_var per gatunek|bin|wariant: [z migracja przez prog L_xy, bez migracji], 2 ziarna."""
import json, sys, numpy as np
sys.path.insert(0, '/home/k/Documents/charm/charm_work/INCAL_2026-11-19/scripts'); import inc_lib as IL
P = '/home/k/Documents/charm/charm_work/INCAL_2026-11-19/runs/lxplus/parts'
BINS = ['3-5', '5-8', '8-30']; res = {}
for sp in ('D0', 'Dplus'):
    T = IL.load(P, 'fd_nd', 'd150', sp, 'sig'); G = IL.geom(T, sp)
    sc0 = IL.score(sp, T['uid'], IL.feats(T, sp))
    thr = {}
    for b in BINS:
        lo, hi = IL.PTB[b]; m = IL.sigmask(T, sp) & (T['pt'] >= lo) & (T['pt'] < hi); thr[b] = float(np.quantile(sc0[m], 0.30))
    V = {'T1.02': lambda r: IL.smear_T(T, sp, 1.02, 1.02, r, G), 'T1.05': lambda r: IL.smear_T(T, sp, 1.05, 1.05, r, G), 'T1.10': lambda r: IL.smear_T(T, sp, 1.10, 1.10, r, G),
         'F1.02': lambda r: IL.smear_F(T, sp, 1.02, 1.02, r, G), 'F1.05': lambda r: IL.smear_F(T, sp, 1.05, 1.05, r, G), 'F1.10': lambda r: IL.smear_F(T, sp, 1.10, 1.10, r, G),
         'dEdx0.95': lambda r: IL.smear_dedx(T, sp, 0.95, 0.0, r), 'dEdx1.05': lambda r: IL.smear_dedx(T, sp, 1.05, 0.0, r), 'dEdxS5': lambda r: IL.smear_dedx(T, sp, 1.0, 0.05, r)}
    for name, fn in V.items():
        for seed in (1, 2):
            O = fn(np.random.default_rng(seed)); sc = IL.score(sp, O['uid'], IL.feats(O, sp))
            for b in BINS:
                e0, n0 = IL.eff(T, sp, sc0, thr[b], b); e1, n1 = IL.eff(O, sp, sc, thr[b], b); e1b, _ = IL.eff(O, sp, sc, thr[b], b, w=T['Lxy'] >= IL.LXY_NOM[sp])
                res.setdefault(f'{sp}|{b}|{name}', []).append([e1 / e0, e1b / e0, n0])
        print(sp, name, ' '.join(f"{b}: {np.mean([x[0] for x in res[f'{sp}|{b}|{name}']]):.4f}/{np.mean([x[1] for x in res[f'{sp}|{b}|{name}']]):.4f} (Δziarna {abs(res[f'{sp}|{b}|{name}'][0][0] - res[f'{sp}|{b}|{name}'][1][0]):.4f})" for b in BINS), flush=True)
json.dump(res, open('/home/k/Documents/charm/charm_work/SFABS_wynik_2026-11-20/out/c_sens.json', 'w'), indent=1)
