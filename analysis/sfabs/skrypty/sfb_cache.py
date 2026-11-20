"""SFABS B/C: skrót tabel par D0<-D* z części INCAL do jednego pliku per próbka (kolumny potrzebne do fitu 2D)."""
import sys, os, numpy as np
sys.path.insert(0, '/home/k/Documents/charm/charm_work/INCAL_2026-11-19/scripts'); import inc_lib as IL
P = '/home/k/Documents/charm/charm_work/INCAL_2026-11-19/runs/lxplus/parts'; OUT = sys.argv[1]
for s in ('mb_nd', 'fd_nd', 'data_427929', 'data_428770'):
    T = IL.load(P, s, 'd150', 'D0', 'pair')
    w = (T['Lxy'] >= 0.25) & T['ok'] & (T['ds_dm'] < 0.156)
    cols = ['m', 'ds_dm', 'pt', 'score', 'cls', 'swap', 'ds_truth', 'uid', 'run'] + [c for c in ('truth', 'anc') if c in T]
    np.savez(os.path.join(OUT, f'pairs_{s}.npz'), **{c: T[c][w] for c in cols})
    print(s, T['_n'], int(w.sum()), {c: str(T[c].dtype) for c in cols if c in ('truth', 'anc')}, np.unique(T['cls'][w], return_counts=True))
