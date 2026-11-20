import numpy as np


def event_mult(ds, ev, chunk=4_000_000):
    """Per wiersz: n_cand (liczba wierszy w bloku kolejnych wierszy o tym samym event_number) i n_trkid (różne trackId córek)."""
    n = len(ev)
    brk = np.r_[0, np.flatnonzero(ev[1:] != ev[:-1]) + 1, n]
    blk = np.repeat(np.arange(len(brk) - 1), np.diff(brk))
    ncand = np.diff(brk)[blk]
    tid = np.stack([ds.fields([f"Dmeson_daughter{k}_trackId"])[:][f"Dmeson_daughter{k}_trackId"] for k in range(3)], 1).astype(np.int64)
    b = np.repeat(blk, 3); t = tid.ravel(); ok = t >= 0
    key = b[ok].astype(np.int64) * (1 << 24) + t[ok]
    u = np.unique(key)
    ntr = np.bincount((u >> 24).astype(np.int64), minlength=len(brk) - 1)
    return {"n_cand": ncand.astype(np.int32), "n_trkid": ntr[blk].astype(np.int32)}
