#!/usr/bin/env python3
"""SIGEMB: tabela dawców — (pT, η) prawdziwych D0←D* z MC fd_nd po selekcji H06 (jak ISO_origin). Uruchamiane LOKALNIE (wejście: B_trig/runs/lxplus/bt_fd_nd.npz)."""
import os, sys
import numpy as np
bt = sys.argv[1]; out = sys.argv[2]
z = np.load(bt)
s = np.isin(z["truth"], [0, 1]) & np.isin(z["truth_ds"], [0, 1]) & (np.abs(z["m"] - 1.86484) < 0.040) & (z["dm"] > 0.1435) & (z["dm"] < 0.1475) & (z["score"] > 0.6) & (z["pt"] >= 3) & (z["pt"] < 30)
np.savez_compressed(out, pt=z["pt"][s].astype(np.float32), eta=z["eta"][s].astype(np.float32))
print("dawców", int(s.sum()), "→", out, os.path.getsize(out))
