"""Identyfikatory śladów córek niezależne od refitu kandydata (PF).

Kopia 1:1 `charmtag.pairs.build._trk_id_d0z0` z charmtag 2.0.3 — tu, żeby `charmbkg.hyp` nie zależał od kodu par
(projekt A). Zmiana definicji ID śladu musi być wprowadzona w obu repo (test zgodności: tests/test_tracks.py).
"""
import numpy as np


def trk_id_d0z0(pt, d0, z0):
    """Identyfikator śladu z parametrów toru względem PV (d0, z0·sinθ) — niezależnych od refitu kandydata (PF):
    na danych v1.1 100 % par z tolerancją 1 % w (pT,η,φ) ma identyczne bity (d0, z0), a 3,7 % par dzieli ślad przy różnicy pT > 1 %."""
    pt = np.asarray(pt, np.float32)
    d0 = np.asarray(d0, np.float32)
    z0 = np.asarray(z0, np.float32)
    # PF3-FIN (05.11): `<< 32` na int64 dawało UJEMNE ID dla d0 < 0 (bit znaku float32) → shared_track/arbitraż (t >= 0)
    # pomijały ≈ połowę torów (weto działało w połowie: D⁺D⁺ 68 % zamiast 85 %). Teraz 63 bity, zawsze ≥ 0: pełne 32 bity d0
    # (ze znakiem) << 31 | z0 bez najmłodszego bitu mantysy (kolizja wymaga bitowo identycznego d0 i z0 różnego o 1 ulp).
    tid = (d0.view(np.uint32).astype(np.int64) << 31) | (z0.view(np.uint32).astype(np.int64) >> 1)
    return np.where(pt > 0, tid, -1)
