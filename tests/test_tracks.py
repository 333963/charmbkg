"""ID śladu (d0, z0) w charmbkg.tracks musi być identyczny z charmtag.pairs.build (projekt A) — kopia, nie nowa definicja."""
import numpy as np
import pytest

from charmbkg.tracks import trk_id_d0z0


def test_trk_id_nonnegative_and_masked():
    rng = np.random.default_rng(3)
    d0, z0 = rng.normal(0, 0.1, 1000).astype(np.float32), rng.normal(0, 50, 1000).astype(np.float32)
    pt = np.where(rng.random(1000) < 0.1, 0.0, 1.0).astype(np.float32)
    t = trk_id_d0z0(pt, d0, z0)
    assert (t[pt > 0] >= 0).all() and (t[pt <= 0] == -1).all()


def test_same_as_charmtag():
    build = pytest.importorskip("charmtag.pairs.build", reason="integracja z charmtag: pip install charmtag")
    rng = np.random.default_rng(4)
    pt, d0, z0 = (rng.normal(0, s, 5000).astype(np.float32) for s in (5, 0.1, 50))
    np.testing.assert_array_equal(trk_id_d0z0(pt, d0, z0), build._trk_id_d0z0(pt, d0, z0))
