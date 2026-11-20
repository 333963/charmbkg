"""samples.pt_range: okno pT kandydata — domyślne 3–30 GeV, nadpisywane przez CHARMBKG_PT_RANGE (dane syntetyczne)."""
import pandas as pd
import pytest

from charmbkg import samples as S

SP = {"mass_range": [1.7, 2.1]}


def _frame():
    return pd.DataFrame({"m": [1.86, 1.86, 1.86, 1.86, 2.5], "pt": [2.5, 3.0, 29.9, 30.0, 5.0]})


def test_default_window(monkeypatch):
    monkeypatch.delenv("CHARMBKG_PT_RANGE", raising=False)
    assert S.pt_range() == (3.0, 30.0)
    assert S.in_range(_frame(), SP).tolist() == [False, True, True, False, False]


def test_env_override(monkeypatch):
    monkeypatch.setenv("CHARMBKG_PT_RANGE", "2,3")
    assert S.pt_range() == (2.0, 3.0)
    assert S.in_range(_frame(), SP).tolist() == [True, False, False, False, False]


@pytest.mark.parametrize("bad", ["5,3", "-1,3", "3"])
def test_env_invalid(monkeypatch, bad):
    monkeypatch.setenv("CHARMBKG_PT_RANGE", bad)
    with pytest.raises(ValueError):
        S.pt_range()
