"""
test_mmpa_multi_time.py

One factorization of (A Δt - cI), applied K times, vs Bateman at t = Δt … KΔt.
"""

import numpy as np
import pytest

from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src import mmpa
from data.four_isotope_chain import lambdas as l4, n0 as n0_4
from data.gd157_chain import CHARACTERISTIC_TIME, lambdas as l_gd, n0 as n0_gd

SIGNIFICANT_THRESHOLD = 1e-6
GATE = 1e-8


def _check_multi_time(lambdas, n0, dt, n_steps):
    a_matrix = build_linear_chain_matrix(lambdas)
    got = mmpa.mmpa_multi_time(a_matrix, dt, n0, n_steps)
    times = dt * np.arange(1, n_steps + 1)
    ref = bateman_linear_chain(lambdas, times, n0_parent=n0[0])
    significant = np.abs(ref) >= SIGNIFICANT_THRESHOLD
    assert significant.any(), "mask is empty; nothing was actually checked"
    rel_err = np.abs(got[:, 1:][significant] - ref[significant]) / np.abs(ref[significant])
    max_rel_err = rel_err.max()
    print(f"\norder 32 mmpa_multi_time vs Bateman: max relative error = {max_rel_err:.3e}")
    assert max_rel_err < GATE


def test_multi_time_four_isotope():
    _check_multi_time(l4, n0_4, dt=0.5, n_steps=20)


def test_multi_time_gd157():
    _check_multi_time(l_gd, n0_gd, dt=0.1 * CHARACTERISTIC_TIME, n_steps=10)


def test_factors_once_for_many_steps(monkeypatch):
    calls = {"n": 0}
    real = mmpa.lu_factor

    def counted(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    monkeypatch.setattr(mmpa, "lu_factor", counted)
    a_matrix = build_linear_chain_matrix(l4)
    mmpa.mmpa_multi_time(a_matrix, 0.5, n0_4, n_steps=20)
    assert calls["n"] == 1


def test_zero_steps_returns_initial():
    a_matrix = build_linear_chain_matrix(l4)
    got = mmpa.mmpa_multi_time(a_matrix, 0.5, n0_4, n_steps=0)
    assert got.shape == (4, 1)
    np.testing.assert_array_equal(got[:, 0], n0_4)


def test_negative_steps_rejected():
    a_matrix = build_linear_chain_matrix(l4)
    with pytest.raises(ValueError, match="n_steps"):
        mmpa.mmpa_multi_time(a_matrix, 0.5, n0_4, n_steps=-1)