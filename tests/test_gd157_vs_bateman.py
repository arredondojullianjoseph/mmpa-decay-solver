"""
test_gd157_vs_bateman.py

Bateman validation for the Gd-157 absorber chain. 
"""

import numpy as np
import pytest

from data.gd157_chain import lambdas, n0, t_grid
from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.mmpa import mmpa_linear_chain

TOLERANCE = 1e-10


def hand_derived_solution(t):
    """
    Closed form for parent to stable daughter.
    """
    lambda1 = lambdas[0]
    n1 = n0[0] * np.exp(-lambda1 * t)
    n2 = n0[0] * (1.0 - np.exp(-lambda1 * t))
    return np.array([n1, n2])


@pytest.mark.parametrize("t", [0.0, 0.0001, 0.001, 0.01, 0.1, 1.0])
def test_matches_hand_derived_closed_form(t):
    result = bateman_linear_chain(lambdas, t, n0_parent=n0[0])[:, 0]
    expected = hand_derived_solution(t)
    mask = expected > 1e-12
    relative_error = np.abs(result - expected) / expected
    assert np.max(relative_error[mask]) < TOLERANCE


@pytest.mark.parametrize("t", [0.0, 0.01, 0.1, 1.0, 10.0])
def test_mass_conservation(t):
    result = bateman_linear_chain(lambdas, t, n0_parent=n0[0])[:, 0]
    assert abs(np.sum(result) - np.sum(n0)) < TOLERANCE


def test_parent_decays_to_zero():
    result = bateman_linear_chain(lambdas, 1.0, n0_parent=n0[0])[:, 0]
    assert result[0] < TOLERANCE


def test_daughter_approaches_parent_initial_density():
    result = bateman_linear_chain(lambdas, 1.0, n0_parent=n0[0])[:, 0]
    assert abs(result[1] - n0[0]) < TOLERANCE


def test_mmpa_matches_bateman_on_gd157_chain():
    a = build_linear_chain_matrix(lambdas)
    ref = bateman_linear_chain(lambdas, t_grid, n0_parent=n0[0])
    got = mmpa_linear_chain(a, t_grid, n0)
    significant = np.abs(ref) >= 1e-6
    assert significant.any()
    np.testing.assert_allclose(got[significant], ref[significant], rtol=1e-8)
