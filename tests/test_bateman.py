"""
test_bateman.py

This file contains the automated tests for our Bateman reference solution.
We're checking it against known expected results—like simple exponential decay,
the standard two-nuclide formula from textbooks, and basic atom conservation.
We also verify that the input validation catches bad inputs properly.
Note: The MMPA stuff isn't included here since it hasn't been built yet.
"""

import numpy as np
import pytest

from src.bateman import bateman_linear_chain
from data.four_isotope_chain import lambdas as chain_lambdas, n0 as chain_n0, t_grid as chain_t_grid

rel_tol = 1e-10

def test_initial_condition_at_t_zero():
    # At time zero, all the atoms should be in the parent nuclide, and the daughters should be empty.
    result = bateman_linear_chain(chain_lambdas, t=0.0, n0_parent=chain_n0[0])
    expected = np.array([[chain_n0[0]], [0.0], [0.0], [0.0]])
    np.testing.assert_allclose(result, expected, atol=1e-12)

def test_parent_is_simple_exponential():
    # The parent nuclide doesn't have anything feeding into it, so it just undergoes
    # standard exponential decay. This is true no matter how long or complex the rest of the chain is.
    lambdas = [0.7, 0.3, 0.1]
    t = np.array([0.0, 0.5, 1.0, 5.0])
    result = bateman_linear_chain(lambdas, t, n0_parent=2.0)
    expected_parent = 2.0 * np.exp(-0.7 * t)
    np.testing.assert_allclose(result[0, :], expected_parent, rtol=rel_tol)

def test_two_nuclide_textbook_formula():
    # Checking against the two-member chain formula 
    l1, l2 = 0.8, 0.2
    n0 = 3.0
    t = np.linspace(0.0, 10.0, 25)
    result = bateman_linear_chain([l1, l2], t, n0_parent=n0)
    expected_n1 = n0 * np.exp(-l1 * t)
    expected_n2 = n0 * l1 / (l2 - l1) * (np.exp(-l1 * t) - np.exp(-l2 * t))
    np.testing.assert_allclose(result[0, :], expected_n1, rtol=rel_tol)
    np.testing.assert_allclose(result[1, :], expected_n2, rtol=rel_tol)

def test_atom_conservation_with_stable_end_member():
    # If there's no branching and the last nuclide is completely stable, the atoms just move
    # straight down the line. Because of this, the total number of atoms across the whole chain shouldn't ever change.
    result = bateman_linear_chain(chain_lambdas, chain_t_grid, n0_parent=chain_n0[0])
    total = result.sum(axis=0)
    np.testing.assert_allclose(total, chain_n0[0], rtol=1e-8)

def test_parent_is_nonincreasing():
    # The parent nuclide should only ever decay, meaning its atom count shouldn't go up over time.
    result = bateman_linear_chain(chain_lambdas, chain_t_grid, n0_parent=chain_n0[0])
    parent = result[0, :]
    assert np.all(np.diff(parent) <= 1e-12)

def test_stable_daughter_is_nondecreasing():
    # The last nuclide (N4) is stable, so it just collects atoms as the rest of the chain decays into it. So, its atom count should never drop.
    result = bateman_linear_chain(chain_lambdas, chain_t_grid, n0_parent=chain_n0[0])
    stable_daughter = result[-1, :]
    assert np.all(np.diff(stable_daughter) >= -1e-12)

def test_single_nuclide_chain():
    # Testing the edge case where the chain only has one nuclide (k=1). It should just be simple exponential decay.
    t = np.array([0.0, 1.0, 2.0])
    result = bateman_linear_chain([0.4], t, n0_parent=5.0)
    expected = 5.0 * np.exp(-0.4 * t)
    np.testing.assert_allclose(result[0, :], expected, rtol=rel_tol)
    assert result.shape == (1, 3)

@pytest.mark.parametrize("lambdas, t", [
    ([0.5, 0.5, 0.1], 1.0),  # Two of the same decay constants (causes a divide-by-zero error)
    ([0.5, -0.1], 1.0),      # Negative decay constant
    ([0.5, 0.1], -1.0),      # Negative time value
    ([], 1.0),               # Empty list of decay constants
])
def test_rejects_invalid_input(lambdas, t):
    # This makes sure our function actually catches bad inputs and throws a ValueError when it should.
    with pytest.raises(ValueError):
        bateman_linear_chain(lambdas, t)
