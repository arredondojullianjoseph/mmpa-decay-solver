"""
test_mmpa_vs_bateman.py

Verifies MMPA against the analytical Bateman solution 
"""

import numpy as np

from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.mmpa import mmpa_linear_chain
from data.four_isotope_chain import lambdas, n0, t_grid

SIGNIFICANT_THRESHOLD = 1e-6
TINY_ABS_TOLERANCE = 1e-6

def test_mmpa_matches_bateman_on_significant_inventory():
  """
  Chekks MMPA agree with the analytical solution wherever the answer is still physically meaningful.
  """

    a_matrix = build_linear_chain_matrix(lambdas)
    bateman_ref = bateman_linear_chain(lambdas, t_grid, n0_parent=1.0)
    mmpa_result = mmpa_linear_chain(a_matrix, t_grid, n0)
    significant = np.abs(bateman_ref) >= SIGNIFICANT_THRESHOLD
    assert significant.any(), "mask is empty; nothing was actually checked"
    np.testing.assert_allclose(
        mmpa_result[significant], bateman_ref[significant], rtol=1e-8
    ) #every MMPA value must be within relative tolerance 10e−8 of its Bateman counterpart

#Wherever Bateman's true value is below 10e-6 this only checks the raw absolute error, not a relative one
    abs_err = np.abs(mmpa_result - bateman_ref)
    tiny_abs_err = abs_err[~significant]
    if tiny_abs_err.size: 
        assert tiny_abs_err.max() < TINY_ABS_TOLERANCE

def test_mmpa_conserves_atoms():

    a_matrix = build_linear_chain_matrix(lambdas)
    mmpa_result = mmpa_linear_chain(a_matrix, t_grid, n0)
    totals = mmpa_result.sum(axis=0)
    np.testing.assert_allclose(totals, n0[0], rtol=1e-8)
