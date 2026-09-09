"""
test_radau_vs_bateman.py

Checks the Radau reference solver against the analytical Bateman solution.
"""

import numpy as np

from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.radau import radau_linear_chain
from data.four_isotope_chain import lambdas as l4, n0 as n0_4, t_grid as t_4
from data.gd157_chain import lambdas as l_gd, n0 as n0_gd, t_grid as t_gd

SIGNIFICANT_THRESHOLD = 1e-6
GATE_RTOL = 1e-8 

def test_radau_matches_bateman_on_four_isotope_chain():
    a_matrix = build_linear_chain_matrix(l4)
    bateman_ref = bateman_linear_chain(l4, t_4, n0_parent=n0_4[0])
    radau_result = radau_linear_chain(a_matrix, t_4, n0_4)
    significant = np.abs(bateman_ref) >= SIGNIFICANT_THRESHOLD
    assert significant.any(), "mask is empty; nothing was actually checked"
    np.testing.assert_allclose(radau_result[significant], bateman_ref[significant], rtol=GATE_RTOL)

def test_radau_matches_bateman_on_gd157_chain():
    a_matrix = build_linear_chain_matrix(l_gd)
    bateman_ref = bateman_linear_chain(l_gd, t_gd, n0_parent=n0_gd[0])
    radau_result = radau_linear_chain(a_matrix, t_gd, n0_gd)
    significant = np.abs(bateman_ref) >= SIGNIFICANT_THRESHOLD
    assert significant.any(), "mask is empty; nothing was actually checked"
    np.testing.assert_allclose(radau_result[significant], bateman_ref[significant], rtol=GATE_RTOL)

def test_radau_conserves_atoms_on_four_isotope_chain():
    a_matrix = build_linear_chain_matrix(l4)
    radau_result = radau_linear_chain(a_matrix, t_4, n0_4)
    totals = radau_result.sum(axis=0)
    np.testing.assert_allclose(totals, n0_4[0], rtol=1e-6)

def test_radau_rejects_non_increasing_t():
    a_matrix = build_linear_chain_matrix(l4)
    with np.testing.assert_raises(ValueError):
        radau_linear_chain(a_matrix, [0.0, 1.0, 0.5], n0_4)
