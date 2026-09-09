"""
test_mmpa_order_comparison.py

Compares MMPA order 16 against order 32 on both chains, using the coefficient tables already in src/mmpa.py.
"""

import numpy as np
import pytest

from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.mmpa import mmpa_linear_chain
from data.four_isotope_chain import lambdas as l4, n0 as n0_4, t_grid as t_4
from data.gd157_chain import lambdas as l_gd, n0 as n0_gd, t_grid as t_gd

SIGNIFICANT_THRESHOLD = 1e-6
ORDER_32_GATE = 1e-8  

CHAINS = {
    "four_isotope": (l4, n0_4, t_4),
    "gd157": (l_gd, n0_gd, t_gd),
}

def _max_rel_err(a_matrix, lambdas, n0, t, order):
    bateman_ref = bateman_linear_chain(lambdas, t, n0_parent=n0[0])
    mmpa_result = mmpa_linear_chain(a_matrix, t, n0, order=order)
    significant = np.abs(bateman_ref) >= SIGNIFICANT_THRESHOLD
    assert significant.any(), "mask is empty; nothing was actually checked"
    rel_err = np.abs(mmpa_result[significant] - bateman_ref[significant]) / np.abs(bateman_ref[significant])
    return rel_err.max()

@pytest.mark.parametrize("chain_name", CHAINS)
def test_mmpa_order_32_meets_gate(chain_name):
    lambdas, n0, t = CHAINS[chain_name]
    a_matrix = build_linear_chain_matrix(lambdas)
    max_rel_err = _max_rel_err(a_matrix, lambdas, n0, t, order=32)
    print(f"\norder 32 vs Bateman, {chain_name} chain: max relative error = {max_rel_err:.3e}")
    assert max_rel_err < ORDER_32_GATE

@pytest.mark.parametrize("chain_name", CHAINS)
def test_mmpa_order_16_is_less_accurate_than_order_32(chain_name):
    
    #Checks the ordering,
    lambdas, n0, t = CHAINS[chain_name]
    a_matrix = build_linear_chain_matrix(lambdas)
    err_16 = _max_rel_err(a_matrix, lambdas, n0, t, order=16)
    err_32 = _max_rel_err(a_matrix, lambdas, n0, t, order=32)
    print(f"\norder 16 vs Bateman, {chain_name} chain: max relative error = {err_16:.3e}")
    print(f"order 32 vs Bateman, {chain_name} chain: max relative error = {err_32:.3e}")
    assert err_16 > err_32 #Order 16 uses fewer terms, so it should be worse than order 32 on the same chain.
