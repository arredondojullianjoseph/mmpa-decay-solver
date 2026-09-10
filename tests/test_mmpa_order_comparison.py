"""
test_mmpa_order_comparison.py

Sweeps every published MMPA order (Tables 2-4 of Chiba, Yamamoto, Nagano
2026) against Bateman on both chains, and checks that accuracy improves
with order.
"""

import numpy as np
import pytest

from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.mmpa import mmpa_linear_chain, MMPA_COEFFS
from data.four_isotope_chain import lambdas as l4, n0 as n0_4, t_grid as t_4
from data.gd157_chain import lambdas as l_gd, n0 as n0_gd, t_grid as t_gd

SIGNIFICANT_THRESHOLD = 1e-6
ORDER_32_GATE = 1e-8

CHAINS = {
    "four_isotope": (l4, n0_4, t_4),
    "gd157": (l_gd, n0_gd, t_gd),
}

PUBLISHED_ORDERS = sorted(MMPA_COEFFS)  

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
def test_mmpa_error_decreases_with_order(chain_name):
    """
    Runs every published order and checks that error drops as order
    goes up. 
    """
    lambdas, n0, t = CHAINS[chain_name]
    a_matrix = build_linear_chain_matrix(lambdas)

    errors = []
    for order in PUBLISHED_ORDERS:
        err = _max_rel_err(a_matrix, lambdas, n0, t, order=order)
        errors.append(err)
        print(f"\norder {order:2d} vs Bateman, {chain_name} chain: max relative error = {err:.3e}")

    for i in range(1, len(errors)):
        assert errors[i] < errors[i - 1], (
            f"order {PUBLISHED_ORDERS[i]} ({errors[i]:.3e}) is not more accurate "
            f"than order {PUBLISHED_ORDERS[i - 1]} ({errors[i - 1]:.3e})"
        )
