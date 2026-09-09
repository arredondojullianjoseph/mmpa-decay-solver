"""
test_mmpa_order_comparison.py

Compares MMPA order 16 against order 32 on both chains.
"""

import numpy as np
import pytest

from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.mmpa import mmpa_linear_chain
from data.four_isotope_chain import lambdas as l4, n0 as n0_4, t_grid as t_4
from data.gd157_chain import lambdas as l_gd, n0 as n0_gd, t_grid as t_gd

SIGNIFICANT_THRESHOLD = 1e-6

CHAINS = {
    "four_isotope": (l4, n0_4, t_4),
    "gd157": (l_gd, n0_gd, t_gd),
}

ORDER_GATES = {16: 1e-2, 32: 1e-8}

@pytest.mark.parametrize("chain_name", CHAINS)
@pytest.mark.parametrize("order", [16, 32])
def test_mmpa_order_accuracy_on_chain(chain_name, order):
    lambdas, n0, t = CHAINS[chain_name]
    a_matrix = build_linear_chain_matrix(lambdas)
    bateman_ref = bateman_linear_chain(lambdas, t, n0_parent=n0[0])
    mmpa_result = mmpa_linear_chain(a_matrix, t, n0, order=order)

    significant = np.abs(bateman_ref) >= SIGNIFICANT_THRESHOLD
    assert significant.any(), "mask is empty; nothing was actually checked"

    rel_err = np.abs(mmpa_result[significant] - bateman_ref[significant]) / np.abs(bateman_ref[significant])
    assert rel_err.max() < ORDER_GATES[order]
