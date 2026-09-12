"""
absorber_sweep.py

Core sweep logic for the Gd-157 absorber interval-length test.
"""

import numpy as np

from data.gd157_chain import lambdas, n0
from src.bateman import bateman_linear_chain
from src.chain import build_linear_chain_matrix
from src.mmpa import DEFAULT_ORDER, mmpa_linear_chain

TOLERANCE = 1e-4
SIGNIFICANT_THRESHOLD = 1e-6  
PLOT_ORDERS = [8, 16, DEFAULT_ORDER]  # low order is step-sensitive

def sweep_interval_lengths(t_values, order=DEFAULT_ORDER):
    """
    Runs one MMPA call per T, evaluated at t = T only, against the Bateman
    reference. Returns the max relative error and lambda*T at each T.
    """
    a = build_linear_chain_matrix(lambdas)
    lambda1 = lambdas[0]

    max_rel_errors = np.zeros(t_values.size)
    for k, big_t in enumerate(t_values):
        ref = bateman_linear_chain(lambdas, big_t, n0_parent=n0[0])[:, 0]
        got = mmpa_linear_chain(a, big_t, n0, order=order)[:, 0]

        significant = np.abs(ref) >= SIGNIFICANT_THRESHOLD
        assert significant.any(), f"mask is empty at T={big_t:.3e}; nothing was checked"
        rel_err = np.abs(got[significant] - ref[significant]) / np.abs(ref[significant])
        max_rel_errors[k] = rel_err.max()

    return max_rel_errors, lambda1 * t_values
