"""
test_gd157_interval_sweep.py

Enforces the absorber-case success criterion of a relative error sbelow 1e-4 at t = T as T is swept across three
decades, for order 32
"""

import numpy as np

from src.absorber_sweep import TOLERANCE, sweep_interval_lengths
from src.mmpa import DEFAULT_ORDER

def test_order_32_holds_tolerance_across_interval_sweep():
    t_values = np.logspace(-3, 0, 31)
    max_rel_errors, lambda_t = sweep_interval_lengths(t_values, order=DEFAULT_ORDER)

    worst_k = int(np.argmax(max_rel_errors))
    print(
        f"\norder {DEFAULT_ORDER} vs Bateman, Gd-157 interval sweep: "
        f"worst T = {t_values[worst_k]:.3e} s (lambda*T = {lambda_t[worst_k]:.3e}), "
        f"max relative error = {max_rel_errors[worst_k]:.3e}"
    )
    assert max_rel_errors[worst_k] < TOLERANCE
