"""
test_gd157_interval_sweep.py

Enforces the absorber-case success criterion of a relative error below 1e-4
at t = T as T is swept so that sigma*phi*T covers three decades, for order 32.
"""

from data.gd157_chain import interval_sweep_times
from src.absorber_sweep import TOLERANCE, sweep_interval_lengths
from src.mmpa import DEFAULT_ORDER


def test_order_32_holds_tolerance_across_interval_sweep():
    t_values = interval_sweep_times()
    max_rel_errors, sigma_phi_t = sweep_interval_lengths(t_values, order=DEFAULT_ORDER)

    worst_k = int(max_rel_errors.argmax())
    print(
        f"\norder {DEFAULT_ORDER} vs Bateman, Gd-157 interval sweep: "
        f"worst T = {t_values[worst_k]:.3e} s "
        f"(sigma*phi*T = {sigma_phi_t[worst_k]:.3e}), "
        f"max relative error = {max_rel_errors[worst_k]:.3e}"
    )
    assert max_rel_errors[worst_k] < TOLERANCE
