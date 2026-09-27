"""
plot_error_vs_order.py

Plots MMPA max relative error vs Bateman against published order for both
test chains, using the same error measure as tests/test_mmpa_order_comparison.py.
Writes plots/error_vs_order.png.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from data.four_isotope_chain import lambdas as l4, n0 as n0_4, t_grid as t_4
from data.gd157_chain import lambdas as l_gd, n0 as n0_gd, t_grid as t_gd
from src.chain import build_linear_chain_matrix
from src.mmpa import MMPA_COEFFS
from tests.test_mmpa_order_comparison import _max_rel_err, ORDER_32_GATE

CHAINS = {
    "Four-isotope chain": (l4, n0_4, t_4),
    "Gd-157 absorber chain": (l_gd, n0_gd, t_gd),
}


def main():
    orders = sorted(MMPA_COEFFS)
    fig, ax = plt.subplots(figsize=(7, 4.5))

    for (label, (lambdas, n0, t)), marker in zip(CHAINS.items(), ["o", "s"]):
        a_matrix = build_linear_chain_matrix(lambdas)
        errors = [_max_rel_err(a_matrix, lambdas, n0, t, order=order) for order in orders]
        ax.semilogy(orders, errors, marker=marker, label=label)

    ax.axhline(ORDER_32_GATE, color="gray", linestyle="--", label="Order-32 gate (1e-8)")
    ax.set_xlabel("MMPA order")
    ax.set_ylabel("Max relative error vs Bateman")
    ax.set_xticks(orders)
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()

    out = Path(__file__).resolve().parent.parent / "plots" / "error_vs_order.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)


if __name__ == "__main__":
    main()
