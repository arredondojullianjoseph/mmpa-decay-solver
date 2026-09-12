"""
absorber_interval_sweep.py

Writes the CSV and figure for the Gd-157 absorber interval-length
sweep. 
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.absorber_sweep import PLOT_ORDERS, TOLERANCE, sweep_interval_lengths
from src.mmpa import DEFAULT_ORDER
from data.gd157_chain import lambdas

OUTPUT_DIR = Path("results")

def write_csv(t_values, errors_by_order, path):
    """
    One error column per order in PLOT_ORDERS, so the CSV alone can
    regenerate the figure.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["T_seconds"] + [f"max_rel_error_order_{o}" for o in PLOT_ORDERS])
        for k, big_t in enumerate(t_values):
            row = [f"{big_t:.6e}"] + [f"{errors_by_order[o][k]:.6e}" for o in PLOT_ORDERS]
            writer.writerow(row)


def write_plot(t_values, errors_by_order, path):
    """
    Overlays one curve per order in PLOT_ORDERS so the figure shows the
    literature effect directly
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    for order in PLOT_ORDERS:
        ax.loglog(t_values, errors_by_order[order], marker="o", label=f"order {order}")
    ax.axhline(TOLERANCE, linestyle="--", color="gray", label=f"tolerance ({TOLERANCE:.0e})")
    ax.set_xlabel("Interval length T (s)")
    ax.set_ylabel("Max relative error at t = T")
    ax.set_title("MMPA accuracy vs. absorber interval length, Gd-157/Gd-158")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)

    lambda1 = lambdas[0]
    top_ax = ax.secondary_xaxis("top", functions=(lambda t: t * lambda1, lambda lt: lt / lambda1))
    top_ax.set_xlabel(r"$\lambda T$")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)

def main():
    t_values = np.logspace(-3, 0, 31)

    errors_by_order = {}
    for order in PLOT_ORDERS:
        errors, _ = sweep_interval_lengths(t_values, order=order)
        errors_by_order[order] = errors

    max_rel_errors = errors_by_order[DEFAULT_ORDER]
    lambda1 = lambdas[0]
    lambda_t = lambda1 * t_values

    csv_path = OUTPUT_DIR / "absorber_interval_sweep.csv"
    png_path = OUTPUT_DIR / "absorber_interval_sweep.png"
    write_csv(t_values, errors_by_order, csv_path)
    write_plot(t_values, errors_by_order, png_path)

    print(f"wrote {csv_path}")
    print(f"wrote {png_path}")
    print()
    print(f"{'T (s)':>12}  {'lambda*T':>10}  {'max rel error (order ' + str(DEFAULT_ORDER) + ')':>24}")
    for k in range(0, t_values.size, 5):
        print(f"{t_values[k]:12.4e}  {lambda_t[k]:10.3e}  {max_rel_errors[k]:24.3e}")

    worst_k = int(np.argmax(max_rel_errors))
    print()
    print(
        f"worst T: {t_values[worst_k]:.4e} s (lambda*T = {lambda_t[worst_k]:.3e}), "
        f"max rel error = {max_rel_errors[worst_k]:.3e}"
    )


if __name__ == "__main__":
    main()
