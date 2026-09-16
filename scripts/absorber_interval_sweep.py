"""
absorber_interval_sweep.py

Writes the CSV and figure for the Gd-157 absorber interval-length
sweep.
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from data.gd157_chain import CAPTURE_RATE, interval_sweep_times
from src.absorber_sweep import PLOT_ORDERS, TOLERANCE, sweep_interval_lengths
from src.mmpa import DEFAULT_ORDER

OUTPUT_DIR = ROOT / "results"


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
    literature effect directly.
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

    top_ax = ax.secondary_xaxis(
        "top",
        functions=(lambda t: t * CAPTURE_RATE, lambda spt: spt / CAPTURE_RATE),
    )
    top_ax.set_xlabel(r"$\sigma\phi\, T$")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    t_values = interval_sweep_times()

    errors_by_order = {}
    for order in PLOT_ORDERS:
        errors, _ = sweep_interval_lengths(t_values, order=order)
        errors_by_order[order] = errors

    max_rel_errors = errors_by_order[DEFAULT_ORDER]
    sigma_phi_t = CAPTURE_RATE * t_values

    csv_path = OUTPUT_DIR / "absorber_interval_sweep.csv"
    png_path = OUTPUT_DIR / "absorber_interval_sweep.png"
    write_csv(t_values, errors_by_order, csv_path)
    write_plot(t_values, errors_by_order, png_path)

    print(f"wrote {csv_path}")
    print(f"wrote {png_path}")
    print()
    print(f"{'T (s)':>12}  {'sigma*phi*T':>12}  {'max rel error (order ' + str(DEFAULT_ORDER) + ')':>24}")
    for k in range(0, t_values.size, 5):
        print(f"{t_values[k]:12.4e}  {sigma_phi_t[k]:12.3e}  {max_rel_errors[k]:24.3e}")

    worst_k = int(np.argmax(max_rel_errors))
    print()
    print(
        f"worst T: {t_values[worst_k]:.4e} s (sigma*phi*T = {sigma_phi_t[worst_k]:.3e}), "
        f"max rel error = {max_rel_errors[worst_k]:.3e}"
    )


if __name__ == "__main__":
    main()
