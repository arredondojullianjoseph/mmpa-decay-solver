"""
gd157_chain.py

Two-nuclide burnup step for the saturable absorber case. 
"""

import numpy as np

names = ["Gd-157", "Gd-158"]

# Order-of-magnitude PWR-like values, not a lattice calculation.
SIGMA_N_GAMMA_CM2 = 2.54e-19  
PHI_N_CM2_S = 1.0e14
CAPTURE_RATE = SIGMA_N_GAMMA_CM2 * PHI_N_CM2_S 
CHARACTERISTIC_TIME = 1.0 / CAPTURE_RATE

# Effective removal rates for the linear-chain helpers. 
lambdas = [CAPTURE_RATE, 0.0]
n0 = [1.0, 0.0]  # Everything starts as Gd-157.

t_grid = np.concatenate(([0.0], CHARACTERISTIC_TIME * np.logspace(-2, 2.0, 49)))

def interval_sweep_times(n=31):
  
    """Endpoint times with sigma*phi*T spanning [0.1, 100]."""
  
    return CHARACTERISTIC_TIME * np.logspace(-1, 2.0, n)
