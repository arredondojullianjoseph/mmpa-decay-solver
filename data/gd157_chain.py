"""
gd157_chain.py

Two-nuclide chain for the saturable absorber case. Gd-157 removes by
(n,gamma) capture into stable Gd-158. 
"""

import numpy as np

names = ["Gd-157", "Gd-158"]

# Not a literal cross section times flux. It's a clean order-of-magnitude separation for testing.
lambdas = [100.0, 0.0]  # Decay constants (1/s). Gd-158 is stable (0.0).
n0 = [1.0, 0.0]  # Everything starts as Gd-157.
t_grid = np.concatenate(([0.0], np.logspace(-4, 0.0, 49)))
