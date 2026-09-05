"""
four_isotope_chain.py

A made-up 4-nuclide decay chain (N1 -> N2 -> N3 -> N4, N4 stable) for
testing the Bateman solver. The decay constants are round numbers, not
real half-lives so we can check the math works before we plug
in the actual U-238 chain.
"""

import numpy as np

names = ["N1", "N2", "N3", "N4"]
lambdas = [1.0, 0.5, 0.2, 0.0]  # Decay constants (1/s). N4 is stable (0.0).
n0 = [1.0, 0.0, 0.0, 0.0]  # Everything starts as N1.

# t=0 plus a log-spaced grid out to 50 s. Log spacing because N1 decays fast (its transient is over in the first couple seconds) 
t_grid = np.concatenate(([0.0], np.logspace(-2, np.log10(50.0), 49)))
