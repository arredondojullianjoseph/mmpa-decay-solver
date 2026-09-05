"""
bateman.py

Closed-form solution for a linear decay chain (Bateman 1910). This is
the exact analytical solution. Only handles the case where every atom starts in the parent nuclide
(N1) and every daughter starts at zero. 
"""

import numpy as np


def bateman_linear_chain(lambdas, t, n0_parent=1.0):
    """
    Calculates the Bateman equations for a linear radioactive decay
    chain: N1 -> N2 -> ... -> Nn. This assumes a straight chain with no
    branching and requires all decay constants to be unique.
    """
    
    lambdas = np.asarray(lambdas, dtype=float)
    t = np.atleast_1d(np.asarray(t, dtype=float))

    n = lambdas.size
    if n == 0:
        raise ValueError("lambdas must contain at least one decay constant")
    if np.any(lambdas < 0.0):
        
        # A negative decay constant would mean the nuclide grows on its own with no source feeding it.
        raise ValueError("decay constants must be non-negative")
    if np.unique(lambdas).size != n:
        
        # Two equal lambdas would divide by zero below
        raise ValueError(
            "decay constants must be distinct - the formula divides by "
            "(lambda_j - lambda_i), so repeats cause a divide-by-zero"
        )
    if np.any(t < 0.0):
        raise ValueError("times must be non-negative")
    n_times = t.size
    density = np.zeros((n, n_times))
    for k in range(1, n + 1):
        prefactor = np.prod(lambdas[: k - 1]) if k > 1 else 1.0

        chain_slice = lambdas[:k]
        term_sum = np.zeros(n_times)
        for i in range(k):
            denom = 1.0
            for j in range(k):
                if j != i:
                    denom *= chain_slice[j] - chain_slice[i]
            term_sum += np.exp(-chain_slice[i] * t) / denom

        density[k - 1, :] = n0_parent * prefactor * term_sum

    return density
