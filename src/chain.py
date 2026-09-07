"""
chain.py
Builds the decay matrix A
"""
import numpy as np

def build_linear_chain_matrix(lambdas):
    lambdas = np.asarray(lambdas, dtype=float)
    n = lambdas.size
    if n == 0:
        raise ValueError("lambdas must contain at least one decay constant")
    if np.any(lambdas < 0.0):
        raise ValueError("decay constants must be non-negative")
    a = np.zeros((n, n))
    for i in range(n):
        a[i, i] = -lambdas[i] #removal of nuclide i by its own decay
        if i + 1 < n:
            a[i + 1, i] = lambdas[i] #production of nuclide i+1 from i
    return a
