"""
radau.py
Reference solver for the decay chain ODE using SciPy's Radau method
"""
import numpy as np
from scipy.integrate import solve_ivp

DEFAULT_RTOL = 1e-10
DEFAULT_ATOL = 1e-12


def radau_linear_chain(a_matrix, t, n0, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """
    Integrates dN/dt = A N with SciPy's Radau method and returns N(t) for every time in t
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    n0 = np.asarray(n0, dtype=float)
    n = a_matrix.shape[0]

    if t.size == 1:
        if t[0] != 0.0:
            raise ValueError("a single-element t must be 0.0; Radau needs a nonzero span to integrate over")
        return n0.reshape(n, 1).copy()

    if np.any(np.diff(t) <= 0.0):
        raise ValueError("t must be strictly increasing") 

    def rhs(_t, y):
        return a_matrix @ y
      
    def jac(_t, _y):
        return a_matrix 

    sol = solve_ivp(
        rhs, t_span=(t[0], t[-1]), y0=n0, method="Radau",
        t_eval=t, rtol=rtol, atol=atol, jac=jac,
    )
    if not sol.success:
        raise RuntimeError(f"Radau integration failed: {sol.message}")
    return sol.y
