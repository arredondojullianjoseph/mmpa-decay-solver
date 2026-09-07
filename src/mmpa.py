"""
mmpa.py
Mini-Max Polynomial Approximation (MMPA) solver 
"""
import numpy as np
from scipy.linalg import lu_factor, lu_solve

# Coefficients transcribed directly from Tables 2-4 of Chiba, Yamamoto, Nagano (2026)
MMPA_COEFFS = {
    16: dict(c=11.55385, a=[
        9.55818297829682e-6, -2.21884266813879e-4, 2.34816021434733e-3,
        -1.48316170495578e-2, 6.20732505511125e-2, -1.79379941734168e-1,
        3.52530860768316e-1, -4.26657273240572e-1, 2.10349191740789e-1,
        1.50152179516416e-1, -2.12029252834189e-1, -3.32538652145767e-2,
        1.22178202678754e-1, 4.35315041682305e-3, -4.55089821578535e-2,
        -1.60727919677468e-4, 8.04903136564369e-3,
    ]),
    32: dict(c=24.13694, a=[
        3.28979735872738e-11, -1.58916844491099e-9, 3.67820454991818e-8,
        -5.42102978872146e-7, 5.70945430126668e-6, -4.56738376635576e-5,
        2.87599439571872e-4, -1.45743552929181e-3, 6.02719933663635e-3,
        -2.04839279890385e-2, 5.72047094494963e-2, -1.30239342084535e-1,
        2.37375193939405e-1, -3.33263630944967e-1, 3.27390685838193e-1,
        -1.52635962819998e-1, -1.14728463172040e-1, 2.41197225880137e-1,
        -8.15420042770902e-2, -1.55339492153136e-1, 1.35971357413276e-1,
        7.33867740449378e-2, -1.15986977045075e-1, -2.79633313486442e-2,
        7.37831013682762e-2, 8.52011363362550e-3, -3.57914361107692e-2,
        -1.94605446533440e-3, 1.24721515959643e-2, 2.92704265189227e-4,
        -2.75700953444178e-3, -2.14229591230169e-5, 2.88145489362067e-4,
    ]),
}

# Order 32 is the default because order 16 does not reach the same accuracy on this solver's test interval.
DEFAULT_ORDER = 32


def expm_mmpa_apply(a_matrix, dt, n0, order=DEFAULT_ORDER):
    """
    Applies exp(a_matrix * dt) to n0 without ever forming the full matrix exponential.
    """
    if order not in MMPA_COEFFS:
        raise ValueError(
            f"no MMPA coefficients transcribed for order {order}; "
            f"available orders are {sorted(MMPA_COEFFS)}" 
        ) # MMPA can't run without knowing the coefficients.
    if dt < 0.0:
        raise ValueError("dt must be non-negative") #A negative timestep has no physical meaning.

  #Pulls out this order's c and coefficient list from the table
    c = MMPA_COEFFS[order]["c"]
    coeffs = MMPA_COEFFS[order]["a"]
  
#Builds x = a*dt - c 
    n = a_matrix.shape[0]
    x_matrix = a_matrix * dt - c * np.eye(n)

    lu_piv = lu_factor(x_matrix)   #factors x
    y = np.asarray(n0, dtype=float).copy() #Copies the initial vector so the loop below can alter y without changing the inital vector.
    result = coeffs[0] * y

  #Each step updates y using the fast solve, then adds its weighted share to the total
    for i in range(1, order + 1):
        y = y + 2.0 * c * lu_solve(lu_piv, y)
        result = result + coeffs[i] * y
    return result


def mmpa_linear_chain(a_matrix, t, n0, order=DEFAULT_ORDER):
    """
    Evaluates exp(A*t_k) n0 for every time in t
    """

    t = np.atleast_1d(np.asarray(t, dtype=float))
    n = a_matrix.shape[0]
    density = np.zeros((n, t.size))

  #Fills in one column of density per time value, then returns the whole grid.
    for k, tk in enumerate(t):
        density[:, k] = expm_mmpa_apply(a_matrix, tk, n0, order=order)
    return density
