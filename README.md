# MMPA-Decay-Solver

Mini-Max Polynomial Approximation (MMPA) solver for matrix exponentials
in nuclear decay chains. Currently in progress.

The analytical Bateman solution for a linear decay chain is implemented
and tested, as are automated tests that cover the parent-only initial condition and atom
conservation. The MMPA matrix-exponential is also implemented and
checked against Bateman on the four-isotope chain. **Next:** SciPy
Radau integration.

**Author:** Jullian J. Arredondo (jjarredondo@liberty.edu)
**Research Supervisor:** Dr. Timo Budarz (tbudarz@liberty.edu)

## Mathematical model

### Governing equation

A linear decay chain
$N_1 \rightarrow N_2 \rightarrow \cdots \rightarrow N_n$
with no branching follows (Bateman 1910)

$$
\frac{\mathrm{d}N_k}{\mathrm{d}t}
= -\lambda_k N_k + \lambda_{k-1} N_{k-1}
$$

for each nuclide $k$, where $\lambda_k$ is its decay constant
($\mathrm{s}^{-1}$) and $\lambda_0 = 0$ (the parent has no source).
A stable end member uses the same equation with that nuclide's
$\lambda_k$ set to $0$.

### Bateman solution

When all atoms start in the parent,
$N_1(0) = n_{0,\mathrm{parent}}$ and $N_k(0) = 0$ for $k \ge 2$,
the Bateman (1910) solution for nuclide $k$ is

$$
N_k(t)
= n_{0,\mathrm{parent}}
  \left( \prod_{j=1}^{k-1} \lambda_j \right)
  \sum_{i=1}^{k}
  \frac{e^{-\lambda_i t}}
       {\displaystyle\prod_{\substack{j=1 \\ j \neq i}}^{k}
        \left( \lambda_j - \lambda_i \right)}.
$$

$\lambda_1,\ldots,\lambda_k$ must be distinct. A repeated value makes
a factor in the product vanish and the formula divides by zero, so
`bateman_linear_chain` rejects repeated decay constants.

### Current implementation

- `bateman_linear_chain` (`src/bateman.py`) evaluates the analytical Bateman solution.
- `expm_mmpa_apply` / `mmpa_linear_chain` (`src/mmpa.py`) evaluates exp(AΔt) via the MMPA method (Kawamoto et al. 2015), using the order-16 and order-32 coefficient tables from Chiba et al. (2026).
- `build_linear_chain_matrix` (`src/chain.py`) assembles the decaymatrix A.

## Repository layout

- `src/` — solver implementation
- `tests/` — verification against analytical solutions
- `data/` — decay constants and reference values

## Next steps

- Extend to the gadolinium-157 absorber case
- Compare MMPA against SciPy Radau

## References

Bateman H. 1910. Solution of a system of differential equations
occurring in the theory of radio-active transformations.
*Proc. Camb. Philos. Soc.* 15:423–427.
- Closed-form solution implemented in `src/bateman.py`.

Kawamoto Y, Chiba G, Tsuji M, Narabayashi T. 2015. Numerical solution
of matrix exponential in burn-up equation using mini-max polynomial
approximation. *Ann. Nucl. Energy* 80:219–224.
- Source of the MMPA method implemented in `src/mmpa.py`.

Chiba G, Yamamoto K, Nagano H. 2026. Revisiting mini-max polynomial
approximation method for nuclear fuel depletion calculation.
*Ann. Nucl. Energy* 227:111948.
- Source of the order-16 and order-32 MMPA coefficient tables used in
  `src/mmpa.py`.
