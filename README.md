# MMPA-Decay-Solver

Mini-Max Polynomial Approximation (MMPA) solver for matrix exponentials in nuclear decay chains. Verifies MMPA against the analytical Bateman solution on a four-isotope decay chain and a two-nuclide Gd-157-style absorber chain, and checks it independently against a SciPy Radau reference solver. Sweeps published MMPA order (4 through 32) on both chains and checks that accuracy improves with order. Automated tests cover all of these checks and print their own measured accuracy.

Scope: verifies MMPA's matrix-exponential accuracy on toy decay chains before scaling to real cross sections and a burnup matrix derived from a low-order neutronics model.

**Author:** Jullian J. Arredondo (jjarredondo@liberty.edu)
(ArredondoJullianJoseph@gmail.com)

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

## Implementation

- **Bateman solution:** `bateman_linear_chain` (`src/bateman.py`) evaluates the analytical solution above.
- **Decay matrix:** `build_linear_chain_matrix` (`src/chain.py`) assembles the decay matrix $A$.
- **MMPA solver:** `expm_mmpa_apply` / `mmpa_linear_chain` (`src/mmpa.py`) evaluates $\exp(A\Delta t)N_0$ via MMPA (Kawamoto et al. 2015), using all ten published coefficient tables from Chiba et al. (2026): orders 4, 6, 8, 10, 12, 16, 20, 24, 28, and 32. Default order 32. One factorization per time point.
- **Radau reference solver:** `radau_linear_chain` (`src/radau.py`) integrates $\mathrm{d}N/\mathrm{d}t = AN$ with SciPy's Radau method, independent of MMPA.
- **Four-isotope chain:** $\lambda = [1.0, 0.5, 0.2, 0.0]$ (`data/four_isotope_chain.py`).
- **Gd-157-style chain:** $\lambda = [100.0, 0.0]$ (`data/gd157_chain.py`). Scaled toy rate, not $\sigma\phi$. Constants are test values by design; real cross sections and flux require the coupled neutronics model, out of scope here.
- **Automated tests:** `tests/` covers all checks below and prints its own measured max relative error when run with `pytest -s`.

## Verification results

### MMPA order 32 vs Bateman

Order 32 is gated at $10^{-8}$ relative error against Bateman on significant inventories ($N \ge 10^{-6}$).

Running `pytest -s tests/test_mmpa_vs_bateman.py tests/test_gd157_vs_bateman.py`:

| Chain | Max relative error |
| --- | --- |
| Four-isotope | $2.657\times10^{-9}$ |
| Gd-157 | $2.624\times10^{-9}$ |

Both agree with Bateman well inside the $10^{-8}$ gate.

### MMPA order sweep

Every published order (4, 6, 8, 10, 12, 16, 20, 24, 28, 32) is checked against Bateman on both chains. A higher order has more polynomial terms, so error must strictly decrease as order goes up; the test fails if any order is not more accurate than the one before it.

Running `pytest -s tests/test_mmpa_order_comparison.py`:

| Order | Four-isotope max relative error | Gd-157 max relative error |
| --- | --- | --- |
| 4 | $5.444\times10^{2}$ | $5.592\times10^{2}$ |
| 6 | $1.739\times10^{2}$ | $7.343\times10^{1}$ |
| 8 | $1.657\times10^{1}$ | $7.472\times10^{0}$ |
| 10 | $3.738\times10^{0}$ | $2.623\times10^{0}$ |
| 12 | $2.960\times10^{-1}$ | $3.016\times10^{-1}$ |
| 16 | $6.915\times10^{-3}$ | $7.036\times10^{-3}$ |
| 20 | $3.184\times10^{-4}$ | $1.233\times10^{-4}$ |
| 24 | $4.430\times10^{-6}$ | $2.316\times10^{-6}$ |
| 28 | $8.671\times10^{-8}$ | $6.698\times10^{-8}$ |
| 32 | $2.657\times10^{-9}$ | $2.624\times10^{-9}$ |

Orders 4 through 12 are far too inaccurate to use (relative errors above $10^{-1}$ for order 12 and worse for lower orders). This matches Chiba et al. (2026): their low-order coefficients target a 1% error on a full PWR pincell burnup matrix. Against this project's $10^{-4}$ success criterion,order 24 is the lowest order that passes on both ($4.430\times10^{-6}$ and $2.316\times10^{-6}$).

### Radau vs Bateman

Radau is gated at the same $10^{-8}$ as MMPA.

Running `pytest -s tests/test_radau_vs_bateman.py`:

| Chain | Max relative error |
| --- | --- |
| Four-isotope | $2.574\times10^{-9}$ |
| Gd-157 | $7.311\times10^{-9}$ |

Both agree with Bateman well inside the $10^{-8}$ gate.

## Limitations

- Toy decay chains, not real cross sections; real cross sections and flux require the coupled neutronics model.
- Multi-time evaluation from one factorization is not implemented; MMPA still factors once per time point.
- No interval-length sweep on the absorber chain yet. The order sweep above varies MMPA order at a fixed time grid; it does not vary the absorber's interval length, which is the actual result the absorber case needs.

## Repository layout

- `src/` — solver implementation
- `tests/` — verification against analytical solutions
- `data/` — decay constants and reference values
- `.github/workflows/tests.yml` — `pytest` on every push

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
- Source of all ten MMPA coefficient tables (orders 4 through 32) used
  in `src/mmpa.py`, and of the order-dependence result in the order
  sweep above.
