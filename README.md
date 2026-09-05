# MMPA-Decay-Solver
Mini-Max Polynomial Approximation (MMPA) solver for matrix exponentials in nuclear decay chains.

**Status:** Early development. Bateman analytical reference is implemented and tested (`src/bateman.py`, `tests/test_bateman.py`). The core MMPA matrix-exponential routine is not started.

**Author:** Jullian J. Arredondo (ArredondoJullianJoseph@gmail.com)
**Research supervisor:** Dr. Timo Budarz (tbudarz@liberty.edu)

## Project scope
- Decay chain verification (4 isotopes) against the Bateman analytical solution
- Absorber case (gadolinium-157) testing multi-time-point evaluation from a single matrix solve
- Performance benchmarking against SciPy Radau
- Speed-accuracy trade-off analysis
- Target success criterion: max relative error below 1e-4 for both cases

## Current state
Bateman analytical reference is implemented and tested. The core MMPA matrix-exponential routine has not been written yet. `src/`, `tests/`, and `data/` hold real code now; `benchmarks/`, `analysis/`, `docs/`, and `slides/` are still placeholders.

## Planned repository layout
- `src/` — solver implementation
- `tests/` — verification against analytical solutions
- `benchmarks/` — performance comparison against SciPy Radau
- `analysis/` — error and convergence analysis
- `data/` — decay constants and reference values
- `docs/` — documentation
- `slides/` — presentation materials

## Next steps
- Push the core MMPA matrix-exponential routine
- Extend to the gadolinium-157 absorber case
- Compare against SciPy Radau

## References

Bateman H. 1910. Solution of a system of differential equations occurring in the theory of radio-active transformations. Proc Camb Philos Soc. 15:423-427.
- Closed-form solution used in `src/bateman.py` as the analytical reference the MMPA solver will be checked against.
