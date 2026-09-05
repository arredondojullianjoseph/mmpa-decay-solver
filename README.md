# MMPA-Decay-Solver
Mini-Max Polynomial Approximation (MMPA) solver for matrix exponentials in nuclear decay chains.

**Status:** Early development. Solver code is in progress and not yet pushed to this repository. Bateman verification has not been run.

## Project scope
- Decay chain verification (4 isotopes) against the Bateman analytical solution
- Absorber case (gadolinium-157) testing multi-time-point evaluation from a single matrix solve
- Performance benchmarking against SciPy Radau
- Speed-accuracy trade-off analysis
- Target success criterion: max relative error below 1e-4 for both cases

## Current state
Solver code is actively being written and has not yet been pushed to this repository. Folders below are placeholders until that code lands.

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
- Verify against the Bateman solution for the 4-isotope decay chain
- Extend to the gadolinium-157 absorber case
- Compare against SciPy Radau
