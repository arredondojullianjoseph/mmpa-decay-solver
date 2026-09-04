# MMPA-Decay-Solver

Mini-Max Polynomial Approximation (MMPA) solver for matrix exponentials in nuclear decay chains. Verified against the Bateman analytical solution for U-238 → Th-234 → Pa-234m → U-234. Also tested on a gadolinium-157 absorber case for multi-time-point evaluation.

**Scope:** Verifies the MMPA method's multi-time-point evaluation mechanism against known analytical solutions. Does not test neutronics acceleration.

## Project scope

- Decay chain verification (4 isotopes) against the Bateman analytical solution
- Absorber case (gadolinium-157) testing multi-time-point evaluation from a single matrix solve
- Performance benchmarking against SciPy Radau
- Speed-accuracy trade-off analysis
- Target success criterion: max relative error below 1e-4 for both cases

## Current state

Repository scaffolding only. Folders below are placeholders. No code has been written yet.

## Planned repository layout

- `src/` — solver implementation
- `tests/` — verification against analytical solutions
- `benchmarks/` — performance comparison against SciPy Radau
- `analysis/` — error and convergence analysis
- `data/` — decay constants and reference values
- `docs/` — documentation
- `slides/` — presentation materials

## Next steps

- Implement the core MMPA matrix-exponential routine
- Verify against the Bateman solution for the 4-isotope decay chain
- Extend to the gadolinium-157 absorber case
- Compare against SciPy Radau
