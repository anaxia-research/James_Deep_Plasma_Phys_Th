# Exact tearing threshold of the periodic current sheet in Hall MHD

Machine checks for the paper by J. M. Lockwood and D. Bhattacharjee.

- `lean/Tearing.lean`: 10 algebraic lemmas (core Lean 4, no Mathlib). Check with `lean lean/Tearing.lean`.
- `julia/`: numerical checks of the integrals, the self-adjoint pencil and the growth-rate threshold (standard library only). Run `julia julia/verify_integrals.jl`, `julia julia/verify_pencil.jl`, `julia julia/verify_threshold.jl`.
- `python/nonlinear_growth.py`: independent check that integrates the full nonlinear four-field equations (numpy only) from small noise and fits the growth rate; it matches the eigenvalues to four digits for k=0.5 and decays for k=1.2.
