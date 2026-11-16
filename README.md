# Exact tearing threshold of the periodic current sheet in Hall MHD

Machine checks for the paper by J. M. Lockwood and D. Bhattacharjee.

- `lean/Tearing.lean`: 10 algebraic lemmas (core Lean 4, no Mathlib). Check with `lean lean/Tearing.lean`.
- `julia/`: numerical checks of the integrals, the self-adjoint pencil and the growth-rate threshold (standard library only). Run `julia julia/verify_integrals.jl`, `julia julia/verify_pencil.jl`, `julia julia/verify_threshold.jl`.
