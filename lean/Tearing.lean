/-!
Finite algebraic facts used in "Exact Tearing Threshold of the Periodic
Current Sheet in Resistive Hall Magnetohydrodynamics".
Lean 4 core only (no Mathlib). Squared wavenumbers are written k^2 = p/q
with natural numbers p, q > 0, and inequalities are cleared of denominators.
-/
namespace Tearing

/-- `n^2 + k^2 - 1 ≥ k^2 n^2` for every integer `n ≠ 0` when `k^2 = p/q ≤ 1`,
cleared of denominators: `p n^2 + q ≤ q n^2 + p`. -/
theorem coef_lower (p q n : Nat) (hpq : p ≤ q) (hn : 1 ≤ n) :
    p * (n * n) + q ≤ q * (n * n) + p := by
  have hnn : 1 ≤ n * n := Nat.mul_le_mul hn hn
  have hsplit : q * (n * n) = p * (n * n) + (q - p) * (n * n) := by
    rw [Nat.sub_mul]; have := Nat.mul_le_mul_right (n * n) hpq; omega
  have hmono : (q - p) * 1 ≤ (q - p) * (n * n) := Nat.mul_le_mul_left _ hnn
  omega

/-- For `k^2 = p/q ≥ 1` every coefficient `n^2 + k^2 - 1` of the
free-energy form is nonnegative: `q ≤ q n^2 + p`. -/
theorem free_energy_coef_nonneg (p q n : Nat) (hqp : q ≤ p) :
    q ≤ q * (n * n) + p := by omega

/-- For `m = n^2 + k^2 ≥ 1` the dissipation coefficient `m (m - 1)` is
nonnegative; with `m = M/q`, `M ≥ q`: `M (M - q) ≥ 0` in the integers. -/
theorem dissipation_coef_nonneg (M q : Int) (h : q ≤ M) (hq : 0 ≤ q) :
    0 ≤ M * (M - q) := by
  have h1 : 0 ≤ M := by omega
  have h2 : 0 ≤ M - q := by omega
  exact Int.mul_nonneg h1 h2

/-- Sign of the single positive Fourier weight: for `k^2 = p/q < 1`,
`1 - k^2 > 0`, i.e. `p < q` gives `0 < q - p`. -/
theorem mean_weight_pos (p q : Nat) (h : p < q) : 0 < q - p := by omega

/-- The weights `(m - 1)/(γ + η m)` for `n ≠ 0` never exceed `1/η`:
cleared of denominators, `η (m - 1) ≤ γ + η m` for `γ, η ≥ 0`. -/
theorem weight_le_inv_eta (γ η m : Int) (hγ : 0 ≤ γ) (hη : 0 ≤ η) :
    η * (m - 1) ≤ γ + η * m := by
  rw [Int.mul_sub, Int.mul_one]; omega

/-- Hall-sign similarity: replacing `(w, b, h)` by `(-w, -b, -h)` leaves the
four linear equations invariant. Each Hall or axial coupling term is a
product of two sign-flipped factors or one flipped factor on a flipped
left side; over a commutative ring this is `(-h) * (-b) = h * b` and
`-(h * x) = (-h) * x`. -/
theorem hall_sign_psi (h b : Int) : (-h) * (-b) = h * b := by
  rw [Int.neg_mul_neg]
theorem hall_sign_b (h x : Int) : -(h * x) = (-h) * x := by
  rw [Int.neg_mul]

/-- Exponent bookkeeping in the lower bound: with `ε = a η^{2/5}` and
`γ₁ = b η^{3/5}`, the powers of `η^{1/5}` satisfy `3·2 = 6`, `6 - 3 = 3`
and `2 - 5 = -3`. -/
theorem exponents : (3 * 2 : Int) = 6 ∧ (6 : Int) - 3 = 3 ∧ (2 : Int) - 5 = -3 := by
  decide

/-- Constant bookkeeping: if `512 b^5 = K` and `256 b^4 a^4 = K'`
with `K' = K` then `a^4 · 256 b^4 = 2 b · 256 b^4`, i.e. `k^2 a^4 = 2 b`
after cancelling `256 b^4`. Checked for the integer identity
`512 = 2 * 256`. -/
theorem const_512 : (512 : Nat) = 2 * 256 := by decide

/-- Box corollary arithmetic: if the along-sheet period satisfies
`L ≤ 2π`, every admissible `k_j = 2πj/L` with `j ≥ 1` has `k_j ≥ 1`;
cleared of `π`: `L ≤ 2π` and `j ≥ 1` give `L ≤ 2π j`, here with
`2π` replaced by any positive natural `P`. -/
theorem box_modes (L P j : Nat) (hL : L ≤ P) (hj : 1 ≤ j) : L ≤ P * j := by
  have : P * 1 ≤ P * j := Nat.mul_le_mul_left _ hj
  omega

end Tearing

#print axioms Tearing.coef_lower
#print axioms Tearing.free_energy_coef_nonneg
#print axioms Tearing.dissipation_coef_nonneg
#print axioms Tearing.mean_weight_pos
#print axioms Tearing.weight_le_inv_eta
#print axioms Tearing.hall_sign_psi
#print axioms Tearing.hall_sign_b
#print axioms Tearing.exponents
#print axioms Tearing.const_512
#print axioms Tearing.box_modes
