# Section 5-7: the real root of λ_max(γ) = γ equals the top eigenvalue of the
# generator; the bounds γ ≤ sqrt((1-k^2)(1+h^2k^2)) and γ ≥ b η^{3/5}.
include(joinpath(@__DIR__, "common.jl"))
k = 0.375
# For large h the pencil has a pole just below the root (Lemma 5.4), so the
# scan needs a fine grid there.
for (η, h, M, np) in ((1e-3, 0.0, 128, 200), (1e-3, 0.5, 128, 200), (1e-3, 3.0, 128, 200), (2e-5, 0.1, 160, 200), (1e-2, 10.0, 96, 6000))
    γp = largest_root(γ -> lam(k, η, γ, M; h = h) - γ, 1e-5, 1.01 * sqrt((1 - k^2) * (1 + h^2 * k^2)); npts = np)
    γf = real(maxgrowth(k, η, h, M))
    check("pencil root = generator eigenvalue (η=$η, h=$h): $γp vs $γf", abs(γp - γf) < 1e-8 * max(1, γf) + 1e-10)
    check("upper bound γ^2 < (1-k^2)(1+h^2k^2) (η=$η, h=$h)", γp^2 < (1 - k^2) * (1 + h^2 * k^2))
end
b = (k^2 * (1 - k^2)^4 / 512)^(1 / 5); a = (1 - k^2) / (4b)
η0 = min(((1 - sqrt(3) / 2) / a)^(5 / 2), (b / (2k^2))^(5 / 2))
println("b = $b, a = $a, η0 = $η0")
for (η, M) in ((3e-3, 96), (1e-3, 128), (3e-4, 192), (1e-4, 256), (3e-5, 384))
    η > η0 && continue
    γ1 = b * η^(3 / 5); ε = a * η^(2 / 5)
    # test function χ = -sin x/(sin^2 x + ε^2) evaluated by quadrature
    x = range(0, 2π, length = 400001)[1:end-1]; s = sin.(x); c = cos.(x)
    χ = -s ./ (s .^ 2 .+ ε^2); f = s .^ 2 ./ (s .^ 2 .+ ε^2)
    dχ = -c .* (ε^2 .- s .^ 2) ./ (s .^ 2 .+ ε^2) .^ 2
    Nχ = sum(dχ .^ 2 .+ k^2 .* χ .^ 2) / length(x)          # N(χ)/(2π)
    f0 = sum(f) / length(x); var = sum((f .- f0) .^ 2) / length(x)
    check("N(χ)/2π ≤ ε^-3 (η=$η)", Nχ <= ε^-3)
    check("Var f ≤ ε (η=$η)", var <= ε)
    qlow = k^2 * ε^3 * ((1 - k^2) * (1 - ε)^2 / (γ1 + η * k^2) - ε / η)
    check("q(γ1) ≥ 2γ1 from the explicit bound (η=$η)", qlow >= 2γ1 * (1 - 1e-12))
    γp = largest_root(γ -> lam(k, η, γ, M) - γ, 1e-6, 2.0; npts = 60)
    check("γ* = $(round(γp, sigdigits=5)) ≥ b η^{3/5} = $(round(γ1, sigdigits=5))", γp >= γ1)
end
println("ALL PASS  $(PASS[]) checks")
