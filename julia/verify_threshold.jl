# Theorem A and B numerically: max Re γ > 0 exactly for 0 < k < 1, any h, η.
include(joinpath(@__DIR__, "common.jl"))
for η in (1e-2, 1e-3), h in (0.0, 0.5, 3.0)
    for k in (0.3, 0.5, 0.9, 0.98)
        g = maxgrowth(k, η, h, 72)
        check("unstable k=$k h=$h η=$η (γ=$(round(real(g), sigdigits=5)), real)", real(g) > 1e-5 && abs(imag(g)) < 1e-8)
    end
    for k in (1.0, 1.02, 1.5, 2.0)
        g = maxgrowth(k, η, h, 72)
        check("stable   k=$k h=$h η=$η (max Re=$(round(real(g), sigdigits=3)))", real(g) < 1e-10)
    end
end
# Reproduction of the reference value of the source manuscript (k=0.375, η=2e-5, h=0.1).
g = maxgrowth(0.375, 2e-5, 0.1, 160)
check("reference value 0.005187 at M=160 (got $(real(g)))", abs(real(g) - 0.005187) < 5e-6)
println("ALL PASS  $(PASS[]) checks")
