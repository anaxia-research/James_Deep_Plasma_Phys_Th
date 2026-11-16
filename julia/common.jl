using LinearAlgebra
# Fourier-Galerkin matrices on x in R/2πZ, modes n = -M..M.
function ops(k, M)
    n = collect(-M:M); N = length(n)
    m = n .^ 2 .+ k^2
    Q = Diagonal(-m .+ 0im)
    L = Q + I
    B = zeros(ComplexF64, N, N)          # multiplication by B(x) = -sin x
    for i in 1:N
        i > 1 && (B[i, i-1] += 0.5im)
        i < N && (B[i, i+1] -= 0.5im)
    end
    return n, m, Matrix(Q), Matrix(L), B
end
# Generator of the linearized four-field system (7.1)-(7.4), unknowns (φ, ψ, w, b).
function generator(k, η, h, M)
    n, m, Q, L, B = ops(k, M); N = length(n); Z = zeros(ComplexF64, N, N)
    Qi = inv(Q)
    [Z  Qi*(im*k*B*L)  Z  Z;
     im*k*B  η*Q  Z  im*k*h*B;
     Z  Z  Z  im*k*B;
     Z  -im*h*k*B*L  im*k*B  η*Q]
end
maxgrowth(k, η, h, M) = (ev = eigvals(generator(k, η, h, M)); ev[argmax(real.(ev))])
# Top of the self-adjoint pencil A_h(γ) = (-Q)^{-1/2} M_h(γ) (-Q)^{-1/2}.
function lam(k, η, γ, M; h = 0.0)
    n, m, Q, L, B = ops(k, M)
    X = Diagonal((γ .+ η .* m) ./ (1 .- m))
    Y = Matrix{ComplexF64}(X)
    if h != 0
        P = γ * I + η * Diagonal(m) + (k^2 / γ) * (B * B)
        Y -= h^2 * k^2 * (B * (P \ B))
    end
    Mh = k^2 * B * (Y \ B)
    D = Diagonal(1 ./ sqrt.(m))
    A = Hermitian((D * Mh * D + (D * Mh * D)') / 2)
    return eigmax(A)
end
function bisect(f, a, b; it = 60)
    fa = f(a)
    for _ in 1:it
        c = (a + b) / 2; fc = f(c)
        if (fc > 0) == (fa > 0); a, fa = c, fc; else; b = c; end
    end
    (a + b) / 2
end
# Largest root of f on [lo, hi]: scan a geometric grid downward from hi for the
# first sign change from negative (right) to positive (left), then bisect.
function largest_root(f, lo, hi; npts = 200)
    g = exp.(range(log(hi), log(lo), length = npts))
    fr = f(g[1])
    for i in 2:npts
        fl = f(g[i])
        if fl > 0 && fr < 0
            return bisect(f, g[i], g[i-1])
        end
        fr = fl
    end
    error("no root")
end
const PASS = Ref(0)
check(name, cond) = (cond || error("FAIL: $name"); PASS[] += 1; println("PASS  ", name))
