using CairoMakie
using Polynomials

p = Polynomial([-1, 0, 0, 1])
println("p(z) = ", p)
println("roots: ", roots(p))

function newton(p, z0; maxiter=100, tol=1e-10)
    dp = derivative(p)
    i = 0
    z = z0
    while i < maxiter
        i += 1
        z = z - p(z) / dp(z)
        abs(p(z)) < tol && break
    end
    return (z, i)
end
