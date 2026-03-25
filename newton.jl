using CairoMakie
using Polynomials

p = Polynomial([-1, 0, 0, 1])
println("p(z) = ", p)
println("roots: ", roots(p))

function newton(p, dp, z0; maxiter=100, tol=1e-10)
    i = 0
    z = z0
    while i < maxiter
        i += 1
        z = z - p(z) / dp(z)
        abs(p(z)) < tol && break
    end
    return (z, i)
end

function basin_grid(p, xrange, yrange, nx, ny; maxiter=100, tol=1e-10)
    dp = derivative(p)
    xs = range(xrange[1], xrange[2], length=nx)
    ys = range(yrange[1], yrange[2], length=ny)
    [newton(p, dp, x + y*im, maxiter=maxiter, tol=tol) for x in xs, y in ys]
end
