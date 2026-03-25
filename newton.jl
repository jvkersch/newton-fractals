using CairoMakie
using Polynomials


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
    grid = [newton(p, dp, x + y*im, maxiter=maxiter, tol=tol) for x in xs, y in ys]
    return(; grid, xs, ys)
end

function classify_root(z, known_roots; tol=1e-6)
    which = argmin(abs.(known_roots .- z))
    abs(known_roots[which] - z) < tol ? which : 0
end

function classify_grid(grid, known_roots; tol=1e-6)
    zs = first.(grid)
    classify_root.(zs, Ref(known_roots); tol=tol)
end

function render_basins(xs, ys, root_indices, iterations; filename="newton.png")
    fig, ax, hm = heatmap(xs, ys, root_indices)
    save(filename, fig)
end


if abspath(PROGRAM_FILE) == @__FILE__
    p = Polynomial([-1, 0, 0, 1])
    println("p(z) = ", p)
    println("roots: ", roots(p))

    grid, xs, ys = basin_grid(p, (-2, 2), (-2, 2), 800, 800)
    indices = classify_grid(grid, roots(p))
    iterations = last.(grid)
    render_basins(xs, ys, indices, iterations)
end
