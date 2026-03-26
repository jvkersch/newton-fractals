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
    # timing single threaded: 339 ms (median) for 1000 x 1000 grid
    # timing with 4 threads: 91 ms (median)
    dp = derivative(p)
    xs = range(xrange[1], xrange[2], length=nx)
    ys = range(yrange[1], yrange[2], length=ny)
    grid = Matrix{Tuple{ComplexF64, Int}}(undef, nx, ny)

    Threads.@threads for j in 1:ny
        for i in 1:nx
            grid[i, j] = newton(p, dp, xs[i] + ys[j]*im, maxiter=maxiter, tol=tol)
        end
    end

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

function build_color_matrix(indices, iterations, maxiter)
    base_colors = [
        RGBf(0.00, 0.00, 0.00),
        RGBf(0.90, 0.20, 0.20),
        RGBf(0.20, 0.70, 0.30),
        RGBf(0.20, 0.40, 0.90)
    ]
    log_iters = log.(iterations)
    scale = (log_iters .- minimum(log_iters)) ./ (maximum(log_iters) - minimum(log_iters))
    colors = base_colors[indices .+ 1] .* (1 .- scale)
    return colors
end

function render_basins(xs, ys, colors; filename="newton.png")
    fig, ax, hm = image((xs[begin], xs[end]), (ys[begin], ys[end]), colors; interpolate=false)
    save(filename, fig)
    return filename
end

if abspath(PROGRAM_FILE) == @__FILE__

    const maxiter = 100
    const xrange = (-2, 2)
    const yrange = (-2, 2)
    const nx = 1000
    const ny = 1000
    const p = Polynomial([-1, 0, 0, 1])

    println("p(z) = ", p)
    println("roots: ", roots(p))

    grid, xs, ys = basin_grid(p, xrange, yrange, nx, ny; maxiter=maxiter)
    indices = classify_grid(grid, roots(p))
    iterations = last.(grid)

    println("min/max iterations: ", extrema(iterations))

    colors = build_color_matrix(indices, iterations, maxiter)
    fname = render_basins(xs, ys, colors)
    println("Figure saved to ", fname)
end
