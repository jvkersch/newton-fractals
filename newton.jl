# run with "julia --project=. -t auto newton.jl"

using CairoMakie
using Colors
using LaTeXStrings
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

function build_color_matrix(nroots, indices, iterations)
    base_colors = vcat(
        [RGBf(0, 0, 0)], [convert(RGB, HSV{Float64}(i/nroots * 360, 0.8, 0.9)) for i in 1:nroots]
    )
    log_iters = log.(iterations)
    scale = (log_iters .- minimum(log_iters)) ./ (maximum(log_iters) - minimum(log_iters))
    colors = base_colors[indices .+ 1] .* (1 .- scale)
    return colors
end

function render_basins(xs, ys, colors, p; filename="newton.png")
    fig = Figure(backgroundcolor=:black, size=(842, 595))
    ax = Axis(fig[1, 1])

    rs = roots(p)
    image!((xs[begin], xs[end]), (ys[begin], ys[end]), colors; interpolate=false)
    scatter!(real.(rs), imag.(rs), strokecolor=:black, strokewidth=1, color=:transparent)
    Label(fig[0, 1], poly_to_latex(p), color=:white, fontsize=24)

    colsize!(fig.layout, 1, Aspect(1, 1))
    hidedecorations!(ax)
    hidespines!(ax)

    save(filename, fig; px_per_unit=4)
    return filename
end

function poly_to_latex(coeffs::AbstractVector{<:Integer})
    terms = String[]
    n = length(coeffs)

    for i in n:-1:1
        a = coeffs[i]
        a == 0 && continue

        p = i - 1
        abs_a = abs(a)

        varpart =
            p == 0 ? "" :
            p == 1 ? "z" :
            "z^{$p}"

        coeffpart =
            p == 0 ? string(abs_a) :
            abs_a == 1 ? "" :
            string(abs_a)

        term = coeffpart * varpart

        if isempty(terms)
            push!(terms, a < 0 ? "-" * term : term)
        else
            push!(terms, a < 0 ? " - " * term : " + " * term)
        end
    end

    LaTeXString(isempty(terms) ? "0" : join(terms))
end

poly_to_latex(p::Polynomial{<:Integer}) = poly_to_latex(coeffs(p))

function poly_to_latex(coeffs::AbstractVector{<:Complex{<:Integer}})
    terms = String[]
    n = length(coeffs)

    for i in n:-1:1
        c = coeffs[i]
        c == 0 && continue

        p = i - 1
        a, b = real(c), imag(c)

        varpart =
            p == 0 ? "" :
            p == 1 ? "z" :
            "z^{$p}"

        if b == 0
            negative = a < 0
            abs_a = abs(a)
            coeffpart =
                p == 0 ? string(abs_a) :
                abs_a == 1 ? "" :
                string(abs_a)
        elseif a == 0
            negative = b < 0
            abs_b = abs(b)
            coeffpart = abs_b == 1 ? "i" : "$(abs_b)i"
        else
            negative = false
            sign_str = b < 0 ? "-" : "+"
            coeffpart = "($a $sign_str $(abs(b))i)"
        end

        term = coeffpart * varpart

        if isempty(terms)
            push!(terms, negative ? "-" * term : term)
        else
            push!(terms, negative ? " - " * term : " + " * term)
        end
    end

    LaTeXString(isempty(terms) ? "0" : join(terms))
end

poly_to_latex(p::Polynomial{<:Complex{<:Integer}}) = poly_to_latex(coeffs(p))


if abspath(PROGRAM_FILE) == @__FILE__

    const maxiter = 100
    const xrange = (-2, 2)
    const yrange = (-2, 2)
    const nx = 4000
    const ny = 4000
    const p = Polynomial([-1, 0, 0, 1])
    # const p = Polynomial([-1, 0, 0, 0, 0, 1])
    # const p = Polynomial([-6+2im, 4+4im, 1-7im, 1+3im, -4-2im,4])
    # const p = Polynomial([1, 1, 4, 4])

    println("p(z) = ", p)
    println("roots: ", roots(p))

    grid, xs, ys = basin_grid(p, xrange, yrange, nx, ny; maxiter=maxiter)
    indices = classify_grid(grid, roots(p))
    iterations = last.(grid)

    println("min/max iterations: ", extrema(iterations))

    colors = build_color_matrix(degree(p), indices, iterations)
    fname = render_basins(xs, ys, colors, p)
    println("Figure saved to ", fname)
end
