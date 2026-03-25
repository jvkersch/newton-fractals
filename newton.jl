using CairoMakie
using Polynomials

p = Polynomial([-1, 0, 0, 1])
println("p(z) = ", p)
println("roots: ", roots(p))
