using SimplicialContinuation
using Printf

# Trace the intersection of a sphere and a cylinder in R³.
#
#   sphere:   x₁² + x₂² + x₃² = 2
#   cylinder: x₁² + x₂²       = 1
#
# The intersection is two circles at x₃ = ±1.
# Starting point (1, 0, 1) satisfies both equations.
#
# f : R³ → R² so continuation traces a 1-manifold in 3D.

f(x, _) = [x[1]^2 + x[2]^2 + x[3]^2 - 2.0,
           x[1]^2 + x[2]^2         - 1.0]

y0 = [1.0, 0.0, 1.0]

path = continuation(f, nothing, y0; grain = 0.15, maxsteps = 300)

println("Tracing sphere ∩ cylinder (grain = 0.15)")
println("  point                                   residual norms")
for (i, pt) in enumerate(path)
    res = f(pt, nothing)
    @printf "  [%+.4f, %+.4f, %+.4f]   [%.2e, %.2e]\n" pt[1] pt[2] pt[3] abs(res[1]) abs(res[2])
    i == 30 && break
end
