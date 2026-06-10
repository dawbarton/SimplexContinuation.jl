using SimplexContinuation
using Printf

# Trace the unit circle: f(x) = x₁² + x₂² - 1 = 0
# Starting point (1, 0) is exactly on the curve.

f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
y0 = [1.0, 0.0]

path = continuation(f, nothing, y0; grain = 0.2, maxsteps = 200)

println("Tracing unit circle (grain = 0.2)")
println("  point                          radius")
for (i, pt) in enumerate(path)
    r = sqrt(pt[1]^2 + pt[2]^2)
    @printf "  [%+.4f, %+.4f]             %.6f\n" pt[1] pt[2] r
    i == 30 && break
end
