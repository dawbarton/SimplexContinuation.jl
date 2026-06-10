using SimplexContinuation

# Trace a circle whose radius is supplied as a parameter.
# f(x, r) = x₁² + x₂² - r = 0
# This demonstrates the SciML-style (f, p, y0) calling convention.

f(x, r) = [x[1]^2 + x[2]^2 - r]

for radius in [1.0, 4.0, 9.0]
    r = Float64(radius)
    y0 = [sqrt(r), 0.0]   # known starting point on the circle
    path = continuation(f, r, y0; grain = 0.3, maxsteps = 100)

    errors = [abs(pt[1]^2 + pt[2]^2 - r) for pt in Iterators.take(path, 20)]
    println("radius = $r  →  max residual over 20 points: $(maximum(errors))")
end
