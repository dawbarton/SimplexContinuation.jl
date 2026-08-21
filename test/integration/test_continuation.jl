using SimplicialContinuation
using Test

@testset "Continuation" begin
    @testset "Linear zero set (diagonal)" begin
        # f: R^2 → R^1, zero set is x1 = x2
        f(x, p) = [x[1] - x[2]]
        y0 = [0.5, 0.5]
        path = continuation(f, nothing, y0; grain = 0.5, maxsteps = 20)

        points = collect(Iterators.take(path, 10))
        @test length(points) == 10

        # Every yielded point should satisfy f ≈ 0
        for pt in points
            @test abs(f(pt, nothing)[1]) < 1e-10
        end

        # Points should be on the line x1 = x2
        for pt in points
            @test abs(pt[1] - pt[2]) < 1e-10
        end
    end

    @testset "Circle in 2D" begin
        # f: R^2 → R^1, zero set is unit circle
        f(x, p) = [x[1]^2 + x[2]^2 - 1.0]
        y0 = [1.0, 0.0]  # known point on zero set
        path = continuation(f, nothing, y0; grain = 0.2, maxsteps = 200)

        points = collect(Iterators.take(path, 50))
        @test length(points) == 50

        # Each point should lie approximately on the unit circle
        for pt in points
            r = sqrt(pt[1]^2 + pt[2]^2)
            # Linear approximation error scales with grain^2; grain=0.2 → ~0.04 tolerance
            @test abs(r - 1.0) < 0.1
        end
    end

    @testset "3D helix: line in 3D" begin
        # f: R^3 → R^2, zero set is the z-axis (x1=0, x2=0)
        f(x, p) = [x[1], x[2]]
        y0 = [0.0, 0.0, 1.0]
        path = continuation(f, nothing, y0; grain = 0.5, maxsteps = 20)

        points = collect(Iterators.take(path, 10))
        @test length(points) == 10

        for pt in points
            res = f(pt, nothing)
            @test abs(res[1]) < 1e-10
            @test abs(res[2]) < 1e-10
        end
    end

    @testset "ContinuationPath properties" begin
        f(x, p) = [x[1] - x[2]]
        path = continuation(f, nothing, [0.5, 0.5]; grain = 0.5, maxsteps = 5)

        @test Base.IteratorSize(typeof(path)) == Base.SizeUnknown()
        @test Base.eltype(typeof(path)) == Vector{Float64}

        # Respects maxsteps
        points = collect(path)
        @test length(points) ≤ 5
    end

    @testset "Vector grain (per-dimension scaling)" begin
        # Ellipse: (x/2)^2 + y^2 = 1, i.e. x ∈ [-2,2], y ∈ [-1,1].
        # Using a scalar grain sized for y would be too coarse for x; a vector
        # grain lets each dimension be scaled independently.
        f(x, _) = [(x[1]/2)^2 + x[2]^2 - 1.0]
        y0 = [2.0, 0.0]
        path = continuation(f, nothing, y0; grain = (0.4, 0.2), maxsteps = 200)

        points = collect(Iterators.take(path, 20))
        @test length(points) == 20
        for pt in points
            @test abs((pt[1]/2)^2 + pt[2]^2 - 1.0) < 0.1
        end

        # Tuple and Vector grain should both be accepted
        path_v = continuation(f, nothing, y0; grain = [0.4, 0.2], maxsteps = 10)
        @test length(collect(path_v)) ≤ 10
    end

    @testset "Parameters are forwarded" begin
        # f: R^2 → R^1 with parameter p (radius)
        f(x, p) = [x[1]^2 + x[2]^2 - p]
        y0 = [2.0, 0.0]   # on circle of radius 4
        path = continuation(f, 4.0, y0; grain = 0.3, maxsteps = 100)

        points = collect(Iterators.take(path, 20))
        for pt in points
            r2 = pt[1]^2 + pt[2]^2
            @test abs(r2 - 4.0) < 0.2
        end
    end
end
