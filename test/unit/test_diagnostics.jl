using SimplicialContinuation
using LinearAlgebra: norm
using Random
using Test

@testset "Diagnostics" begin
    @testset "check_start_point" begin
        f(x, _) = [x[1]^2 + x[2]^2 - 1.0]

        on_curve = check_start_point(f, nothing, [1.0, 0.0]; grain = 0.2)
        @test on_curve.norm == 0.0
        @test on_curve.ok

        off_curve = check_start_point(f, nothing, [1.5, 0.0]; grain = 0.2)
        @test off_curve.norm ≈ 1.25
        @test !off_curve.ok

        # Vector grain uses the smallest per-dimension value for the threshold.
        near_curve = check_start_point(f, nothing, [1.0 + 1.0e-6, 0.0]; grain = (0.2, 1.0e-8))
        @test near_curve.threshold ≈ 1.0e-8 / 100
        @test !near_curve.ok
    end

    @testset "check_scaling" begin
        # Near (1,0), ∂f/∂x = 2, ∂f/∂y ≈ 0, so x should dominate the sensitivity.
        f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
        cs = check_scaling(f, nothing, [1.0, 0.0]; grain = 0.01)
        @test cs.idx_max == 1
        @test cs.idx_min == 2
        @test cs.sens_ref > 0
        @test cs.ratio ≥ 1

        # A well-scaled, symmetric example should report :ok.
        f_sym(x, _) = [x[1] - x[2]]
        cs_sym = check_scaling(f_sym, nothing, [0.5, 0.5]; grain = 0.1)
        @test cs_sym.status == :ok
        @test cs_sym.sensitivities[1] ≈ cs_sym.sensitivities[2]
    end

    @testset "score_start_point" begin
        f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
        rng = MersenneTwister(1)

        good = score_start_point(f, nothing, [1.0, 0.0], 0.2; ntrials = 30, rng)
        @test good.ntrials == 30
        @test good.hits == round(Int, good.score * good.ntrials)
        @test good.score > 0.5   # a point on a smooth curve should thread reliably

        far = score_start_point(f, nothing, [10.0, 10.0], 0.2; ntrials = 30, rng)
        @test far.score == 0.0   # far from the zero set, no simplex should thread it
    end

    @testset "measure_noise" begin
        f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
        y0 = [1.0, 0.0]

        @test_throws ArgumentError measure_noise(f, nothing, y0; n = 1)

        # Deterministic residual: zero noise, and downstream stats become trivial.
        mn = measure_noise(f, nothing, y0; n = 5)
        @test mn.sigma == 0.0
        @test mn.snr === missing
        @test mn.fom_floor === missing
        @test mn.suggested_grain === missing

        cs = check_scaling(f, nothing, y0; grain = 0.2)
        mn2 = measure_noise(f, nothing, y0; n = 5, sens_ref = cs.sens_ref)
        @test mn2.snr == Inf
        @test mn2.fom_floor == 0.0
        @test mn2.suggested_grain === missing   # sigma == 0, no grain suggestion possible

        # Injected per-call noise gives a nonzero, finite SNR.
        rng = MersenneTwister(2)
        f_noisy(x, _) = [x[1]^2 + x[2]^2 - 1.0 + 1.0e-3 * randn(rng)]
        mn3 = measure_noise(f_noisy, nothing, y0; n = 50, sens_ref = cs.sens_ref, grain = 0.2)
        @test mn3.sigma > 0
        @test isfinite(mn3.snr)
        @test mn3.suggested_grain isa Float64
    end

    @testset "improve_start_point" begin
        f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
        y_bad = [1.2, 0.1]
        r_before = norm(f(y_bad, nothing))
        y_improved = improve_start_point(f, nothing, y_bad; grain = 0.2)
        r_after = norm(f(y_improved, nothing))

        @test r_after < r_before
        @test y_improved[2] == y_bad[2]   # last coordinate held fixed

        # A point already on the curve should not be moved (much) further off it.
        y_good = [1.0, 0.0]
        y_still_good = improve_start_point(f, nothing, y_good; grain = 0.2)
        @test norm(f(y_still_good, nothing)) < 1.0e-8
    end

    @testset "curr_fom" begin
        f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
        y = [1.1, 0.0]

        raw = curr_fom(f, nothing, y)
        @test raw.fom ≈ 0.21
        @test raw.residual ≈ [0.21]

        cs = check_scaling(f, nothing, [1.0, 0.0]; grain = 0.2)
        normalized = curr_fom(f, nothing, y; sens_ref = cs.sens_ref)
        @test normalized.fom ≈ raw.fom / cs.sens_ref
    end
end
