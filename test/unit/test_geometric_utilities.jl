using SimplexContinuation
using Test

@testset "Geometric Utilities" begin
    @testset "get_facet" begin
        s = SimplexContinuation.Simplex([[0, 0], [1, 0], [0, 1]])

        # Drop vertex 1: should give vertices 2 and 3
        f1 = get_facet(s, 1)
        @test f1.n == 1
        @test f1.dims == 2
        @test size(f1.vertices, 2) == 2
        @test f1.vertices[:, 1] == [1, 0]
        @test f1.vertices[:, 2] == [0, 1]

        # Drop vertex 2: should give vertices 1 and 3
        f2 = get_facet(s, 2)
        @test f2.vertices[:, 1] == [0, 0]
        @test f2.vertices[:, 2] == [0, 1]

        # Drop vertex 3: should give vertices 1 and 2
        f3 = get_facet(s, 3)
        @test f3.vertices[:, 1] == [0, 0]
        @test f3.vertices[:, 2] == [1, 0]

        # 3D tetrahedron
        s3 = SimplexContinuation.Simplex([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
        f = get_facet(s3, 2)
        @test f.n == 2
        @test f.dims == 3
        @test size(f.vertices, 2) == 3

        # Error on bad index
        @test_throws ArgumentError get_facet(s, 0)
        @test_throws ArgumentError get_facet(s, 4)
    end

    @testset "barycenter" begin
        # Equilateral-ish triangle: exact centroid
        s = SimplexContinuation.Simplex([[0, 0], [3, 0], [0, 3]])
        bc = barycenter(s)
        @test bc ≈ [1.0, 1.0]

        # Unit Freudenthal 2-simplex: vertices (0,0),(1,0),(1,1)
        sf = freudenthal_initial_simplex(Float64, 2)
        bc_f = barycenter(sf)
        @test bc_f ≈ [2/3, 1/3]

        # 1D line segment
        s1 = SimplexContinuation.Simplex([[0.0], [4.0]])
        @test barycenter(s1) ≈ [2.0]

        # 3D tetrahedron
        s3 = SimplexContinuation.Simplex{Float64}([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
        bc3 = barycenter(s3)
        @test bc3 ≈ [0.25, 0.25, 0.25]

        # Type preservation: integer simplex → rational-style output
        si = SimplexContinuation.Simplex([[0, 0], [2, 0], [0, 2]])
        @test barycenter(si) == [2/3, 2/3]
    end
end
