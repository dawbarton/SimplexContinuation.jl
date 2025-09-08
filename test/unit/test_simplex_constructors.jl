using SimplexContinuation
using Test

@testset "Simplex Constructors" begin
    @testset "Basic Construction from Vertices" begin
        # Test 2D simplex construction
        vertices_2d = [[0, 0], [1, 0], [0, 1]]
        s2d = SimplexContinuation.Simplex(vertices_2d)
        @test s2d.n == 2
        @test s2d.dims == 2
        @test length(s2d.vertices) == 3
        @test s2d.vertices[1] == [0, 0]
        @test s2d.vertices[2] == [1, 0]
        @test s2d.vertices[3] == [0, 1]

        # Test 3D simplex construction
        vertices_3d = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        s3d = SimplexContinuation.Simplex(vertices_3d)
        @test s3d.n == 3
        @test s3d.dims == 3
        @test length(s3d.vertices) == 4

        # Test mixed integer types get promoted
        mixed_vertices = [[0, 0], [1.0, 0], [0, 1]]
        s_mixed = SimplexContinuation.Simplex(mixed_vertices)
        @test eltype(s_mixed) == Float64
    end

    @testset "Type-Specific Construction" begin
        # Test explicit type specification
        vertices = [[0, 0], [1, 0], [0, 1]]
        s_int = SimplexContinuation.Simplex{Int}(vertices)
        @test eltype(s_int) == Int
        @test s_int.vertices[1] == [0, 0]

        s_float = SimplexContinuation.Simplex{Float64}(vertices)
        @test eltype(s_float) == Float64
        @test s_float.vertices[1] == [0.0, 0.0]

        # Test construction with different numeric types
        s_int32 = SimplexContinuation.Simplex{Int32}(vertices)
        @test eltype(s_int32) == Int32
    end

    @testset "Uninitialized Construction" begin
        # Test uninitialized construction with square dimensions
        s_uninit = SimplexContinuation.Simplex{Float64}(undef, 2, 2)
        @test s_uninit.n == 2
        @test s_uninit.dims == 2
        @test length(s_uninit.vertices) == 3
        @test length(s_uninit.vertices[1]) == 2

        # Test uninitialized construction with different dimensions
        s_uninit_rect = SimplexContinuation.Simplex{Int}(undef, 2, 3)
        @test s_uninit_rect.n == 2
        @test s_uninit_rect.dims == 3
        @test length(s_uninit_rect.vertices) == 3
        @test length(s_uninit_rect.vertices[1]) == 3
    end

    @testset "Copy Constructor" begin
        # Test copying between same types
        original = SimplexContinuation.Simplex([[1, 2], [3, 4], [5, 6]])
        copied = SimplexContinuation.Simplex(original)
        @test copied.n == original.n
        @test copied.dims == original.dims
        @test copied.vertices == original.vertices
        @test eltype(copied) == eltype(original)

        # Test type conversion during copy
        s_int = SimplexContinuation.Simplex([[1, 2], [3, 4], [5, 6]])
        s_float = SimplexContinuation.Simplex{Float64}(s_int)
        @test eltype(s_float) == Float64
        @test s_float.vertices[1] == [1.0, 2.0]
        @test s_float.vertices[2] == [3.0, 4.0]
        @test s_float.vertices[3] == [5.0, 6.0]
    end

    @testset "Error Handling" begin
        # Test mismatched vertex dimensions
        bad_vertices = [[0, 0], [1, 0, 0], [0, 1]]  # Mixed 2D and 3D
        @test_throws ArgumentError SimplexContinuation.Simplex(bad_vertices)

        # Test empty vertices
        @test_throws BoundsError SimplexContinuation.Simplex([])
    end
end
