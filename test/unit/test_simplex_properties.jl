using SimplexContinuation
using Test

@testset "Simplex Properties and Utility Functions" begin
    @testset "Dimension Functions" begin
        # Test simplex_dimension
        s2d = SimplexContinuation.Simplex([[0, 0], [1, 0], [0, 1]])
        @test simplex_dimension(s2d) == 2
        @test s2d.n == 2

        s3d = SimplexContinuation.Simplex([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        @test simplex_dimension(s3d) == 3
        @test s3d.n == 3

        # Test space_dimension
        @test space_dimension(s2d) == 2
        @test s2d.dims == 2
        @test space_dimension(s3d) == 3
        @test s3d.dims == 3

        # Test non-square case (simplex in higher dimensional space)
        s_rect = SimplexContinuation.Simplex{Int}(undef, 2, 3)
        @test simplex_dimension(s_rect) == 2
        @test space_dimension(s_rect) == 3
    end

    @testset "Element Type Function" begin
        s_int = SimplexContinuation.Simplex{Int}([[0, 0], [1, 0], [0, 1]])
        @test eltype(s_int) == Int

        s_float = SimplexContinuation.Simplex{Float64}([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        @test eltype(s_float) == Float64

        s_rational = SimplexContinuation.Simplex{Rational{Int}}([[0 // 1, 0 // 1], [1 // 1, 0 // 1], [0 // 1, 1 // 1]])
        @test eltype(s_rational) == Rational{Int}
    end

    @testset "Sorting Function" begin
        # Test that sort! sorts vertices and returns the simplex
        vertices_unsorted = [[2, 1], [0, 0], [1, 2]]
        s = SimplexContinuation.Simplex(vertices_unsorted)

        # Sort and check return value
        result = sort!(s)
        @test result === s  # Should return the same object

        # Check that vertices are sorted
        @test s.vertices[1] == [0, 0]
        @test s.vertices[2] == [1, 2]
        @test s.vertices[3] == [2, 1]

        # Test with different data types
        vertices_float = [[2.5, 1.0], [0.0, 0.0], [1.0, 2.5]]
        s_float = SimplexContinuation.Simplex(vertices_float)
        sort!(s_float)
        @test s_float.vertices[1] == [0.0, 0.0]
        @test s_float.vertices[2] == [1.0, 2.5]
        @test s_float.vertices[3] == [2.5, 1.0]

        # Test sorting with negative numbers
        vertices_neg = [[-1, -2], [-3, 0], [1, 1]]
        s_neg = SimplexContinuation.Simplex(vertices_neg)
        sort!(s_neg)
        @test s_neg.vertices[1] == [-3, 0]
        @test s_neg.vertices[2] == [-1, -2]
        @test s_neg.vertices[3] == [1, 1]

        # Test sorting with duplicate vertices
        vertices_dup = [[1, 1], [0, 0], [1, 1]]
        s_dup = SimplexContinuation.Simplex(vertices_dup)
        sort!(s_dup)
        @test s_dup.vertices[1] == [0, 0]
        @test s_dup.vertices[2] == [1, 1]
        @test s_dup.vertices[3] == [1, 1]
    end

    @testset "Type Consistency" begin
        # Test that dimension and property functions work with different types
        for T in [Int8, Int16, Int32, Int64, Float32, Float64]
            if T <: Integer
                vertices = [T.([0, 0]), T.([1, 0]), T.([0, 1])]
            else
                vertices = [[T(0), T(0)], [T(1), T(0)], [T(0), T(1)]]
            end
            s = SimplexContinuation.Simplex(vertices)

            @test eltype(s) == T
            @test simplex_dimension(s) == 2
            @test space_dimension(s) == 2

            # Test sorting preserves type
            sort!(s)
            @test eltype(s) == T
        end
    end

    @testset "Memory Safety" begin
        # Test that property functions don't mutate inputs
        original_vertices = [[0, 0], [1, 0], [0, 1]]
        s = SimplexContinuation.Simplex(copy.(original_vertices))

        # These functions should not mutate the original simplex
        simplex_dimension(s)
        space_dimension(s)
        eltype(s)

        # Verify original is unchanged (except sort! which is explicitly mutating)
        @test s.vertices == original_vertices

        # Test that sort! does mutate
        s_for_sort = SimplexContinuation.Simplex([[2, 1], [0, 0], [1, 2]])
        original_order = copy(s_for_sort.vertices)
        sort!(s_for_sort)
        @test s_for_sort.vertices != original_order
    end
end
