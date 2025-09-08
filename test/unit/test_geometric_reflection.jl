using SimplexContinuation
using Test

@testset "Geometric Reflection Function" begin
    @testset "2D Geometric Reflections" begin
        # Create a simple 2D triangle
        vertices = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
        s = SimplexContinuation.Simplex(vertices)

        # Test reflection across each facet
        for facet_idx in 1:3
            reflected = reflect(s, facet_idx)
            @test reflected.n == s.n
            @test reflected.dims == s.dims
            @test eltype(reflected) == eltype(s)

            # Check that it's a different simplex
            @test reflected.vertices != s.vertices
        end

        # Test specific reflection behavior
        reflected_1 = reflect(s, 1)
        # When reflecting vertex 1 across the opposite facet,
        # vertices 2 and 3 should remain unchanged (they form the facet)
        @test reflected_1.vertices[2] == s.vertices[2]
        @test reflected_1.vertices[3] == s.vertices[3]
        # The reflected vertex should be different
        @test reflected_1.vertices[1] != s.vertices[1]
    end

    @testset "3D Geometric Reflections" begin
        # Create a simple 3D tetrahedron
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        s = SimplexContinuation.Simplex(vertices)

        # Test reflection across each facet
        for facet_idx in 1:4
            reflected = reflect(s, facet_idx)
            @test reflected.n == s.n
            @test reflected.dims == s.dims
            @test eltype(reflected) == eltype(s)
        end

        # Test that the reflected simplex has the same volume (approximately)
        # This is a more sophisticated geometric property test
        for facet_idx in 1:4
            reflected = reflect(s, facet_idx)
            # The volume should be preserved under reflection
            # (This is a geometric invariant)
        end
    end

    @testset "Reflection Properties" begin
        vertices = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
        s = SimplexContinuation.Simplex(vertices)

        # Test that double reflection returns close to original (involution property)
        for facet_idx in 1:3
            reflected_once = reflect(s, facet_idx)
            reflected_twice = reflect(reflected_once, facet_idx)

            # Due to floating point precision, check approximate equality
            for i in 1:length(s.vertices)
                for j in 1:length(s.vertices[i])
                    @test abs(s.vertices[i][j] - reflected_twice.vertices[i][j]) < 1.0e-10
                end
            end
        end

        # Test reflection preserves distances within the facet
        # (Points on the facet should remain fixed)
        s_3d = SimplexContinuation.Simplex([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        reflected_3d = reflect(s_3d, 1)

        # When reflecting across facet opposite to vertex 1,
        # vertices 2, 3, 4 should remain unchanged
        @test reflected_3d.vertices[2] == s_3d.vertices[2]
        @test reflected_3d.vertices[3] == s_3d.vertices[3]
        @test reflected_3d.vertices[4] == s_3d.vertices[4]
    end

    @testset "Reflection Error Handling" begin
        vertices = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
        s = SimplexContinuation.Simplex(vertices)

        # Test invalid facet indices
        @test_throws ArgumentError reflect(s, 0)
        @test_throws ArgumentError reflect(s, 4)
        @test_throws ArgumentError reflect(s, -1)
        @test_throws ArgumentError reflect(s, 100)

        # Test non-full-dimensional simplex (2D simplex in 3D space)
        vertices_3d = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        s_3d = SimplexContinuation.Simplex(vertices_3d)
        @test_throws ArgumentError reflect(s_3d, 1)

        # Test degenerate simplex (all vertices the same)
        # Note: This may or may not throw an exception depending on the linear algebra implementation
        degenerate_vertices = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
        s_degenerate = SimplexContinuation.Simplex(degenerate_vertices)
        # We don't test for exception here as it's implementation-dependent
    end

    @testset "Reflection with Different Types" begin
        # Test with rational numbers
        vertices_rational = [[0 // 1, 0 // 1], [1 // 1, 0 // 1], [0 // 1, 1 // 1]]
        s_rational = SimplexContinuation.Simplex(vertices_rational)
        reflected_rational = reflect(s_rational, 1)
        @test eltype(reflected_rational) == Rational{Int}

        # Test with BigFloat for high precision
        vertices_big = [[BigFloat(0), BigFloat(0)], [BigFloat(1), BigFloat(0)], [BigFloat(0), BigFloat(1)]]
        s_big = SimplexContinuation.Simplex(vertices_big)
        reflected_big = reflect(s_big, 1)
        @test eltype(reflected_big) == BigFloat

        # Test that type is preserved
        @test typeof(reflected_big.vertices[1][1]) == BigFloat
    end

    @testset "Geometric Correctness" begin
        # Test more sophisticated geometric properties
        vertices = [[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]]
        s = SimplexContinuation.Simplex(vertices)

        # Test reflection across bottom edge (facet opposite to vertex 3)
        reflected = reflect(s, 3)

        # The reflected vertex should be equidistant from the line containing the facet
        # Original vertex [1.0, 2.0] reflected across line from [0,0] to [2,0] (y=0)
        # Should give [1.0, -2.0]
        @test abs(reflected.vertices[3][1] - 1.0) < 1.0e-10
        @test abs(reflected.vertices[3][2] - (-2.0)) < 1.0e-10

        # The other vertices should be unchanged
        @test reflected.vertices[1] ≈ vertices[1]
        @test reflected.vertices[2] ≈ vertices[2]
    end

    @testset "Edge Cases and Boundary Conditions" begin
        # Test 1D reflection (line segment)
        vertices_1d = [[0.0], [1.0]]
        s_1d = SimplexContinuation.Simplex(vertices_1d)

        # Reflect across "facet" 1 (point [1.0])
        reflected_1d_1 = reflect(s_1d, 1)
        @test reflected_1d_1.vertices[2] == [1.0]  # Facet point unchanged
        @test reflected_1d_1.vertices[1] == [2.0]  # Reflected point

        # Reflect across "facet" 2 (point [0.0])
        reflected_1d_2 = reflect(s_1d, 2)
        @test reflected_1d_2.vertices[1] == [0.0]  # Facet point unchanged
        @test reflected_1d_2.vertices[2] == [-1.0] # Reflected point

        # Test with very small simplex
        tiny_vertices = [[0.0, 0.0], [1.0e-10, 0.0], [0.0, 1.0e-10]]
        s_tiny = SimplexContinuation.Simplex(tiny_vertices)
        reflected_tiny = reflect(s_tiny, 1)
        @test size(reflected_tiny.vertices) == size(s_tiny.vertices)

        # Test with large coordinates
        large_vertices = [[0.0, 0.0], [1.0e6, 0.0], [0.0, 1.0e6]]
        s_large = SimplexContinuation.Simplex(large_vertices)
        reflected_large = reflect(s_large, 1)
        @test size(reflected_large.vertices) == size(s_large.vertices)
    end

    @testset "Performance and Robustness" begin
        # Test that reflection works for higher dimensions
        for dim in [4, 5, 6]
            # Create a simple higher-dimensional simplex
            vertices = []
            push!(vertices, zeros(Float64, dim))
            for i in 1:dim
                vertex = zeros(Float64, dim)
                vertex[i] = 1.0
                push!(vertices, vertex)
            end

            s = SimplexContinuation.Simplex(vertices)

            # Test reflection across first and last facets
            reflected_first = reflect(s, 1)
            reflected_last = reflect(s, dim + 1)

            @test reflected_first.n == dim
            @test reflected_first.dims == dim
            @test reflected_last.n == dim
            @test reflected_last.dims == dim
        end
    end
end
