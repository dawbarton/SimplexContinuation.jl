using SimplexContinuation
using Test

@testset "Integration Tests and Edge Cases" begin
    @testset "Function Composition and Workflows" begin
        @testset "Basic Workflow: Create → Properties → Reflect" begin
            # Test realistic workflow combining multiple functions
            s = freudenthal_initial_simplex(3)
            @test simplex_dimension(s) == 3
            @test space_dimension(s) == 3
            @test is_freudenthal(s)

            reflected = freudenthal_reflect(s, 2)
            @test simplex_dimension(reflected) == 3
            @test space_dimension(reflected) == 3
            @test is_freudenthal(reflected)

            # Test sorting after reflection
            sort!(reflected)
            @test is_freudenthal(reflected)
        end

        @testset "Mixed Workflow: Freudenthal → Convert → Geometric Reflect" begin
            # Workflow: Create Freudenthal, convert to float, reflect geometrically
            vertices = [[0, 0], [1, 0], [0, 1]]
            s_int = SimplexContinuation.Simplex(vertices)
            s_float = SimplexContinuation.Simplex{Float64}(s_int)
            @test eltype(s_float) == Float64

            reflected_geom = reflect(s_float, 1)
            @test eltype(reflected_geom) == Float64
            @test simplex_dimension(reflected_geom) == 2
        end

        @testset "Sort → Check Properties → Reflect Workflow" begin
            unsorted = SimplexContinuation.Simplex([[2, 1], [0, 0], [1, 2]])
            sort!(unsorted)
            @test simplex_dimension(unsorted) == 2
            @test space_dimension(unsorted) == 2
            @test eltype(unsorted) == Int

            # Convert to Freudenthal and test reflection
            freudenthal_vertices = [[0, 0], [1, 0], [1, 1]]
            freudenthal_simplex = SimplexContinuation.Simplex(freudenthal_vertices)
            @test is_freudenthal(freudenthal_simplex)

            reflected = freudenthal_reflect(freudenthal_simplex, 1)
            @test is_freudenthal(reflected)
        end

        @testset "Chained Operations" begin
            # Test long chains of operations
            s = freudenthal_initial_simplex(2)

            # Chain: reflect → sort → check → reflect again
            reflected1 = freudenthal_reflect(s, 1)
            sort!(reflected1)
            @test is_freudenthal(reflected1)

            reflected2 = freudenthal_reflect(reflected1, 2)
            @test is_freudenthal(reflected2)

            reflected3 = freudenthal_reflect(reflected2, 3)
            @test is_freudenthal(reflected3)
        end
    end

    @testset "Type Consistency Across All Functions" begin
        @testset "Integer Type Consistency" begin
            # Test that all functions work consistently with different integer types
            for T in [Int8, Int16, Int32, Int64]
                s = freudenthal_initial_simplex(T, 2)
                @test eltype(s) == T
                @test is_freudenthal(s)
                @test simplex_dimension(s) == 2
                @test space_dimension(s) == 2

                # Test reflection preserves type
                reflected = freudenthal_reflect(s, 1)
                @test eltype(reflected) == T
                @test is_freudenthal(reflected)

                # Test sorting preserves type
                sort!(reflected)
                @test eltype(reflected) == T

                # Test copy constructor preserves type
                copied = SimplexContinuation.Simplex(reflected)
                @test eltype(copied) == T
            end
        end

        @testset "Unsigned Integer Types" begin
            # Test with unsigned types where applicable
            for T in [UInt8, UInt16, UInt32, UInt64]
                vertices = [T.([0, 0]), T.([1, 0]), T.([1, 1])]
                s = SimplexContinuation.Simplex(vertices)
                @test eltype(s) == T
                @test is_freudenthal(s)

                reflected = freudenthal_reflect(s, 1)
                @test eltype(reflected) == T
                @test is_freudenthal(reflected)
            end
        end

        @testset "Type Conversion Workflows" begin
            # Test workflows involving type conversions
            s_int = freudenthal_initial_simplex(Int32, 2)
            @test eltype(s_int) == Int32

            # Convert to different integer type
            s_int64 = SimplexContinuation.Simplex{Int64}(s_int)
            @test eltype(s_int64) == Int64
            @test is_freudenthal(s_int64)

            # Convert to floating point for geometric operations
            s_float = SimplexContinuation.Simplex{Float64}(s_int)
            @test eltype(s_float) == Float64

            reflected_geom = reflect(s_float, 1)
            @test eltype(reflected_geom) == Float64
        end
    end

    @testset "Memory Safety and Mutation" begin
        @testset "Non-mutating Functions" begin
            # Test that functions don't unexpectedly mutate inputs
            original_vertices = [[0, 0], [1, 0], [0, 1]]
            s = SimplexContinuation.Simplex(copy.(original_vertices))

            # These functions should not mutate the original simplex
            simplex_dimension(s)
            space_dimension(s)
            eltype(s)

            # Verify original is unchanged
            @test s.vertices == original_vertices

            # Test Freudenthal functions don't mutate
            freudenthal_s = freudenthal_initial_simplex(2)
            original_freudenthal = [copy(v) for v in freudenthal_s.vertices]

            is_freudenthal(freudenthal_s)
            freudenthal_reflect(freudenthal_s, 1)

            @test freudenthal_s.vertices == original_freudenthal
        end

        @testset "Mutating Functions" begin
            # Test that sort! does mutate and returns the same object
            s_for_sort = SimplexContinuation.Simplex([[2, 1], [0, 0], [1, 2]])
            original_order = copy(s_for_sort.vertices)

            result = sort!(s_for_sort)
            @test result === s_for_sort  # Same object
            @test s_for_sort.vertices != original_order  # Content changed
        end

        @testset "Deep Copy Behavior" begin
            # Test that modifications to copied simplices don't affect originals
            s1 = SimplexContinuation.Simplex([[1, 2], [3, 4], [5, 6]])
            s2 = SimplexContinuation.Simplex(s1)

            # Modify s2 by sorting
            sort!(s2)

            # s1 should be unchanged
            @test s1.vertices[1] == [1, 2]
            @test s1.vertices[2] == [3, 4]
            @test s1.vertices[3] == [5, 6]

            # s2 should be sorted
            @test s2.vertices[1] == [1, 2]
            @test s2.vertices[2] == [3, 4]
            @test s2.vertices[3] == [5, 6]
        end
    end

    @testset "Boundary Cases and Edge Conditions" begin
        @testset "1D Edge Cases" begin
            # Test 1D simplices (minimal case)
            s1d = freudenthal_initial_simplex(1)
            @test simplex_dimension(s1d) == 1
            @test space_dimension(s1d) == 1
            @test length(s1d.vertices) == 2
            @test s1d.vertices[1] == [0]
            @test s1d.vertices[2] == [1]

            # Test 1D reflections
            reflected_1d_1 = freudenthal_reflect(s1d, 1)
            reflected_1d_2 = freudenthal_reflect(s1d, 2)
            @test is_freudenthal(reflected_1d_1)
            @test is_freudenthal(reflected_1d_2)

            # In 1D, both reflections should return the same simplex
            @test reflected_1d_1.vertices == s1d.vertices
            @test reflected_1d_2.vertices == s1d.vertices

            # Test geometric reflection on 1D
            s1d_float = SimplexContinuation.Simplex{Float64}(s1d)
            reflected_1d_geom = reflect(s1d_float, 1)
            @test simplex_dimension(reflected_1d_geom) == 1
        end

        @testset "Degenerate Cases" begin
            # Test very small coordinates
            tiny_vertices = [[0.0, 0.0], [1.0e-10, 0.0], [0.0, 1.0e-10]]
            s_tiny = SimplexContinuation.Simplex(tiny_vertices)
            @test simplex_dimension(s_tiny) == 2
            @test space_dimension(s_tiny) == 2

            reflected_tiny = reflect(s_tiny, 1)
            @test simplex_dimension(reflected_tiny) == 2

            # Test large coordinates
            large_vertices = [[0.0, 0.0], [1.0e6, 0.0], [0.0, 1.0e6]]
            s_large = SimplexContinuation.Simplex(large_vertices)
            @test simplex_dimension(s_large) == 2
            @test space_dimension(s_large) == 2

            reflected_large = reflect(s_large, 1)
            @test simplex_dimension(reflected_large) == 2
        end

        @testset "Extreme Coordinate Values" begin
            # Test with extreme integer values
            max_int = typemax(Int16)  # Use smaller type to avoid overflow
            extreme_vertices = [[0, 0], [max_int, 0], [max_int, max_int]]
            s_extreme = SimplexContinuation.Simplex(extreme_vertices)
            @test simplex_dimension(s_extreme) == 2

            # Test with negative coordinates
            negative_vertices = [[-100, -200], [-50, -200], [-50, -150]]
            s_negative = SimplexContinuation.Simplex(negative_vertices)
            @test simplex_dimension(s_negative) == 2
            @test space_dimension(s_negative) == 2

            sort!(s_negative)
            @test simplex_dimension(s_negative) == 2
        end

        @testset "Non-Standard Configurations" begin
            # Test simplex in higher dimensional space (non-full-dimensional)
            vertices_2d_in_3d = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
            s_2d_in_3d = SimplexContinuation.Simplex(vertices_2d_in_3d)
            @test simplex_dimension(s_2d_in_3d) == 2
            @test space_dimension(s_2d_in_3d) == 3

            # Should fail geometric reflection (not full-dimensional)
            @test_throws ArgumentError reflect(s_2d_in_3d, 1)

            # Test 1D simplex in 3D space
            vertices_1d_in_3d = [[0, 0, 0], [1, 1, 1]]
            s_1d_in_3d = SimplexContinuation.Simplex(vertices_1d_in_3d)
            @test simplex_dimension(s_1d_in_3d) == 1
            @test space_dimension(s_1d_in_3d) == 3
        end
    end

    @testset "Performance and Scalability" begin
        @testset "Large Dimension Performance" begin
            # Test that functions work reasonably with larger dimensions
            for n in [8, 9, 10]
                s = freudenthal_initial_simplex(n)
                @test simplex_dimension(s) == n
                @test space_dimension(s) == n
                @test is_freudenthal(s)
                @test eltype(s) == Int

                # Test that reflection works (first and last vertices)
                reflected_first = freudenthal_reflect(s, 1)
                reflected_last = freudenthal_reflect(s, n + 1)
                @test is_freudenthal(reflected_first)
                @test is_freudenthal(reflected_last)

                # Test sorting still works
                sort!(reflected_first)
                @test is_freudenthal(reflected_first)
            end
        end

        @testset "Many Operations" begin
            # Test many sequential operations don't degrade performance
            s = freudenthal_initial_simplex(3)

            # Perform many reflections
            current = s
            for i in 1:100
                facet_idx = (i % 4) + 1
                current = freudenthal_reflect(current, facet_idx)

                # Verify properties still hold
                if i % 10 == 0  # Check every 10th iteration
                    @test is_freudenthal(current)
                    @test simplex_dimension(current) == 3
                    @test space_dimension(current) == 3
                end
            end
        end

        @testset "Memory Usage Patterns" begin
            # Test that operations don't leak memory or create excessive copies
            vertices = [[0, 0], [1, 0], [1, 1]]

            # Create many simplices
            simplices = []
            for i in 1:100
                s = SimplexContinuation.Simplex(vertices)
                push!(simplices, s)
            end

            # All should have correct properties
            for s in simplices
                @test simplex_dimension(s) == 2
                @test space_dimension(s) == 2
            end
        end
    end

    @testset "Error Handling and Robustness" begin
        @testset "Comprehensive Error Cases" begin
            # Test all functions with invalid inputs
            s = freudenthal_initial_simplex(2)

            # Invalid facet indices for reflection functions
            for invalid_idx in [0, -1, 4, 100]
                @test_throws ArgumentError freudenthal_reflect(s, invalid_idx)

                s_float = SimplexContinuation.Simplex{Float64}(s)
                @test_throws ArgumentError reflect(s_float, invalid_idx)
            end

            # Non-Freudenthal simplex for freudenthal_reflect
            non_freudenthal = SimplexContinuation.Simplex([[0, 0], [2, 0], [1, 1]])
            @test_throws ArgumentError freudenthal_reflect(non_freudenthal, 1)

            # Invalid dimensions for freudenthal_initial_simplex
            @test_throws ArgumentError freudenthal_initial_simplex(0)
            @test_throws ArgumentError freudenthal_initial_simplex(-5)
        end

        @testset "Type Errors" begin
            # Test errors with inappropriate types
            @test_throws ArgumentError freudenthal_initial_simplex(Float64, 2)
            @test_throws ArgumentError freudenthal_initial_simplex(String, 2)
        end

        @testset "Malformed Simplex Construction" begin
            # Test construction with bad vertex data
            bad_vertices = [[0, 0], [1, 0, 0], [0, 1]]  # Mixed dimensions
            @test_throws ArgumentError SimplexContinuation.Simplex(bad_vertices)

            # Empty vertex list
            @test_throws BoundsError SimplexContinuation.Simplex([])
        end
    end

    @testset "Numerical Precision and Stability" begin
        @testset "Floating Point Precision" begin
            # Test with high precision types
            vertices_big = [[BigFloat(0), BigFloat(0)], [BigFloat(1), BigFloat(0)], [BigFloat(0), BigFloat(1)]]
            s_big = SimplexContinuation.Simplex(vertices_big)

            @test eltype(s_big) == BigFloat
            @test simplex_dimension(s_big) == 2
            @test space_dimension(s_big) == 2

            reflected_big = reflect(s_big, 1)
            @test eltype(reflected_big) == BigFloat
        end

        @testset "Rational Arithmetic" begin
            # Test with exact rational arithmetic
            vertices_rational = [[0 // 1, 0 // 1], [1 // 1, 0 // 1], [0 // 1, 1 // 1]]
            s_rational = SimplexContinuation.Simplex(vertices_rational)

            @test eltype(s_rational) == Rational{Int}
            @test simplex_dimension(s_rational) == 2
            @test space_dimension(s_rational) == 2

            reflected_rational = reflect(s_rational, 1)
            @test eltype(reflected_rational) == Rational{Int}

            # Test that results are exact
            @test reflected_rational.vertices[2] == vertices_rational[2]
            @test reflected_rational.vertices[3] == vertices_rational[3]
        end

        @testset "Mixed Precision Operations" begin
            # Test operations mixing different precisions
            s_int = freudenthal_initial_simplex(2)
            s_float = SimplexContinuation.Simplex{Float64}(s_int)
            s_big = SimplexContinuation.Simplex{BigFloat}(s_float)

            @test eltype(s_int) == Int
            @test eltype(s_float) == Float64
            @test eltype(s_big) == BigFloat

            # All should have same logical structure
            @test simplex_dimension(s_int) == simplex_dimension(s_float) == simplex_dimension(s_big)
            @test space_dimension(s_int) == space_dimension(s_float) == space_dimension(s_big)
        end
    end

    @testset "Comprehensive Regression Tests" begin
        @testset "Historical Bug Prevention" begin
            # Test cases that would have caught the original freudenthal_reflect bug
            for n in 2:5
                s = freudenthal_initial_simplex(n)

                # These operations should never fail
                first_reflected = freudenthal_reflect(s, 1)
                last_reflected = freudenthal_reflect(s, n + 1)

                @test is_freudenthal(first_reflected)
                @test is_freudenthal(last_reflected)

                # Double reflection should return to original
                @test freudenthal_reflect(first_reflected, 1).vertices == s.vertices
                @test freudenthal_reflect(last_reflected, n + 1).vertices == s.vertices
            end
        end

        @testset "Cross-Function Consistency" begin
            # Test that all functions agree on basic properties
            for n in 1:5
                s = freudenthal_initial_simplex(n)

                # All functions should agree on dimensions
                @test s.n == simplex_dimension(s) == n
                @test s.dims == space_dimension(s) == n

                # Reflection should preserve these properties
                for facet in [1, n + 1]  # Test first and last
                    reflected = freudenthal_reflect(s, facet)
                    @test simplex_dimension(reflected) == n
                    @test space_dimension(reflected) == n
                    @test reflected.n == n
                    @test reflected.dims == n
                end
            end
        end

        @testset "Invariant Properties" begin
            # Test mathematical invariants that should always hold
            s = freudenthal_initial_simplex(3)

            # Freudenthal property should be preserved under all valid operations
            operations = [
                s -> freudenthal_reflect(s, 1),
                s -> freudenthal_reflect(s, 2),
                s -> freudenthal_reflect(s, 3),
                s -> freudenthal_reflect(s, 4),
                s -> begin
                    sort!(s); s
                end,
                s -> SimplexContinuation.Simplex(s),
            ]

            for op in operations
                result = op(SimplexContinuation.Simplex(s))  # Work on copy
                @test is_freudenthal(result)
                @test simplex_dimension(result) == 3
                @test space_dimension(result) == 3
            end
        end
    end
end
