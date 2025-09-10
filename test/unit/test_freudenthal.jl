using SimplexContinuation
using Test

@testset "Freudenthal Functions" begin
    @testset "Freudenthal Initial Simplex" begin
        @testset "Basic Construction" begin
            # Test 1D case
            s1d = freudenthal_initial_simplex(1)
            @test s1d.n == 1
            @test s1d.dims == 1
            @test size(s1d.vertices, 2) == 2
            @test s1d.vertices[:, 1] == [0]
            @test s1d.vertices[:, 2] == [1]
            @test is_freudenthal(s1d)

            # Test 2D case
            s2d = freudenthal_initial_simplex(2)
            @test s2d.n == 2
            @test s2d.dims == 2
            @test size(s2d.vertices, 2) == 3
            @test s2d.vertices[:, 1] == [0, 0]
            @test s2d.vertices[:, 2] == [1, 0]
            @test s2d.vertices[:, 3] == [1, 1]
            @test is_freudenthal(s2d)

            # Test 3D case
            s3d = freudenthal_initial_simplex(3)
            @test s3d.n == 3
            @test s3d.dims == 3
            @test size(s3d.vertices, 2) == 4
            @test s3d.vertices[:, 1] == [0, 0, 0]
            @test s3d.vertices[:, 2] == [1, 0, 0]
            @test s3d.vertices[:, 3] == [1, 1, 0]
            @test s3d.vertices[:, 4] == [1, 1, 1]
            @test is_freudenthal(s3d)

            # Test 4D case
            s4d = freudenthal_initial_simplex(4)
            @test s4d.n == 4
            @test s4d.dims == 4
            @test size(s4d.vertices, 2) == 5
            @test s4d.vertices[:, 1] == [0, 0, 0, 0]
            @test s4d.vertices[:, 2] == [1, 0, 0, 0]
            @test s4d.vertices[:, 3] == [1, 1, 0, 0]
            @test s4d.vertices[:, 4] == [1, 1, 1, 0]
            @test s4d.vertices[:, 5] == [1, 1, 1, 1]
            @test is_freudenthal(s4d)
        end

        @testset "Type Specification" begin
            # Test with specific type
            s_int = freudenthal_initial_simplex(Int, 2)
            @test eltype(s_int) == Int
            @test is_freudenthal(s_int)

            s_int8 = freudenthal_initial_simplex(Int8, 2)
            @test eltype(s_int8) == Int8
            @test is_freudenthal(s_int8)

            s_int32 = freudenthal_initial_simplex(Int32, 3)
            @test eltype(s_int32) == Int32
            @test is_freudenthal(s_int32)

            # Test that all coordinates are correct type
            s_typed = freudenthal_initial_simplex(Int16, 2)
            @test all(x -> typeof(x) == Int16, s_typed.vertices)
        end

        @testset "Edge Cases" begin
            # Test error for invalid dimension
            @test_throws ArgumentError freudenthal_initial_simplex(0)
            @test_throws ArgumentError freudenthal_initial_simplex(-1)
        end

        @testset "Large Dimensions" begin
            # Test that construction works for larger dimensions
            for n in 5:8
                s = freudenthal_initial_simplex(n)
                @test s.n == n
                @test s.dims == n
                @test size(s.vertices, 2) == n + 1
                @test is_freudenthal(s)

                # Verify staircase pattern
                for i in 1:(n + 1)
                    vertex = s.vertices[:, i]
                    ones_count = sum(vertex)
                    @test ones_count == i - 1
                end
            end
        end
    end

    @testset "is_freudenthal Function" begin
        @testset "Valid Freudenthal Simplices" begin
            # Test initial simplices
            for n in 1:6
                s = freudenthal_initial_simplex(n)
                @test is_freudenthal(s)
            end

            # Test translated Freudenthal simplex
            s2d = freudenthal_initial_simplex(2)
            offset = [3, 5]
            translated_vertices = [s2d.vertices[:, i] .+ offset for i in 1:size(s2d.vertices, 2)]
            translated_simplex = SimplexContinuation.Simplex(translated_vertices)
            @test is_freudenthal(translated_simplex)

            # Test other valid 2D Freudenthal simplex (y-then-x order)
            valid_2d = SimplexContinuation.Simplex([[0, 0], [0, 1], [1, 1]])
            @test is_freudenthal(valid_2d)

            # Test other valid 3D Freudenthal simplices
            valid_3d_1 = SimplexContinuation.Simplex([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 1, 1]])
            @test is_freudenthal(valid_3d_1)

            valid_3d_2 = SimplexContinuation.Simplex([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1]])
            @test is_freudenthal(valid_3d_2)

            valid_3d_3 = SimplexContinuation.Simplex([[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]])
            @test is_freudenthal(valid_3d_3)

            valid_3d_4 = SimplexContinuation.Simplex([[0, 0, 0], [1, 0, 0], [1, 0, 1], [1, 1, 1]])
            @test is_freudenthal(valid_3d_4)

            valid_3d_5 = SimplexContinuation.Simplex([[0, 0, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1]])
            @test is_freudenthal(valid_3d_5)
        end

        @testset "Invalid Freudenthal Simplices" begin
            # Test non-full-dimensional simplex
            vertices_2d_in_3d = [[0, 0, 0], [1, 0, 0], [1, 1, 0]]
            simplex_2d_in_3d = SimplexContinuation.Simplex(vertices_2d_in_3d)
            @test !is_freudenthal(simplex_2d_in_3d)

            # Test arbitrary simplex
            arbitrary_vertices = [[0, 0], [1, 1], [2, 0]]
            arbitrary_simplex = SimplexContinuation.Simplex(arbitrary_vertices)
            @test !is_freudenthal(arbitrary_simplex)

            # Test simplex with wrong coordinate values
            wrong_coords = [[0, 0], [2, 0], [2, 1]]
            wrong_simplex = SimplexContinuation.Simplex(wrong_coords)
            @test !is_freudenthal(wrong_simplex)

            # Test simplex with negative coordinates
            negative_coords = [[0, 0], [-1, 0], [-1, 1]]
            negative_simplex = SimplexContinuation.Simplex(negative_coords)
            @test !is_freudenthal(negative_simplex)

            # Test degenerate case
            degenerate_vertices = [[0, 0], [1, 0], [0, 0]]
            degenerate_simplex = SimplexContinuation.Simplex(degenerate_vertices)
            @test !is_freudenthal(degenerate_simplex)

            # Test simplex with fractional coordinates
            # Note: is_freudenthal only works with Integer types, so we can't test this directly
            # fractional_coords = [[0, 0], [1, 0], [1, 0.5]]
            # fractional_simplex = SimplexContinuation.Simplex(fractional_coords)
            # This would cause a MethodError since is_freudenthal requires Integer types

            # Test wrong increment pattern
            wrong_pattern = [[0, 0, 0], [1, 1, 0], [1, 1, 1], [1, 1, 1]]
            wrong_simplex = SimplexContinuation.Simplex(wrong_pattern)
            @test !is_freudenthal(wrong_simplex)
        end

        @testset "Different Number Types" begin
            # Test with different number types
            vertices_int8 = [[Int8(0), Int8(0)], [Int8(0), Int8(1)], [Int8(1), Int8(1)]]
            simplex_int8 = SimplexContinuation.Simplex(vertices_int8)
            @test is_freudenthal(simplex_int8)

            vertices_uint = [[UInt(0), UInt(0)], [UInt(0), UInt(1)], [UInt(1), UInt(1)]]
            simplex_uint = SimplexContinuation.Simplex(vertices_uint)
            @test is_freudenthal(simplex_uint)

            vertices_int64 = [[Int64(0), Int64(0)], [Int64(1), Int64(0)], [Int64(1), Int64(1)]]
            simplex_int64 = SimplexContinuation.Simplex(vertices_int64)
            @test is_freudenthal(simplex_int64)

            vertices_float64 = [[Float64(0), Float64(0)], [Float64(1), Float64(0)], [Float64(1), Float64(1)]]
            simplex_float64 = SimplexContinuation.Simplex(vertices_float64)
            @test is_freudenthal(simplex_float64)

            vertices_float64eps = [[Float64(1.0e-12), Float64(-1.0e-12)], [Float64(1 + 1.0e-12), Float64(1.0e-12)], [Float64(1 - 1.0e-12), Float64(1 + 1.0e-12)]]
            simplex_float64eps = SimplexContinuation.Simplex(vertices_float64eps)
            @test is_freudenthal(simplex_float64eps)
        end

        @testset "Translated Freudenthal Simplices" begin
            # Test that translation preserves Freudenthal property
            base_simplex = freudenthal_initial_simplex(3)

            # Test various translation vectors
            translations = [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [2, 3, 1],
                [-1, -2, -3],
                [10, -5, 7],
            ]

            for translation in translations
                translated_vertices = [base_simplex.vertices[:, i] .+ translation for i in 1:size(base_simplex.vertices, 2)]
                translated_simplex = SimplexContinuation.Simplex(translated_vertices)
                @test is_freudenthal(translated_simplex)
            end
        end
    end

    @testset "freudenthal_reflect Function" begin
        @testset "2D Reflections" begin
            s2d = freudenthal_initial_simplex(2)

            # Test all possible reflections
            for facet_idx in 1:3
                reflected = freudenthal_reflect(s2d, facet_idx)
                @test is_freudenthal(reflected)
                @test reflected.n == s2d.n
                @test reflected.dims == s2d.dims
                @test eltype(reflected) == eltype(s2d)
            end

            # Test specific reflection results
            # Reflecting around facet 1 moves to next cube
            reflected_1 = freudenthal_reflect(s2d, 1)
            expected_1 = SimplexContinuation.Simplex([[1, 0], [1, 1], [2, 1]])
            @test reflected_1.vertices == expected_1.vertices

            # Test double reflection returns to original only for internal facets
            reflected_2 = freudenthal_reflect(s2d, 2)
            double_reflected_2 = freudenthal_reflect(reflected_2, 2)
            @test double_reflected_2.vertices == s2d.vertices

            # Test reflection around facet 2 (internal facet - swaps increments)
            expected_2 = SimplexContinuation.Simplex([[0, 0], [0, 1], [1, 1]])
            @test reflected_2.vertices == expected_2.vertices

            # Test reflection around facet 3 moves to previous cube
            reflected_3 = freudenthal_reflect(s2d, 3)
            expected_3 = SimplexContinuation.Simplex([[0, -1], [0, 0], [1, 0]])
            @test reflected_3.vertices == expected_3.vertices
        end

        @testset "3D Reflections" begin
            s3d = freudenthal_initial_simplex(3)

            # Test all possible reflections
            for facet_idx in 1:4
                reflected = freudenthal_reflect(s3d, facet_idx)
                @test is_freudenthal(reflected)
                @test reflected.n == s3d.n
                @test reflected.dims == s3d.dims
                @test eltype(reflected) == eltype(s3d)
            end

            # Test specific reflection results
            reflected_1 = freudenthal_reflect(s3d, 1)
            expected_1 = SimplexContinuation.Simplex([[1, 0, 0], [1, 1, 0], [1, 1, 1], [2, 1, 1]])
            @test reflected_1.vertices == expected_1.vertices

            reflected_2 = freudenthal_reflect(s3d, 2)
            expected_2 = SimplexContinuation.Simplex([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 1, 1]])
            @test reflected_2.vertices == expected_2.vertices

            reflected_3 = freudenthal_reflect(s3d, 3)
            expected_3 = SimplexContinuation.Simplex([[0, 0, 0], [1, 0, 0], [1, 0, 1], [1, 1, 1]])
            @test reflected_3.vertices == expected_3.vertices

            reflected_4 = freudenthal_reflect(s3d, 4)
            expected_4 = SimplexContinuation.Simplex([[0, 0, -1], [0, 0, 0], [1, 0, 0], [1, 1, 0]])
            @test reflected_4.vertices == expected_4.vertices
        end

        @testset "Higher Dimension Reflections" begin
            # Test 4D reflections
            s4d = freudenthal_initial_simplex(4)
            for facet_idx in 1:5
                reflected = freudenthal_reflect(s4d, facet_idx)
                @test is_freudenthal(reflected)
                @test reflected.n == s4d.n
                @test reflected.dims == s4d.dims
            end

            # Test even higher dimensions
            for n in 5:7
                s = freudenthal_initial_simplex(n)
                # Test first and last vertex reflections
                reflected_first = freudenthal_reflect(s, 1)
                reflected_last = freudenthal_reflect(s, n + 1)
                @test is_freudenthal(reflected_first)
                @test is_freudenthal(reflected_last)
            end
        end

        @testset "Reflection Properties" begin
            s3d = freudenthal_initial_simplex(3)

            # Test that reflection is involutive only for internal facets (2 to n)
            # Boundary facets (1 and n+1) move between cubes and are not involutive
            for facet_idx in 2:3  # Internal facets for 3D
                reflected = freudenthal_reflect(s3d, facet_idx)
                double_reflected = freudenthal_reflect(reflected, facet_idx)
                @test double_reflected.vertices == s3d.vertices
            end

            # Test that all reflections produce valid Freudenthal simplices
            for facet_idx in 1:4
                reflected = freudenthal_reflect(s3d, facet_idx)
                @test is_freudenthal(reflected)
            end

            # Test reflection of translated simplex
            offset = [2, 3, 1]
            translated_vertices = [s3d.vertices[:, i] .+ offset for i in 1:size(s3d.vertices, 2)]
            translated_simplex = SimplexContinuation.Simplex(translated_vertices)
            @test is_freudenthal(translated_simplex)

            reflected_translated = freudenthal_reflect(translated_simplex, 1)
            @test is_freudenthal(reflected_translated)

            # Test that translation commutes with reflection
            reflected_then_translated = freudenthal_reflect(s3d, 1)
            reflected_then_translated_vertices = [reflected_then_translated.vertices[:, i] .+ offset for i in 1:size(reflected_then_translated.vertices, 2)]
            reflected_then_translated_simplex = SimplexContinuation.Simplex(reflected_then_translated_vertices)
            @test reflected_translated.vertices == reflected_then_translated_simplex.vertices
        end

        @testset "Type Preservation" begin
            # Test that reflection preserves element type
            for T in [Int8, Int16, Int32, Int64]
                s = freudenthal_initial_simplex(T, 2)
                reflected = freudenthal_reflect(s, 1)
                @test eltype(reflected) == T
                @test is_freudenthal(reflected)
            end
        end

        @testset "Error Handling" begin
            s2d = freudenthal_initial_simplex(2)

            # Test invalid facet indices
            @test_throws ArgumentError freudenthal_reflect(s2d, 0)
            @test_throws ArgumentError freudenthal_reflect(s2d, 4)
            @test_throws ArgumentError freudenthal_reflect(s2d, -1)
            @test_throws ArgumentError freudenthal_reflect(s2d, 100)

            # Test non-full-dimensional simplex
            vertices_2d_in_3d = [[0, 0, 0], [1, 0, 0], [1, 1, 0]]
            simplex_2d_in_3d = SimplexContinuation.Simplex(vertices_2d_in_3d)
            @test_throws ArgumentError freudenthal_reflect(simplex_2d_in_3d, 1)
        end

        @testset "Original Bug Regression Test" begin
            # This test specifically addresses the original issue:
            # "A simplex passed to freudenthal_reflect with a facet_index of
            # either the first or last vertex, then fails the is_freudenthal function test"

            # Test 2D case
            s2d = freudenthal_initial_simplex(2)

            # Before the fix, these would fail the is_freudenthal test
            reflected_first_2d = freudenthal_reflect(s2d, 1)  # First vertex
            reflected_last_2d = freudenthal_reflect(s2d, 3)   # Last vertex

            @test is_freudenthal(reflected_first_2d)
            @test is_freudenthal(reflected_last_2d)

            # Test 3D case
            s3d = freudenthal_initial_simplex(3)

            reflected_first_3d = freudenthal_reflect(s3d, 1)  # First vertex
            reflected_last_3d = freudenthal_reflect(s3d, 4)   # Last vertex

            @test is_freudenthal(reflected_first_3d)
            @test is_freudenthal(reflected_last_3d)

            # Test 4D case
            s4d = freudenthal_initial_simplex(4)

            reflected_first_4d = freudenthal_reflect(s4d, 1)  # First vertex
            reflected_last_4d = freudenthal_reflect(s4d, 5)   # Last vertex

            @test is_freudenthal(reflected_first_4d)
            @test is_freudenthal(reflected_last_4d)
        end

        @testset "Adjacency and Combinatorial Properties" begin
            # Test that reflected simplices share the correct facet
            s2d = freudenthal_initial_simplex(2)
            reflected = freudenthal_reflect(s2d, 2)  # Reflect around facet 2

            # These should share a facet (2 vertices in 2D)
            shared_vertices = 0
            for i in 1:size(s2d.vertices, 2)
                for j in 1:size(reflected.vertices, 2)
                    if s2d.vertices[:, i] == reflected.vertices[:, j]
                        shared_vertices += 1
                    end
                end
            end
            @test shared_vertices == 2  # Should share exactly 2 vertices (an edge)

            # Test combinatorial structure preservation for 3D
            s3d = freudenthal_initial_simplex(3)

            # Test internal facets which should share n vertices (the facet)
            for facet_idx in 2:3  # Internal facets for 3D
                reflected = freudenthal_reflect(s3d, facet_idx)
                shared_count = 0
                for i in 1:size(s3d.vertices, 2)
                    for j in 1:size(reflected.vertices, 2)
                        if s3d.vertices[:, i] == reflected.vertices[:, j]
                            shared_count += 1
                        end
                    end
                end
                @test shared_count == 3  # Should share exactly 3 vertices (a face)
            end

            # Test that boundary facets (1 and n+1) also share exactly 3 vertices (a face)
            for facet_idx in [1, 4]
                reflected = freudenthal_reflect(s3d, facet_idx)
                shared_count = 0
                for i in 1:size(s3d.vertices, 2)
                    for j in 1:size(reflected.vertices, 2)
                        if s3d.vertices[:, i] == reflected.vertices[:, j]
                            shared_count += 1
                        end
                    end
                end
                @test shared_count == 3  # Should share exactly 3 vertices (a face)
            end

            # Test that all reflections produce valid structures
            for facet_idx in 1:4
                reflected = freudenthal_reflect(s3d, facet_idx)
                @test is_freudenthal(reflected)
                @test reflected.n == s3d.n
                @test reflected.dims == s3d.dims
            end
        end
    end

    @testset "Freudenthal Integration Tests" begin
        @testset "Multiple Reflections" begin
            s3d = freudenthal_initial_simplex(3)

            # Test chain of reflections
            current = s3d
            for i in 1:10
                facet_idx = (i % 4) + 1
                current = freudenthal_reflect(current, facet_idx)
                @test is_freudenthal(current)
            end

            # Test different reflection sequences lead to different simplices
            path1 = freudenthal_reflect(freudenthal_reflect(s3d, 1), 2)
            path2 = freudenthal_reflect(freudenthal_reflect(s3d, 2), 1)
            @test is_freudenthal(path1)
            @test is_freudenthal(path2)
        end

        @testset "Complex Workflows" begin
            # Workflow 1: Create, reflect, check properties
            s = freudenthal_initial_simplex(3)
            @test simplex_dimension(s) == 3
            @test space_dimension(s) == 3
            @test is_freudenthal(s)

            reflected = freudenthal_reflect(s, 2)
            @test simplex_dimension(reflected) == 3
            @test space_dimension(reflected) == 3
            @test is_freudenthal(reflected)

            # Workflow 2: Multiple reflections with property checks
            current_simplex = freudenthal_initial_simplex(2)
            for _ in 1:5
                for facet in 1:3
                    current_simplex = freudenthal_reflect(current_simplex, facet)
                    @test is_freudenthal(current_simplex)
                    @test simplex_dimension(current_simplex) == 2
                    @test space_dimension(current_simplex) == 2
                end
            end
        end

        @testset "Performance with Large Dimensions" begin
            # Test that Freudenthal functions scale reasonably
            for n in [8, 9, 10]
                s = freudenthal_initial_simplex(n)
                @test is_freudenthal(s)

                # Test a few reflections
                reflected_first = freudenthal_reflect(s, 1)
                reflected_last = freudenthal_reflect(s, n + 1)

                @test is_freudenthal(reflected_first)
                @test is_freudenthal(reflected_last)
            end
        end
    end
end
