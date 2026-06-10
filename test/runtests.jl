using SimplexContinuation
using Test

# Main test runner for SimplexContinuation.jl
# This file includes all unit tests organized into separate files

@testset "SimplexContinuation.jl" begin
    include("unit/test_simplex_constructors.jl")
    include("unit/test_simplex_properties.jl")
    include("unit/test_geometric_reflection.jl")
    include("unit/test_freudenthal.jl")
    include("unit/test_geometric_utilities.jl")
    include("integration/test_continuation.jl")
end
