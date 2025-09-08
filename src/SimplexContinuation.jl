module SimplexContinuation

using FixedSizeArrays: FixedSizeArray, FixedSizeMatrixDefault, FixedSizeVector, FixedSizeVectorDefault
using LinearAlgebra: qr, qr!, nullspace, dot

const FSVector{T} = FixedSizeVectorDefault{T}

export Simplex
export simplex_dimension, space_dimension, reflect, is_freudenthal, freudenthal_initial_simplex, freudenthal_reflect

"""
    Simplex{T}(vertices)

A representation of an `n`-simplex that lives in dimension `dims` which may or
may not be the same as `n`. Both `n` and `dims` are inferred from the size of
`vertices` and each corresponding vertex. If not specifed, `T` will be inferred.
"""
struct Simplex{T}
    n::Int
    dims::Int
    vertices::FixedSizeMatrixDefault{T}

    function Simplex(vertices)
        T = Union{}
        for vertex in vertices
            T = promote_type(T, eltype(vertex))
        end
        return Simplex{T}(vertices)
    end
    function Simplex{T}(vertices) where {T}
        n = length(vertices) - 1
        dims = length(first(vertices))
        new_vertices = FixedSizeMatrixDefault{T}(undef, dims, n + 1)
        for (i, vertex) in enumerate(vertices)
            if length(vertex) != dims
                throw(ArgumentError("All vertices must have the same dimension"))
            else
                new_vertices[:, i] .= vertex
            end
        end
        return new{T}(n, dims, new_vertices)
    end
    function Simplex{T}(::UndefInitializer, n::Integer, dims::Integer = n) where {T}
        new_vertices = FixedSizeMatrixDefault{T}(undef, dims, n + 1)
        return new{T}(n, dims, new_vertices)
    end
end

function Simplex{NT}(simplex::Simplex{T}) where {NT, T}
    new_simplex = Simplex{NT}(undef, simplex.n, simplex.dims)
    copyto!(new_simplex.vertices, simplex.vertices)
    return new_simplex
end
(Simplex(simplex::Simplex{T}) where {T}) = Simplex{T}(simplex)

Base.sort!(simplex::Simplex) = (copyto!(simplex.vertices, sortslices(simplex.vertices, dims = 2)); return simplex)
(Base.eltype(simplex::Simplex{T}) where {T}) = T

"""
    simplex_dimension(simplex::Simplex)

Returns the dimension `n` of the simplex (i.e., an `n`-simplex has `n+1` vertices).
"""
simplex_dimension(simplex::Simplex) = simplex.n

"""
    space_dimension(simplex::Simplex)

Returns the dimension `dims` of the space in which the simplex lives.
"""
space_dimension(simplex::Simplex) = simplex.dims

"""
    reflect(simplex::Simplex, facet_index)

Returns a new simplex that is `simplex` geometrically reflected around one of its facets.
"""
function reflect(simplex::Simplex{T}, facet_index) where {T}
    if (facet_index < 1) || (facet_index > simplex.n + 1)
        throw(ArgumentError("Facet index must be between 1 and $(simplex.n + 1)"))
    end
    if simplex.n != simplex.dims
        throw(ArgumentError("Reflection only implemented for full-dimensional simplices"))
    end

    # Take the first vertex not equal to the facet index to be the origin and find the edge vectors
    origin_index = facet_index == 1 ? 2 : 1
    origin_vertex = simplex.vertices[:, origin_index]
    edge_vectors = FixedSizeMatrixDefault{T}(undef, simplex.dims, simplex.n - 1)
    j = 1
    for i in 1:(simplex.n + 1)
        if i != facet_index && i != origin_index
            edge_vectors[:, j] .= simplex.vertices[:, i] .- origin_vertex
            j += 1
        end
    end

    # The vertex to be reflected
    vertex_to_reflect = simplex.vertices[:, facet_index]

    # Find normal to the hyperplane using QR decomposition
    # The normal is orthogonal to all edge vectors
    Q, _ = qr(edge_vectors)  # this will produce an InexactError if not full rank

    # Normal vector is the last column of Q (orthogonal to the span)
    normal = Q[:, end]

    # Ensure normal points away from the vertex to reflect
    direction_to_vertex = vertex_to_reflect - origin_vertex
    if dot(normal, direction_to_vertex) < 0
        normal = -normal
    end

    # Reflect the vertex across the hyperplane
    # Distance from vertex to hyperplane
    distance = dot(vertex_to_reflect - origin_vertex, normal)

    # Reflected vertex
    reflected_vertex = vertex_to_reflect - 2 * distance * normal

    # Create new simplex with reflected vertex
    new_simplex = Simplex(simplex)
    if T <: Integer
        new_simplex.vertices[:, facet_index] .= round.(T, reflected_vertex)
    else
        new_simplex.vertices[:, facet_index] .= reflected_vertex
    end

    return new_simplex
end

"""
    is_freudenthal(simplex::Simplex{<: Integer})

Check if a simplex follows the structure of a Freudenthal triangulation.
A simplex is Freudenthal if it can be obtained by permuting the coordinates
of the initial Freudenthal simplex and then translating by an integer vector.
"""
function is_freudenthal(simplex::Simplex{T}) where {T <: Integer}
    if simplex.n != simplex.dims
        return false
    end

    # Translate simplex so that the first vertex is at origin
    origin = simplex.vertices[:, 1]
    translated_vertices = simplex.vertices[:, 2:end] .- origin

    # Check if the translated simplex has the staircase pattern
    # when coordinates are sorted appropriately
    n = simplex.n

    # The translated vertices should form a pattern where each vertex
    # has exactly one more coordinate equal to 1 than the previous
    for i in 1:n
        vertex = translated_vertices[:, i]
        expected_ones = i

        # Count how many coordinates are 1 and how many are 0
        ones_count = 0
        zeros_count = 0

        for coord in vertex
            if coord == 1
                ones_count += 1
            elseif coord == 0
                zeros_count += 1
            else
                return false
            end
        end

        if ones_count != expected_ones || zeros_count != (n - expected_ones)
            return false
        end
    end

    # Additional check: verify the vertices form a valid simplex structure
    # The difference between consecutive vertices should be a unit vector
    for i in 2:n
        diff = translated_vertices[:, i] - translated_vertices[:, i - 1]
        unit_vector_count = 0
        for coord in diff
            if coord == 1
                unit_vector_count += 1
            elseif coord != 0
                return false
            end
        end
        if unit_vector_count != 1
            return false
        end
    end

    return true
end


"""
    freudenthal_initial_simplex([T=Int], n)

Generate the initial simplex of a Freudenthal triangulation in n dimensions. The
vertices form a staircase pattern: (0,0,…,0), (1,0,…,0), (1,1,0,…,0), …,
(1,1,…,1)
"""
function freudenthal_initial_simplex(T, n)
    if n < 1
        throw(ArgumentError("Dimension must be at least 1"))
    end
    if !(T <: Integer)
        throw(ArgumentError("Type T must be a subtype of Integer"))
    end
    simplex = Simplex{T}(undef, n, n)
    for j in 1:(n + 1)
        for i in 1:n
            simplex.vertices[i, j] = (i < j) ? one(T) : zero(T)
        end
    end
    return simplex
end
freudenthal_initial_simplex(n) = freudenthal_initial_simplex(Int, n)
"""
    freudenthal_reflect(simplex::Simplex{<:Integer}, facet_index)

Reflect a simplex in a Freudenthal triangulation around one of its facets. Note
that this preserves the Freudenthal structure and so is not a general
reflection. Returns a new simplex.
"""
function freudenthal_reflect(simplex::Simplex{T}, facet_index) where {T <: Integer}
    if (facet_index < 1) || (facet_index > simplex.n + 1)
        throw(ArgumentError("Facet index must be between 1 and $(simplex.n + 1)"))
    end
    if simplex.n != simplex.dims
        throw(ArgumentError("Reflection only implemented for full-dimensional simplices"))
    end

    # For Freudenthal triangulation, reflection across a facet corresponds to
    # swapping the order of coordinate increments in the staircase pattern

    # Find which coordinate is incremented at each step
    increment_order = zeros(Int, simplex.n)
    for i in 2:(simplex.n + 1)
        for j in 1:simplex.dims
            diff = simplex.vertices[j, i] - simplex.vertices[j, i - 1]
            if diff == 1
                increment_order[i - 1] = j
                break
            end
        end
        if increment_order[i - 1] == 0
            throw(ArgumentError("Invalid Freudenthal simplex structure"))
        end
    end

    # Create new increment order by swapping the position corresponding to facet_index
    if facet_index == 1
        # Reflecting first vertex: swap first two increments
        if simplex.n >= 2
            increment_order[1], increment_order[2] = increment_order[2], increment_order[1]
        end
    elseif facet_index == simplex.n + 1
        # Reflecting last vertex: swap last two increments
        if simplex.n >= 2
            n = simplex.n
            increment_order[n - 1], increment_order[n] = increment_order[n], increment_order[n - 1]
        end
    else
        # Reflecting middle vertex: swap the increment that creates this vertex with the next one
        idx = facet_index - 1  # Convert from vertex index to increment index
        increment_order[idx], increment_order[idx + 1] = increment_order[idx + 1], increment_order[idx]
    end

    # Construct new simplex from the new increment order
    new_simplex = Simplex(simplex)
    current_vertex = simplex.vertices[:, 1]  # slicing creates a copy
    for (i, inc) in enumerate(increment_order)
        current_vertex[inc] += 1
        new_simplex.vertices[:, i + 1] .= current_vertex
    end

    return new_simplex
end

end  # module
