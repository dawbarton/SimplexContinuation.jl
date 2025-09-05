module SimplexContinuation

using FixedSizeArrays: FixedSizeArray, FixedSizeVector, FixedSizeVectorDefault
using LinearAlgebra: qr, qr!, nullspace, dot

const FSVector{T} = FixedSizeVectorDefault{T}

export Simplex
export simplex_dimension, space_dimension, reflect, freudenthal_initial_simplex, freudenthal_reflect

"""
    Simplex{T}(vertices)

A representation of an `n`-simplex that lives in dimension `dims` which may or
may not be the same as `n`. Both `n` and `dims` are inferred from the size of
`vertices` and each corresponding vertex. If not specifed, `T` will be inferred.

Vertices are stored in lexicographical order by default. However, since the
vectors stored are mutable, the lexicographical order may be violated. To
restore the order, call `sort!` on the simplex.
"""
struct Simplex{T}
    n::Int
    dims::Int
    vertices::FSVector{FSVector{T}}

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
        new_vertices = FSVector{FSVector{T}}(undef, n + 1)
        for (i, vertex) in enumerate(vertices)
            if length(vertex) != dims
                throw(ArgumentError("All vertices must have the same dimension"))
            else
                new_vertices[i] = FSVector{T}(vertex)
            end
        end
        return new{T}(n, dims, sort!(new_vertices))
    end
    function Simplex{T}(::UndefInitializer, n::Integer, dims::Integer = n) where {T}
        new_vertices = FSVector{FSVector{T}}(undef, n + 1)
        for i in eachindex(new_vertices)
            new_vertices[i] = FSVector{T}(undef, dims)
        end
        return new{T}(n, dims, new_vertices)
    end
end

function Simplex{NT}(simplex::Simplex{T}) where {NT, T}
    new_simplex = Simplex{NT}(undef, simplex.n, simplex.dims)
    for i in eachindex(simplex.vertices)
        new_simplex.vertices[i] .= simplex.vertices[i]
    end
    return new_simplex
end
(Simplex(simplex::Simplex{T}) where {T}) = Simplex{T}(simplex)

Base.sort!(simplex::Simplex) = (sort!(simplex.vertices); return simplex)
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
    origin_vertex = simplex.vertices[origin_index]
    edge_vectors = FixedSizeArray{T}(undef, (simplex.dims, simplex.n - 1))
    j = 1
    for i in (origin_index + 1):length(simplex.vertices)
        if i != facet_index
            edge_vectors[:, j] .= simplex.vertices[i] .- origin_vertex
            j += 1
        end
    end

    # The vertex to be reflected
    vertex_to_reflect = simplex.vertices[facet_index]

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
        new_simplex.vertices[facet_index] .= round.(T, reflected_vertex)
    else
        new_simplex.vertices[facet_index] .= reflected_vertex
    end

    return sort!(new_simplex)
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
    simplex = Simplex{T}(undef, n, n)
    for (j, vertex) in enumerate(simplex.vertices)
        for i in eachindex(vertex)
            vertex[i] = (i < j) ? one(T) : zero(T)
        end
    end
    return simplex
end
freudenthal_initial_simplex(n) = freudenthal_initial_simplex(Int, n)

"""
    freudenthal_reflect(simplex::Simplex, facet_index)

Reflect a simplex in a Freudenthal triangulation around one of its facets. Note
that this preserves the Freudenthal structure and so is not a general
reflection. Returns a new simplex.
"""
function freudenthal_reflect(simplex::Simplex, facet_index)
    if (facet_index < 1) || (facet_index > simplex.n + 1)
        throw(ArgumentError("Facet index must be between 1 and $(simplex.n + 1)"))
    end
    if simplex.n != simplex.dims
        throw(ArgumentError("Reflection only implemented for full-dimensional simplices"))
    end

    new_simplex = Simplex(simplex)
    vertices = simplex.vertices
    new_vertices = new_simplex.vertices
    if facet_index == 1
        new_vertices[facet_index] .= vertices[end] .- vertices[facet_index] .+ vertices[facet_index + 1]
    elseif facet_index == simplex.n + 1
        new_vertices[facet_index] .= vertices[facet_index - 1] .- vertices[facet_index] .+ vertices[1]
    else
        new_vertices[facet_index] .= vertices[facet_index - 1] .- vertices[facet_index] .+ vertices[facet_index + 1]
    end
    return sort!(new_simplex)
end

end  # module
