module SimplicialContinuation

using FixedSizeArrays: FixedSizeArray, FixedSizeMatrixDefault, FixedSizeVector, FixedSizeVectorDefault
using LinearAlgebra: qr, qr!, nullspace, dot

const FSVector{T} = FixedSizeVectorDefault{T}

export Simplex
export simplex_dimension, space_dimension, reflect, is_freudenthal, freudenthal_initial_simplex, freudenthal_reflect
export get_facet, barycenter
export continuation, ContinuationPath

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
    simplex = Simplex{T}(undef, n, dims)
    for (i, vertex) in enumerate(vertices)
        if length(vertex) != dims
            throw(ArgumentError("All vertices must have the same dimension"))
        else
            simplex.vertices[:, i] .= vertex
        end
    end
    return simplex
end

function Simplex(vertices::AbstractMatrix{T}) where {T}
    dims = size(vertices, 1)
    n = size(vertices, 2)
    simplex = Simplex{T}(undef, n - 1, dims)
    copyto!(simplex.vertices, vertices)
    return simplex
end

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

_isapprox(a::T, b) where {T} = isapprox(a, b; atol = sqrt(eps(T)), rtol = sqrt(eps(T)))
_isapprox(a::Integer, b) = (a == b)

"""
    is_freudenthal(simplex::Simplex{T}; [isapprox])

Check if a simplex follows the structure of a Freudenthal triangulation.
`isapprox` can be specified for non-integer simplicies (defaults to
`Base.isapprox` with `atol=rtol=sqrt(eps)`.)
"""
function is_freudenthal(simplex::Simplex; isapprox = _isapprox)
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
            if isapprox(coord, 1)
                ones_count += 1
            elseif isapprox(coord, 0)
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
            if isapprox(coord, 1)
                unit_vector_count += 1
            elseif !isapprox(coord, 0)
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
    freudenthal_reflect(simplex::Simplex, facet_index)

Reflect a simplex in a Freudenthal triangulation around one of its facets. Note
that this preserves the Freudenthal structure and so is not a general
reflection. Returns a new simplex.
"""
function freudenthal_reflect(simplex::Simplex{T}, facet_index) where {T}
    if facet_index < 1 || facet_index > simplex.n + 1
        throw(ArgumentError("Facet index must be between 1 and $(simplex.n + 1)"))
    end
    if simplex.n != simplex.dims
        throw(ArgumentError("Only full-dimensional simplices supported"))
    end
    new_simplex = Simplex{T}(undef, simplex.n, simplex.dims)
    if facet_index == 1
        @views copyto!(new_simplex.vertices[:, 1:(end - 1)], simplex.vertices[:, 2:end])
        new_simplex.vertices[:, end] .= simplex.vertices[:, end] .- simplex.vertices[:, facet_index] .+ simplex.vertices[:, facet_index + 1]
    elseif facet_index == simplex.n + 1
        @views copyto!(new_simplex.vertices[:, 2:end], simplex.vertices[:, 1:(end - 1)])
        new_simplex.vertices[:, 1] .= simplex.vertices[:, facet_index - 1] .- simplex.vertices[:, facet_index] .+ simplex.vertices[:, 1]
    else
        copyto!(new_simplex.vertices, simplex.vertices)
        new_simplex.vertices[:, facet_index] .= simplex.vertices[:, facet_index - 1] .- simplex.vertices[:, facet_index] .+ simplex.vertices[:, facet_index + 1]
    end
    return new_simplex
end

# ── Phase 1: Geometric utilities ─────────────────────────────────────────────

"""
    get_facet(simplex, facet_index)

Return the facet of `simplex` opposite vertex `facet_index` as an `(n-1)`-simplex.
The returned simplex has the same ambient dimension but one fewer vertex.
"""
function get_facet(simplex::Simplex{T}, facet_index) where {T}
    if facet_index < 1 || facet_index > simplex.n + 1
        throw(ArgumentError("Facet index must be between 1 and $(simplex.n + 1)"))
    end
    facet = Simplex{T}(undef, simplex.n - 1, simplex.dims)
    j = 1
    for i in 1:(simplex.n + 1)
        if i != facet_index
            facet.vertices[:, j] .= simplex.vertices[:, i]
            j += 1
        end
    end
    return facet
end

"""
    barycenter(simplex)

Return the barycenter (centroid) of `simplex` as a `Vector`.
"""
function barycenter(simplex::Simplex{T}) where {T}
    result = zeros(T, simplex.dims)
    for i in 1:(simplex.n + 1)
        result .+= simplex.vertices[:, i]
    end
    return result ./ (simplex.n + 1)
end

# ── Phase 2: Labelled simplex ─────────────────────────────────────────────────

# Internal: pairs a Simplex with residuals f(v_k, p) at each vertex.
# `residuals` has shape (codim) × (n+1) where codim = dims - 1 for full-dimensional use.
struct LabeledSimplex{T, R}
    simplex::Simplex{T}
    residuals::FixedSizeMatrixDefault{R}
end

function label_simplex(simplex::Simplex{T}, f, p) where {T}
    r1 = f(view(simplex.vertices, :, 1), p)
    R = eltype(r1)
    m = length(r1)
    n_verts = simplex.n + 1
    residuals = FixedSizeMatrixDefault{R}(undef, m, n_verts)
    residuals[:, 1] .= r1
    for i in 2:n_verts
        residuals[:, i] .= f(view(simplex.vertices, :, i), p)
    end
    return LabeledSimplex{T, R}(simplex, residuals)
end

# ── Phase 3: Transversality ───────────────────────────────────────────────────

# Build and solve the m×m transversality system for a facet.
# `facet_residuals` is (m-1)×m: residuals at the m vertices of the facet.
# Returns (β, is_transversal).
function _facet_transversality(facet_residuals::AbstractMatrix)
    R = eltype(facet_residuals)
    m = size(facet_residuals, 2)
    A = Matrix{R}(undef, m, m)
    A[1:m-1, :] .= facet_residuals
    for j in 1:m
        A[m, j] = one(R)
    end
    b = zeros(R, m)
    b[m] = one(R)
    β = try
        A \ b
    catch
        return zeros(R, m), false
    end
    tol = 1e-8
    return β, all(x -> isfinite(x) && -tol ≤ x ≤ 1 + tol, β)
end

# Check transversality of the facet OPPOSITE vertex `facet_index` in a LabeledSimplex.
function check_transversality(ls::LabeledSimplex{T, R}, facet_index) where {T, R}
    n = ls.simplex.n
    m_res = size(ls.residuals, 1)
    facet_residuals = Matrix{R}(undef, m_res, n)
    j = 1
    for i in 1:(n + 1)
        if i != facet_index
            facet_residuals[:, j] .= view(ls.residuals, :, i)
            j += 1
        end
    end
    return _facet_transversality(facet_residuals)
end

# Find the exit facet (first transversal facet that is not the entry facet).
# Returns (facet_index, β) or (nothing, nothing).
function find_exit_facet(ls::LabeledSimplex, entry_facet_index)
    for i in 1:(ls.simplex.n + 1)
        i == entry_facet_index && continue
        β, ok = check_transversality(ls, i)
        ok && return i, β
    end
    return nothing, nothing
end

# Find any transversal facet (used for initialisation).
function _find_any_transversal_facet(ls::LabeledSimplex)
    for i in 1:(ls.simplex.n + 1)
        β, ok = check_transversality(ls, i)
        ok && return i, β
    end
    return nothing, nothing
end

# Compute the barycentric interpolation of the zero crossing on the given facet.
function facet_zero_point(ls::LabeledSimplex{T}, facet_index, β) where {T}
    result = zeros(T, ls.simplex.dims)
    j = 1
    for i in 1:(ls.simplex.n + 1)
        if i != facet_index
            result .+= β[j] .* view(ls.simplex.vertices, :, i)
            j += 1
        end
    end
    return result
end

# ── Phase 4: Pivot ────────────────────────────────────────────────────────────

# After freudenthal_reflect(simplex, exit_facet_index), the entry facet in the new
# simplex has the same index when exit_facet_index is interior, but shifts for the
# boundary cases because vertices are renumbered.
function _pivot_entry_facet(exit_facet_index, n)
    exit_facet_index == 1     && return n + 1
    exit_facet_index == n + 1 && return 1
    return exit_facet_index
end

# Pivot across the exit facet: reflect geometry, transfer residuals, evaluate new vertex.
# Returns (new_labeled_simplex, entry_facet_in_new_simplex).
function pivot(ls::LabeledSimplex{T, R}, exit_facet_index, f, p) where {T, R}
    n = ls.simplex.n
    new_simplex = freudenthal_reflect(ls.simplex, exit_facet_index)
    entry_facet = _pivot_entry_facet(exit_facet_index, n)

    new_residuals = FixedSizeMatrixDefault{R}(undef, size(ls.residuals, 1), n + 1)

    if exit_facet_index == 1
        # Old vertices 2..n+1 → new positions 1..n; new vertex at position n+1.
        copyto!(view(new_residuals, :, 1:n), view(ls.residuals, :, 2:(n + 1)))
        new_residuals[:, n + 1] .= f(view(new_simplex.vertices, :, n + 1), p)
    elseif exit_facet_index == n + 1
        # Old vertices 1..n → new positions 2..n+1; new vertex at position 1.
        copyto!(view(new_residuals, :, 2:(n + 1)), view(ls.residuals, :, 1:n))
        new_residuals[:, 1] .= f(view(new_simplex.vertices, :, 1), p)
    else
        # Only vertex at exit_facet_index changes; all others stay in place.
        copyto!(new_residuals, ls.residuals)
        new_residuals[:, exit_facet_index] .= f(view(new_simplex.vertices, :, exit_facet_index), p)
    end

    return LabeledSimplex{T, R}(new_simplex, new_residuals), entry_facet
end

# ── Phase 5: Initialisation ───────────────────────────────────────────────────

_grain_vec(grain::Number, n, ::Type{T}) where {T} = fill(T(grain), n)
_grain_vec(grain, n, ::Type{T}) where {T} = T.(grain)

# Find an initial LabeledSimplex transversal to the zero set of f near y0.
# The simplex is a scaled and centred Freudenthal simplex.
# Returns (labeled_simplex, entry_facet_index).
function find_transverse(f, p, y0, grain)
    n = length(y0)
    T = float(eltype(y0))
    gv = _grain_vec(grain, n, T)

    unit_simplex = freudenthal_initial_simplex(T, n)
    bc = barycenter(unit_simplex)
    translation = T.(y0) .- bc .* gv

    scaled = Simplex{T}(undef, n, n)
    for i in 1:(n + 1)
        scaled.vertices[:, i] .= view(unit_simplex.vertices, :, i) .* gv .+ translation
    end

    ls = label_simplex(scaled, f, p)
    facet_idx, _ = _find_any_transversal_facet(ls)
    !isnothing(facet_idx) && return ls, facet_idx

    # Try immediate Freudenthal neighbours.
    for k in 1:(n + 1)
        neighbour = freudenthal_reflect(scaled, k)
        nls = label_simplex(neighbour, f, p)
        facet_idx, _ = _find_any_transversal_facet(nls)
        !isnothing(facet_idx) && return nls, facet_idx
    end

    throw(ArgumentError("Could not find a transversal simplex near y0; try adjusting grain."))
end

# ── Phase 6: Iterator API ─────────────────────────────────────────────────────

"""
    ContinuationPath{T, R, F, P}

Iterator over the zero curve of `f(x, p) = 0`. Each `iterate` call yields the
approximate coordinates of the zero crossing on the current exit facet and
advances to the next simplex. Create via `continuation`.
"""
struct ContinuationPath{T, R, F, P}
    f::F
    p::P
    initial_ls::LabeledSimplex{T, R}
    initial_entry_facet::Int
    maxsteps::Int
end

"""
    continuation(f, p, y0; grain=1.0, maxsteps=1000)

Return a `ContinuationPath` iterator that traces the zero curve of `f(x, p) = 0`
starting near `y0`. `f` must have the signature `f(x, p) -> res` (SciML convention).
`grain` sets the simplex step size; smaller values give higher spatial resolution.
`grain` may be a scalar (uniform scaling) or a vector/tuple of length `n` for
per-dimension scaling, useful when the state variables have different natural scales.
"""
function continuation(f, p, y0; grain = 1.0, maxsteps = 1000)
    ls, entry_facet = find_transverse(f, p, y0, grain)
    T = eltype(ls.simplex)
    R = eltype(ls.residuals)
    return ContinuationPath{T, R, typeof(f), typeof(p)}(f, p, ls, entry_facet, maxsteps)
end

Base.IteratorSize(::Type{<:ContinuationPath}) = Base.SizeUnknown()
Base.eltype(::Type{ContinuationPath{T, R, F, P}}) where {T, R, F, P} = Vector{T}

function Base.iterate(path::ContinuationPath)
    ls = path.initial_ls
    entry = path.initial_entry_facet
    exit_facet, β = find_exit_facet(ls, entry)
    isnothing(exit_facet) && return nothing
    point = facet_zero_point(ls, exit_facet, β)
    new_ls, new_entry = pivot(ls, exit_facet, path.f, path.p)
    return point, (new_ls, new_entry, 1)
end

function Base.iterate(path::ContinuationPath, state)
    ls, entry, step = state
    step ≥ path.maxsteps && return nothing
    exit_facet, β = find_exit_facet(ls, entry)
    isnothing(exit_facet) && return nothing
    point = facet_zero_point(ls, exit_facet, β)
    new_ls, new_entry = pivot(ls, exit_facet, path.f, path.p)
    return point, (new_ls, new_entry, step + 1)
end

end  # module
