module SimplicialContinuation

using FixedSizeArrays: FixedSizeArray, FixedSizeMatrixDefault, FixedSizeVector, FixedSizeVectorDefault
using LinearAlgebra: qr, qr!, nullspace, dot, lu, norm, I
using Random: default_rng, AbstractRNG, randn
using Statistics: mean, var

const FSVector{T} = FixedSizeVectorDefault{T}

export Simplex
export simplex_dimension, space_dimension, reflect, is_freudenthal, freudenthal_initial_simplex, freudenthal_reflect
export get_facet, barycenter
export continuation, ContinuationPath
export check_start_point, check_scaling, score_start_point, measure_noise, improve_start_point, curr_fom

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

# Evaluate `g(vertex, p)` at every vertex of `simplex`, returning a
# `length(g(...))` × `(n+1)` matrix. Shared by residual labelling and
# auxiliary-function interpolation (see `continuation`'s `aux` keyword).
function _label_vertices(simplex::Simplex, g, p)
    v1 = g(view(simplex.vertices, :, 1), p)
    G = eltype(v1)
    m = length(v1)
    n_verts = simplex.n + 1
    values = FixedSizeMatrixDefault{G}(undef, m, n_verts)
    values[:, 1] .= v1
    for i in 2:n_verts
        values[:, i] .= g(view(simplex.vertices, :, i), p)
    end
    return values
end

function label_simplex(simplex::Simplex{T}, f, p) where {T}
    residuals = _label_vertices(simplex, f, p)
    return LabeledSimplex{T, eltype(residuals)}(simplex, residuals)
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

# Return `M` with column `idx` removed (dense copy). Used to extract a
# facet's vertex coordinates or residual labels from the owning simplex's
# full data matrix.
function _drop_column(M::AbstractMatrix{R}, idx) where {R}
    n = size(M, 2)
    out = Matrix{R}(undef, size(M, 1), n - 1)
    j = 1
    for i in 1:n
        if i != idx
            out[:, j] .= view(M, :, i)
            j += 1
        end
    end
    return out
end

# Check transversality of the facet OPPOSITE vertex `facet_index` in a LabeledSimplex.
function check_transversality(ls::LabeledSimplex, facet_index)
    return _facet_transversality(_drop_column(ls.residuals, facet_index))
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

# Barycentric combination sum_j β[j] * data[:, i_j], over the columns of
# `data` other than `facet_index`. Shared by the zero-crossing point
# (vertex coordinates) and auxiliary-function interpolation (aux values).
function _facet_combination(data::AbstractMatrix{R}, facet_index, β) where {R}
    result = zeros(R, size(data, 1))
    j = 1
    for i in 1:size(data, 2)
        if i != facet_index
            result .+= β[j] .* view(data, :, i)
            j += 1
        end
    end
    return result
end

# Compute the barycentric interpolation of the zero crossing on the given facet.
facet_zero_point(ls::LabeledSimplex, facet_index, β) = _facet_combination(ls.simplex.vertices, facet_index, β)

# Refine β by up to `steps` modified-Newton corrections against the true
# residual f, evaluated at the barycentric point on the facet's own vertex
# coordinates. The facet's transversality system M (residuals stacked on a
# row of ones) is factored once and reused for every step, since it is only
# a secant approximation of the Jacobian on the facet — refactoring would
# not improve it. Falls back to the last accepted β as soon as a step fails
# to reduce ‖f‖ (the noise/curvature floor) or would leave the facet
# (β ∉ [0,1]^m); each step costs one extra evaluation of f.
function _corrector_refine(facet_vertices::AbstractMatrix, facet_residuals::AbstractMatrix{R}, β, f, p, steps) where {R}
    steps ≤ 0 && return β
    m = length(β)
    A = Matrix{R}(undef, m, m)
    A[1:(m - 1), :] .= facet_residuals
    A[m, :] .= one(R)
    fact = lu(A)
    tol = 1.0e-8
    rn_prev = typemax(R)
    for _ in 1:steps
        y1 = facet_vertices * β
        r = f(y1, p)
        rn = norm(r)
        rn ≥ rn_prev && break
        rn_prev = rn
        Δβ = fact \ vcat(-r, zero(R))
        β_trial = β .+ Δβ
        all(x -> -tol ≤ x ≤ 1 + tol, β_trial) || break
        β = β_trial
    end
    return β
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

# Shift the per-vertex data matrix `data` ([m x (n+1)]) the way pivoting
# renumbers simplex vertices when reflecting across facet `exit_facet_index`
# of an n-simplex, and fill the new vertex's column by evaluating `g` at the
# corresponding vertex of `new_simplex`. Shared by residual and
# auxiliary-function pivoting.
function _pivot_data(data::AbstractMatrix{R}, new_simplex::Simplex, exit_facet_index, n, g, p) where {R}
    new_data = FixedSizeMatrixDefault{R}(undef, size(data, 1), n + 1)
    if exit_facet_index == 1
        # Old vertices 2..n+1 → new positions 1..n; new vertex at position n+1.
        copyto!(view(new_data, :, 1:n), view(data, :, 2:(n + 1)))
        new_data[:, n + 1] .= g(view(new_simplex.vertices, :, n + 1), p)
    elseif exit_facet_index == n + 1
        # Old vertices 1..n → new positions 2..n+1; new vertex at position 1.
        copyto!(view(new_data, :, 2:(n + 1)), view(data, :, 1:n))
        new_data[:, 1] .= g(view(new_simplex.vertices, :, 1), p)
    else
        # Only vertex at exit_facet_index changes; all others stay in place.
        copyto!(new_data, data)
        new_data[:, exit_facet_index] .= g(view(new_simplex.vertices, :, exit_facet_index), p)
    end
    return new_data
end

# Pivot across the exit facet: reflect geometry, transfer residuals, evaluate new vertex.
# Returns (new_labeled_simplex, entry_facet_in_new_simplex).
function pivot(ls::LabeledSimplex{T, R}, exit_facet_index, f, p) where {T, R}
    n = ls.simplex.n
    new_simplex = freudenthal_reflect(ls.simplex, exit_facet_index)
    entry_facet = _pivot_entry_facet(exit_facet_index, n)
    new_residuals = _pivot_data(ls.residuals, new_simplex, exit_facet_index, n, f, p)
    return LabeledSimplex{T, R}(new_simplex, new_residuals), entry_facet
end

# ── Phase 5: Initialisation ───────────────────────────────────────────────────

_grain_vec(grain::Number, n, ::Type{T}) where {T} = fill(T(grain), n)
_grain_vec(grain, n, ::Type{T}) where {T} = T.(grain)

# Find an initial LabeledSimplex transversal to the zero set of f near y0.
# The simplex is a scaled and centred Freudenthal simplex. `aux`, if not
# `nothing`, is evaluated at the same vertices as f (see `continuation`'s
# `aux` keyword). Returns (labeled_simplex, aux_values_or_nothing, entry_facet_index).
function find_transverse(f, p, y0, grain, aux = nothing)
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
    if !isnothing(facet_idx)
        aux_vals = isnothing(aux) ? nothing : _label_vertices(scaled, aux, p)
        return ls, aux_vals, facet_idx
    end

    # Try immediate Freudenthal neighbours.
    for k in 1:(n + 1)
        neighbour = freudenthal_reflect(scaled, k)
        nls = label_simplex(neighbour, f, p)
        facet_idx, _ = _find_any_transversal_facet(nls)
        if !isnothing(facet_idx)
            aux_vals = isnothing(aux) ? nothing : _label_vertices(neighbour, aux, p)
            return nls, aux_vals, facet_idx
        end
    end

    throw(ArgumentError("Could not find a transversal simplex near y0; try adjusting grain."))
end

# ── Phase 6: Iterator API ─────────────────────────────────────────────────────

"""
    ContinuationPath{T, R, F, P, AUX, AR}

Iterator over the zero curve of `f(x, p) = 0`. Each `iterate` call yields the
approximate coordinates of the zero crossing on the current exit facet (or,
when an `aux` function was supplied, a `(point, aux_value)` tuple) and
advances to the next simplex. Create via `continuation`.
"""
struct ContinuationPath{T, R, F, P, AUX, AR}
    f::F
    p::P
    aux::AUX
    initial_ls::LabeledSimplex{T, R}
    initial_aux_values::AR
    initial_entry_facet::Int
    maxsteps::Int
    corrector_steps::Int
end

"""
    continuation(f, p, y0; grain=1.0, maxsteps=1000, corrector_steps=0, aux=nothing)

Return a `ContinuationPath` iterator that traces the zero curve of `f(x, p) = 0`
starting near `y0`. `f` must have the signature `f(x, p) -> res` (SciML convention).
`grain` sets the simplex step size; smaller values give higher spatial resolution.
`grain` may be a scalar (uniform scaling) or a vector/tuple of length `n` for
per-dimension scaling, useful when the state variables have different natural scales.

`corrector_steps` applies up to that many modified-Newton corrections to each
exit point, refining it beyond the piecewise-linear approximation (each step
costs one extra evaluation of `f`; see `_corrector_refine`). It amplifies
measurement noise, so keep it at `0` (the default) for noisy/experimental
residuals unless the improved accuracy is worth the extra evaluations.

`aux(x, p) -> value` is an optional extra function, evaluated at the same
simplex vertices as `f` and interpolated onto each exit point the same way
`f`'s zero crossing is. Useful when a quantity of interest is nearly free to
obtain alongside the residual (e.g. from the same measurement or solve) and
cheaper to interpolate than to re-evaluate at the crossing point. When `aux`
is given, `iterate` yields `(point, aux_value)` tuples instead of bare
`point` vectors.
"""
function continuation(f, p, y0; grain = 1.0, maxsteps = 1000, corrector_steps = 0, aux = nothing)
    ls, aux_vals, entry_facet = find_transverse(f, p, y0, grain, aux)
    T = eltype(ls.simplex)
    R = eltype(ls.residuals)
    return ContinuationPath{T, R, typeof(f), typeof(p), typeof(aux), typeof(aux_vals)}(
        f, p, aux, ls, aux_vals, entry_facet, maxsteps, corrector_steps
    )
end

Base.IteratorSize(::Type{<:ContinuationPath}) = Base.SizeUnknown()
function Base.eltype(::Type{ContinuationPath{T, R, F, P, AUX, AR}}) where {T, R, F, P, AUX, AR}
    return AUX === Nothing ? Vector{T} : Tuple{Vector{T}, Vector{eltype(AR)}}
end

# Find the exit facet, optionally apply the Newton corrector, and return the
# (possibly corrected) zero-crossing point together with the exit facet
# index and barycentric coordinates β (needed to interpolate `aux`).
function _advance(ls::LabeledSimplex, entry, f, p, corrector_steps)
    exit_facet, β = find_exit_facet(ls, entry)
    isnothing(exit_facet) && return nothing
    if corrector_steps > 0
        facet_vertices = _drop_column(ls.simplex.vertices, exit_facet)
        facet_residuals = _drop_column(ls.residuals, exit_facet)
        β = _corrector_refine(facet_vertices, facet_residuals, β, f, p, corrector_steps)
    end
    point = facet_zero_point(ls, exit_facet, β)
    return point, exit_facet, β
end

function _advance_state(path::ContinuationPath, ls, aux_vals, entry)
    advanced = _advance(ls, entry, path.f, path.p, path.corrector_steps)
    isnothing(advanced) && return nothing
    point, exit_facet, β = advanced
    new_ls, new_entry = pivot(ls, exit_facet, path.f, path.p)
    if path.aux === nothing
        return point, (new_ls, nothing, new_entry)
    else
        u = _facet_combination(aux_vals, exit_facet, β)
        new_aux_vals = _pivot_data(aux_vals, new_ls.simplex, exit_facet, ls.simplex.n, path.aux, path.p)
        return (point, u), (new_ls, new_aux_vals, new_entry)
    end
end

function Base.iterate(path::ContinuationPath)
    result = _advance_state(path, path.initial_ls, path.initial_aux_values, path.initial_entry_facet)
    isnothing(result) && return nothing
    value, (new_ls, new_aux_vals, new_entry) = result
    return value, (new_ls, new_aux_vals, new_entry, 1)
end

function Base.iterate(path::ContinuationPath, state)
    ls, aux_vals, entry, step = state
    step ≥ path.maxsteps && return nothing
    result = _advance_state(path, ls, aux_vals, entry)
    isnothing(result) && return nothing
    value, (new_ls, new_aux_vals, new_entry) = result
    return value, (new_ls, new_aux_vals, new_entry, step + 1)
end

# ── Phase 7: Start-point diagnostics ─────────────────────────────────────────
#
# These are standalone diagnostics for choosing/validating `y0` and `grain`
# before (or independent of) calling `continuation`. They matter most for
# *experimental* continuation, where `f` is a noisy measurement rather than a
# clean numerical residual: see `measure_noise` and `score_start_point` in
# particular. Each function takes `f`, `p`, and a point/`grain` explicitly and
# returns a `NamedTuple`; none of them hold state, so e.g. `sens_ref` from
# `check_scaling` must be passed explicitly to `measure_noise`/`curr_fom` if
# its normalisation is wanted there.

"""
    check_start_point(f, p, y0; grain=1.0)

Diagnostic: evaluate `‖f(y0, p)‖` and compare it against `grain`. For
`continuation` to have a reasonable chance of finding a transversal simplex
near `y0`, this norm should be well below `grain` (as a rule of thumb,
`grain/100`, using the smallest per-dimension grain if `grain` is a vector).
If this fails, `improve_start_point` is unlikely to be enough on its own.

Returns `(; norm, threshold, ok)`.
"""
function check_start_point(f, p, y0; grain = 1.0)
    n = length(y0)
    T = float(eltype(y0))
    nrm = norm(f(y0, p))
    threshold = minimum(_grain_vec(grain, n, T)) / 100
    return (; norm = nrm, threshold, ok = nrm < threshold)
end

"""
    check_scaling(f, p, y0; grain=1.0)

Diagnostic: perturb each coordinate of `y0` in turn by `grain` (or its
per-dimension value) and report `‖f(y0 + grain*e_j, p) - f(y0, p)‖` for each
`j`. These sensitivities should be comparable in magnitude; a large max/min
ratio signals coordinates that are mis-scaled relative to `grain` and to
each other — rescale the low-sensitivity coordinate up (or the
high-sensitivity one down) in a wrapper around `f`.

Also returns `sens_ref`, the mean sensitivity across coordinates: the
reference scale used by `measure_noise` and `curr_fom` to normalise their
figures of merit into "residual change per grain-sized step" units.

Returns `(; sensitivities, sens_ref, ratio, idx_min, idx_max, status)` with
`status ∈ (:ok, :caution, :warning)` (thresholds at ratio 10 and 100).
"""
function check_scaling(f, p, y0; grain = 1.0)
    n = length(y0)
    T = float(eltype(y0))
    gv = _grain_vec(grain, n, T)
    y0T = T.(y0)
    r0 = f(y0T, p)
    sensitivities = Vector{float(eltype(r0))}(undef, n)
    for j in 1:n
        yj = copy(y0T)
        yj[j] += gv[j]
        sensitivities[j] = norm(f(yj, p) .- r0)
    end
    idx_max = argmax(sensitivities)
    idx_min = argmin(sensitivities)
    ratio = sensitivities[idx_min] > 0 ? sensitivities[idx_max] / sensitivities[idx_min] : Inf
    status = ratio > 100 ? :warning : (ratio > 10 ? :caution : :ok)
    return (; sensitivities, sens_ref = mean(sensitivities), ratio, idx_min, idx_max, status)
end

# Regular m-simplex (m+1 vertices) in R^m, as columns of an m × (m+1)
# matrix: project the standard basis of R^{m+1} onto the hyperplane
# orthogonal to (1,...,1). Simpler than porting the bespoke Gram-Schmidt
# construction in codim1.cc; any equilateral simplex works equally well for
# `improve_start_point`'s local secant Jacobian, regardless of orientation.
function _equilateral_simplex(::Type{T}, m) where {T}
    Qfull = qr(ones(T, m + 1, 1)).Q * Matrix{T}(I, m + 1, m + 1)
    P = Qfull[:, 2:end]  # (m+1) × m, orthonormal basis of the complement of (1,...,1)
    return Matrix(P')    # m × (m+1); columns are the simplex vertices
end

"""
    improve_start_point(f, p, y0; grain=1.0)

Diagnostic/corrector: move `y0` closer to the zero set of `f` by one Newton
step, using a local secant Jacobian built from an equilateral simplex of
radius `π * grain` centred on `y0`. Returns the improved point.

The *last* coordinate of `y0` is held fixed, following the convention that
it is a continuation/bifurcation parameter and the other coordinates are
the state being solved for. `grain` here is always a scalar (unlike
`continuation`'s per-dimension `grain`): the equilateral simplex used for
the Jacobian estimate has no preferred coordinate direction to scale
anisotropically. Calling `check_start_point`/`improve_start_point` more
than once can help, but if there is no improvement after the first call, do
not expect more.
"""
function improve_start_point(f, p, y0; grain = 1.0)
    n = length(y0)
    m = n - 1
    T = float(eltype(y0))
    y0T = T.(y0)

    verts = _equilateral_simplex(T, m)  # m × n
    facet = Matrix{T}(undef, n, n)
    facet[1:m, :] .= verts .* (T(pi) * grain)
    facet[n, :] .= zero(T)
    facet .+= y0T

    residuals = Matrix{T}(undef, m, n)
    for k in 1:n
        residuals[:, k] .= f(view(facet, :, k), p)
    end
    M = Matrix{T}(undef, n, n)
    M[1:m, :] .= residuals
    M[n, :] .= one(T)

    rhs = vcat(f(y0T, p), zero(T))
    dalpha = M \ rhs
    alpha = fill(one(T) / n, n) .- dalpha
    return facet * alpha
end

# A_n root-lattice simplex: n+1 vertices in R^n, centred at the origin
# (v_k = -e_k for k = 1:n, v_{n+1} = ones(n)). Used by `score_start_point`
# for an isotropic probe simplex — unlike the Freudenthal simplex used for
# stepping, it has no preferred axes to bias the random-rotation trials.
function _an_simplex(::Type{T}, n) where {T}
    V = zeros(T, n, n + 1)
    for k in 1:n
        V[k, k] = -one(T)
    end
    V[:, n + 1] .= one(T)
    return V
end

_random_orthogonal(::Type{T}, n, rng::AbstractRNG) where {T} = Matrix(qr(randn(rng, T, n, n)).Q)

"""
    score_start_point(f, p, y0, grain; ntrials=20, rng=Random.default_rng())

Diagnostic: `ntrials` times, build a randomly-rotated, isotropic simplex of
size `grain` centred (with a small random jitter, to avoid the rare
vertex-on-zero-curve degeneracy) on `y0`, and check whether the zero curve
threads it transversally (exactly two transversal facets). Returns
`(; score, hits, ntrials)` where `score = hits / ntrials`.

A score near 1 means `y0` and `grain` are consistent with the zero curve. A
low score means either `y0` needs improving (`improve_start_point`) or
`grain` is too small relative to the resolution/noise floor of `f`: in
simulation `grain` can be arbitrarily small, but for a measured/noisy `f`
the score can collapse well before `grain` reaches machine precision (see
`measure_noise` for where that floor is).
"""
function score_start_point(f, p, y0, grain; ntrials::Int = 20, rng::AbstractRNG = default_rng())
    n = length(y0)
    T = float(eltype(y0))
    gv = _grain_vec(grain, n, T)
    y0T = T.(y0)
    base = _an_simplex(T, n)

    hits = 0
    for _ in 1:ntrials
        Q = _random_orthogonal(T, n, rng)
        jitter = T(0.1) .* gv .* (rand(rng, T, n) .- T(0.5))
        simplex = Simplex{T}(undef, n, n)
        for k in 1:(n + 1)
            simplex.vertices[:, k] .= (Q * view(base, :, k)) .* gv .+ y0T .+ jitter
        end
        ls = label_simplex(simplex, f, p)
        n_trans = count(k -> check_transversality(ls, k)[2], 1:(n + 1))
        n_trans == 2 && (hits += 1)
    end
    return (; score = hits / ntrials, hits, ntrials)
end

"""
    measure_noise(f, p, y; n=10, sens_ref=nothing, grain=nothing)

Diagnostic: evaluate `f(y, p)` `n` times without moving `y`, and report the
spread. Returns `sigma = sqrt(sum_j var_j)`, the expected magnitude of the
noise on one residual evaluation — directly comparable to `‖f(y,p)‖` and to
`sens_ref` from `check_scaling`.

If `sens_ref` is given, also returns the Jacobian SNR (`sens_ref / sigma`;
the divided-difference Jacobian's signal is `sens_ref` and its noise is
`sigma`, and since `sens_ref` falls roughly linearly with `grain` while
`sigma` does not, SNR sets the smallest usable `grain` — aim for ~20, treat
below ~5 as unreliable) and the FOM floor (`sigma / sens_ref`, the smallest
value `curr_fom` can report, since no corrector can drive the residual
below the noise). If `grain` is also given, `suggested_grain` estimates the
grain needed to reach an SNR of 20.

Returns `(; sigma, snr, fom_floor, suggested_grain)`, with the last three
`missing` when `sens_ref` (and, for the last, `grain`) are not supplied.
"""
function measure_noise(f, p, y; n::Int = 10, sens_ref = nothing, grain = nothing)
    n ≥ 2 || throw(ArgumentError("n must be at least 2"))
    r1 = f(y, p)
    m = length(r1)
    samples = Matrix{float(eltype(r1))}(undef, m, n)
    samples[:, 1] .= r1
    for k in 2:n
        samples[:, k] .= f(y, p)
    end
    sigma = sqrt(sum(var(view(samples, j, :)) for j in 1:m))

    sens_ref === nothing && return (; sigma, snr = missing, fom_floor = missing, suggested_grain = missing)
    snr = sigma > 0 ? sens_ref / sigma : Inf
    fom_floor = sigma > 0 ? sigma / sens_ref : zero(sigma)
    suggested_grain = (grain === nothing || sigma <= 0) ? missing : grain * 20 / snr
    return (; sigma, snr, fom_floor, suggested_grain)
end

"""
    curr_fom(f, p, y; sens_ref=nothing)

Figure of merit for the current point: `‖f(y, p)‖ / sens_ref` (or the raw
norm if `sens_ref` is not given), expressed in units of the residual change
produced by one grain-sized coordinate step (see `check_scaling`).

The PL solution is an exact zero only of the piecewise-linear interpolant
of `f` over the current simplex, never of `f` itself, so this is O(1) at
best and never exactly zero. What matters is its trend along the curve:
steady and O(1) or below is healthy; a slow climb means curvature is
outrunning `grain`; a jump of orders of magnitude means the tracked
solution has been lost (branch jump, instrument fault, ...). It is
necessary but not sufficient: a small value confirms *a* zero is being
tracked, not that it is the intended one.

Returns `(; fom, residual)`.
"""
function curr_fom(f, p, y; sens_ref = nothing)
    r = f(y, p)
    rn = norm(r)
    fom = sens_ref === nothing ? rn : rn / sens_ref
    return (; fom, residual = r)
end

end  # module
