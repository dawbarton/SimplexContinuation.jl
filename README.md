# SimplicialContinuation.jl

[![Build Status](https://github.com/dawbarton/SimplicialContinuation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dawbarton/SimplicialContinuation.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/dawbarton/SimplicialContinuation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dawbarton/SimplicialContinuation.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Derivative-free, simplex-based codimension-1 continuation. Traces the zero set of a smooth map `f : Rⁿ → Rⁿ⁻¹` by marching through a Freudenthal triangulation without computing any derivatives or solving nonlinear systems. The underlying piecewise-linear continuation method originates with:

> Allgower, E. L. and Georg, K. (1980). Simplicial and Continuation Methods for Approximating Fixed Points and Solutions to Systems of Equations. *SIAM Review*, 22(1), 28–85.
>
> Allgower, E. L. and Georg, K. (2003). *Introduction to Numerical Continuation Methods.* SIAM Classics in Applied Mathematics 45.

This package's specific approach follows the simplex-marching formulation described in:

> Henderson, M. E. and Melville, R. (2023). *Piecewise Linear Continuation: Derivative-free Manifold Generation.*

## Overview

Continuation methods trace a curve (or more generally a manifold) defined implicitly as the zero set of a system of equations. SimplicialContinuation.jl does this without derivatives by tiling space with simplices and tracking which facet of each simplex the zero curve passes through.

The map `f` follows the SciML convention:

```julia
f(x, p) -> residual
```

- `x` is a vector containing **all** variables involved in the continuation — both state variables and any parameters being varied. From the perspective of `f`, there is no distinction.
- `p` is an optional collection of **fixed** parameters that are passed through unchanged. Use `nothing` if there are no fixed parameters.
- The return value is a vector of length `n - 1` (one fewer than `length(x)`), so that the zero set is a 1-manifold (curve).

## Quick start

```julia
using SimplicialContinuation

# Trace the unit circle: x₁² + x₂² = 1
f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
y0 = [1.0, 0.0]   # a known point on the curve

path = continuation(f, nothing, y0; grain = 0.2, maxsteps = 200)

for pt in path
    println(pt)
end
```

`continuation` returns a lazy iterator; points are computed on demand. Use standard Julia iteration tools to consume it:

```julia
# Collect all points (up to maxsteps)
points = collect(path)

# Take only the first 50
points = collect(Iterators.take(path, 50))
```

## Interface

### `continuation(f, p, y0; grain, maxsteps)`

| Argument | Description |
|---|---|
| `f` | Residual function `f(x, p) -> res` with `length(res) == length(x) - 1` |
| `p` | Fixed parameters passed to `f` unchanged; use `nothing` if not needed |
| `y0` | Initial point near the zero curve (`Vector`, any numeric type) |
| `grain` | Step size (scalar or vector; see below) |
| `maxsteps` | Maximum number of continuation steps (default `1000`) |

Returns a `ContinuationPath` iterator whose elements are `Vector{T}` points lying approximately on the zero curve.

### The `grain` parameter

`grain` controls the size of the simplices used to tile space.

- **Scalar** — uniform step size in every dimension: `grain = 0.2`
- **Vector or tuple** — independent step size per dimension: `grain = (0.5, 0.1, 0.1)`. Use this when the continuation variables have different natural scales.

Smaller grain gives points closer together and smaller approximation error, at the cost of more steps to cover the same arc length.

## Examples

### Circle (2D)

```julia
f(x, _) = [x[1]^2 + x[2]^2 - 1.0]
path = continuation(f, nothing, [1.0, 0.0]; grain = 0.2, maxsteps = 200)
points = collect(path)
```

### Fixed parameters

When part of `f` depends on a parameter that is not being continued, pass it via `p`:

```julia
# Trace a circle of given radius; radius is a fixed parameter, not a continuation variable
f(x, r) = [x[1]^2 + x[2]^2 - r]
path = continuation(f, 4.0, [2.0, 0.0]; grain = 0.3, maxsteps = 200)
```

### Curve in 3D

`f : R³ → R²` gives a curve in 3D — the intersection of two surfaces:

```julia
# Sphere x₁² + x₂² + x₃² = 2 intersected with cylinder x₁² + x₂² = 1
f(x, _) = [x[1]^2 + x[2]^2 + x[3]^2 - 2.0,
           x[1]^2 + x[2]^2         - 1.0]
path = continuation(f, nothing, [1.0, 0.0, 1.0]; grain = 0.15, maxsteps = 300)
```

### Per-dimension grain

```julia
# Ellipse (x/2)² + y² = 1 — x spans [-2,2], y spans [-1,1]
f(x, _) = [(x[1]/2)^2 + x[2]^2 - 1.0]
path = continuation(f, nothing, [2.0, 0.0]; grain = (0.4, 0.2), maxsteps = 200)
```

## Accuracy

Points on the returned path lie on simplex facets and satisfy `f(x) ≈ 0` only up to the piecewise-linear approximation error, which scales as O(`grain`). Reduce `grain` for higher accuracy.

## Algorithm

The method triangulates Rⁿ using the Freudenthal triangulation — a regular tiling by n-simplices with a staircase vertex structure. Starting from an initial simplex near `y0`, the algorithm:

1. Finds the facet through which the zero curve enters the current simplex.
2. Locates the exit facet by solving a small linear system at each facet (the barycentric coordinates of the zero crossing).
3. Yields the crossing point on the exit facet.
4. Pivots into the adjacent simplex sharing that facet and repeats.

No derivatives or nonlinear solves are required; each step costs one function evaluation and one small linear solve (n×n).
