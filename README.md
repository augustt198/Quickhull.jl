# Quickhull.jl

The quickhull algorithm in pure Julia for finding
convex hulls, Delaunay triangulations, and Voronoi diagrams in N dimensions.

```julia
julia> using Quickhull

julia> hull = quickhull(randn(3, 500))
Hull of 500 points in 3 dimensions
  - Point type: StaticArraysCore.SVector{3, Float64}
  - Kernel type: Quickhull.HyperplaneKernelExact_A{3, Float64, 3, 9}
  - 26 Hull vertices: Int32[133, 308  …  307, 414]
  - 48 Hull facets: StaticArraysCore.SVector{3, Int32}[[63, 493, 308], [493, 80, 86]  …  [419, 174, 365], [419, 274, 365]]

julia> using GLMakie
julia> wireframe(GLMakie.GeometryBasics.Mesh(hull))
julia> scatter!(hull.pts, color=:black)
```

<p align="center"><img src="img/hull.png" width="50%"></p>

## QHull Comparison

Quickhull.jl is competetive with QHull's performance even
when exact arithmetic is used, although although it has fewer features.

<p align="center"><img src="img/benchmark.png" width="80%"></p>

## Robustness

`quickhull` can be run with various hyperplane kernels. A hyperplane
kernel is a method of calculating hyperplane-point distances. By default,
an exact kernel is used (i.e. the sign of the distance is always correct)
to ensure robustness. Robustness can be traded for speed by choosing an inexact
kernel, for instance:

```julia
quickhull(pts, Quickhull.Options(kernel = Quickhull.HyperplaneKernelInexact))
```

It should be noted that if an inexact kernel is used – particularly
on inputs with coplanar or nearly coplanar points – the topology of the
hull can be become corrupted and an error will probably occur.
