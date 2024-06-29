# Quickhull.jl

The quickhull algorithm in pure Julia for finding
convex hulls, Delaunay triangulations, and Voronoi diagrams in N dimensions.

```julia
julia> using Quickhull

julia> hull = quickhull(randn(3, 500))
Hull of 500 points in 3 dimensions
  - 31 Hull vertices: Int32[297, 438  …  147, 376]
  - 58 Hull facets: TriangleFace{Int32}[TriangleFace(139, 249, 243)  …  TriangleFace(104, 147, 243)]

julia> using GLMakie, GeometryBasics
julia> wireframe(GeometryBasics.Mesh(hull))
julia> scatter!(hull.pts, color=:black)
```

<p align="center"><img src="img/hull.png" width="50%"></p>

```julia
julia> tri = delaunay(rand(2, 100));
julia> f = Figure()
julia> wireframe!(Axis(f[1,1]), GeometryBasics.Mesh(tri))
julia> linesegments!(Axis(f[1,2]), voronoi_edge_points(tri), color=:red)
```

<p align="center"><img src="img/tri.png" width="80%"></p>


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
