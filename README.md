# Quickhull.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://augustt198.github.io/Quickhull.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://augustt198.github.io/Quickhull.jl/dev)

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
julia> scatter!(points(hull), color=:black)
```

<p align="center"><img src="img/hull.png" width="50%"></p>

```julia
julia> tri = delaunay(rand(2, 100));
julia> f = Figure()
julia> wireframe!(Axis(f[1,1]), GeometryBasics.Mesh(tri))
julia> linesegments!(Axis(f[1,2]), voronoi_edge_points(tri), color=:red)
```

<p align="center"><img src="img/tri.png" width="80%"></p>


## Qhull Comparison

Quickhull.jl is competitive with [Qhull](http://www.qhull.org/)'s performance even
when exact arithmetic is used, although it has fewer features.

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
hull can become corrupted, and an error will probably occur.

## Related Packages
- [QHull.jl](https://github.com/JuliaPolyhedra/QHull.jl), [DirectQHull.jl](https://github.com/JuhaHeiskala/DirectQhull.jl/) – wrappers for [Qhull](http://www.qhull.org/)
- [DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl) – 2D Delaunay triangulation and Voronoi tesselation with many features
- [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl) – wrapper for 3D Delaunay tetrahedralization
library [TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1)
