# Quickhull.jl

The quickhull algorithm in pure Julia for finding
convex hulls, Delaunay triangulations, and Voronoi diagrams in N dimensions.

```julia
using Quickhull

hull = quickhull(randn(100, 3))
```

### QHull Comparison

Quickhull.jl is competetive with QHull's performance,
although although it has fewer features.

Benchmarks go here
