# A point lifted into D+1 dimensions
struct LiftedPoint{D, P, T} <: StaticVector{D, T}
    point::P
    lifted_coord::T

    LiftedPoint(p, lc) = new{length(p)+1, typeof(p), typeof(lc)}(p, lc)

    function LiftedPoint{D, P, T}(tup::Tuple) where {D, P, T}
        head, last = tup[1:end-1], tup[end]
        return LiftedPoint(P(head), last)
    end
end

Base.Tuple(lp::LiftedPoint) = (lp.point..., lp.lifted_coord)

function StaticArrays.similar_type(::Type{LiftedPoint{D, P, T}}, et::Type{NewElType}, ::Size{S}) where {D, P, T, NewElType, S}
    if length(S) == 1
        # I'm thinking maybe we should just not deal with tuples
        # acting as points and reinterpret them instead...
        if P <: Tuple
            Pnew = NTuple{S[1]-1, NewElType}
            return LiftedPoint{S[1], Pnew, NewElType}
        else
            return LiftedPoint{S[1], similar_type(P, et, Size(S[1]-1)), NewElType}
        end
    else
        error("unimplemented")
    end
end

function Base.getindex(lp::LiftedPoint{D, P, T}, idx::Int) where {D, P, T}
    @boundscheck !(1 <= idx <= D) && throw(BoundsError(lp, idx))

    if idx < D
        return @inbounds lp.point[idx]
    else
        return lp.lifted_coord
    end
end

# A vector of LiftedPoint, the lifted coordinates are stored
# separately so the points being lifted don't need to be reallocated
struct LiftedPoints{D, V, P, T} <: AbstractVector{LiftedPoint{D, P, T}}
    points::V
    lifted_coords::Vector{T}

    function LiftedPoints(points::AbstractVector{P}, lift=(pt) -> dot(pt, pt)) where {P}
        D = pointsdim(points) + 1
        lps = new{D, typeof(points), P, eltype(P)}(points, lift.(points))

        # If the lifted coordinates are very similar, the interior
        # point computed by averaging the initial simplex vertices
        # vertices may not actually be inside the simplex. To avoid
        # this we can scale the lifted coordinates.
        min, max = extrema(lps.lifted_coords)
        c = 2^12
        if (max - min) < c * eps(min)
            lps.lifted_coords .= c .* (lps.lifted_coords .- min)
        end

        return lps
    end

    function LiftedPoints(points::AbstractVector{P}, lifted_coords::Vector{T}) where {P, T}
        D = pointsdim(points) + 1
        new{D, typeof(points), P, T}(points, lifted_coords)
    end
end

function Base.similar(lps::LiftedPoints, eltype::Type{LiftedPoint{D, Q, S}}, dims::Dims{1}) where {D, Q, S}
    points = similar(lps.points, Q, dims)
    lifted_coords = similar(lps.lifted_coords, S, dims)
    LiftedPoints(points, lifted_coords)
end

Base.getindex(lps::LiftedPoints, idx::Integer) = LiftedPoint(lps.points[idx], lps.lifted_coords[idx])

function Base.setindex!(lps::LiftedPoints, lp, i)
    lps.points[i] = lp.point
    lps.lifted_coords[i] = lp.lifted_coord
end

Base.size(lps::LiftedPoints) = size(lps.points)

struct DelaunayHull{D, T, I, H <: Hull{D, T, I}} <: AbstractHull{D, T, I}
    hull::H
end

# Filter out facets not on the bottom of the convex hull
function delaunay_facets(dhull::DelaunayHull)
    hull = dhull.hull
    maxlift = maximum(last, hull.pts)

    return filter(hull.facets.arr) do facet
        D = pointsdim(hull.pts)
        above_pt = sum(hull.pts[i] for i in facet.plane.point_indices) / D
        above_pt = setindex(above_pt, above_pt[end] + 2maxlift, D)
        hyperplane_dist(facet.plane, above_pt, hull.pts) < 0
    end
end

points(dhull::DelaunayHull) = mappedarray(lp -> Point(lp.point), dhull.hull.pts)
vertices(dhull::DelaunayHull) = vertices(dhull.hull)
facets(dhull::DelaunayHull) = mappedarray(f -> NgonFace(f.plane.point_indices), delaunay_facets(dhull))


"""
    delaunay(points, options=Quickhull.Options())

Compute the d-dimensional Delaunay triangulation of `points`.
`points` can be a vector of point-like objects (e.g. `Tuple` or
`StaticVector`) or a (D, N)-sized matrix of numbers.

The triangulation is found by lifting into d+1 dimensions
and taking the convex hull.
"""
function delaunay(pts, opts=Options())
    if !(opts.subdivide isa NoSubdivide)
        throw(ArgumentError("Using `delaunay` with subdivision is strictly slower than without subdivision."))
    end

    return DelaunayHull(quickhull(LiftedPoints(pts), opts))
end

delaunay(pts::Matrix, opts=Options()) = delaunay(matrix2points(pts), opts)
