# A point lifted into D+1 dimensions
struct LiftedPoint{P, T} <: AbstractVector{T}
    point::P
    lifted_coord::T
end

Base.eltype(::Type{LiftedPoint{P, T}}) where {P, T} = T
Base.length(::Type{LiftedPoint{P, T}}) where {P, T} = length(P) + 1
Base.size(::Type{LiftedPoint{P, T}}) where {P, T} = (length(P) + 1,)
Base.size(pt::LiftedPoint) = size(typeof(pt))

function Base.getindex(lp::LiftedPoint{P, T}, idx::Integer) where {P, T}
    @boundscheck !(1 <= idx <= length(P) + 1) && throw(BoundsError(lp, idx))

    if idx <= length(P)
        return @inbounds lp.point[idx]
    else
        return lp.lifted_coord
    end
end

# A vector of LiftedPoint, the lifted coordinates are stored
# separately so the points being lifted don't need to be reallocated
struct LiftedPoints{V, P, T} <: AbstractVector{LiftedPoint{P, T}}
    points::V
    lifted_coords::Vector{T}

    function LiftedPoints(points::AbstractVector{P}, lift=(pt) -> dot(pt, pt)) where {P}
        lps = new{typeof(points), P, eltype(P)}(points, lift.(points))

        # If the lifted coordinates are very similar, the interior
        # point computed by averaging the initial simplex vertices
        # vertices may not actually be inside the simplex. To avoid
        # this we can scale the lifted coordinates.
        min, max = extrema(lps.lifted_coords)
        c = 2^12
        if (min - max) < c * eps(min)
            lps.lifted_coords .= c .* (lps.lifted_coords .- min)
        end

        return lps
    end
end

Base.getindex(lps::LiftedPoints, idx::Integer) = LiftedPoint(lps.points[idx], lps.lifted_coords[idx])
Base.size(lps::LiftedPoints) = size(lps.points)

"""
    delaunay(points)

Compute the d-dimensional Delaunay triangulation of `points`.
`points` can be a vector of point-like objects (e.g. `Tuple` or
`StaticVector`) or a (D, N) sized matrix of numbers.

The triangulation is found by lifting into (d+1) dimensions
and taking the convex hull.
"""
delaunay(pts, opts=Options()) = quickhull(LiftedPoints(pts), opts)

delaunay(pts::Matrix, opts=Options()) = delaunay(matrix2points(pts), opts)
