module Quickhull

using StaticArrays
using GeometryBasics
using LinearAlgebra
using Combinatorics
using Base.Iterators
using MacroTools
using MappedArrays

include("utils.jl")
include("predicate.jl")
include("hyperplane.jl")
include("smallvec.jl")
include("hull.jl")
include("delaunay.jl")
include("voronoi.jl")

export quickhull, delaunay,
    facets, vertices, points,
    voronoi_centers, voronoi_edges, voronoi_edge_points

@kwdef struct Options{K <: HyperplaneKernel, I <: Integer}
    # numerical kernel used for plane calculations
    kernel::Type{K} = HyperplaneKernelExact_A

    # integer type used as indices
    indextype::Type{I} = Int32

    # whether to joggle input points and by how much
    joggle::Bool = false
    joggle_amount::Float64 = 100.0

    # whether to record statistics
    statistics::Bool = false

    # mostly just have animations in mind
    iteration_callback::Any = nothing
end

function Base.show(io::IO, ::MIME"text/plain", opts::Options)
    println(io, "Quickhull Options:")
    for prop in propertynames(opts)
        println(io, "  - $(prop): $(getproperty(opts, prop))")
    end
end

# create a Hull from point indices defining a simplex
function makesimplexhull(pts::V, simp, ::Val{D}, I, K) where {V, D}
    T = eltype(eltype(V))
    K′ = typeof(make_kernel(K, SVector{D}(simp[1:end-1]), pts)) # better way to do this??

    hull = Hull(pts, zeros(SVector{D, T}), I, K′)

    data = IterData{D, T, I, K′}()

    simplex_facets = Facet{D, I, K′}[]
    for (i, pt_idx) in enumerate(simp)
        apex = pts[pt_idx]

        indices = simp[(1:i-1) ∪ (i+1:end)]
        plane = Hyperplane(SVector{D, I}(indices), pts, K)
        plane = hyperplane_awayfrom(plane, apex, pts)

        facet = Facet(plane)
        facet.adj = (1:i-1) ∪ (i+1:D+1)
        push!(simplex_facets, facet)
    end

    # initialize above point sets for the simplex facets
    mark_above(simplex_facets, 1:size(pts, 1), pts, data, false)
    foreach(f -> push_hull_facet!(hull.facets, f), simplex_facets)

    center = sum(pts[simp]) / length(simp)
    hull.interior_pt = center

    for f in hull.facets
        dist = hyperplane_dist(f.plane, hull.interior_pt, pts)
        if dist >= 0
            error("Interior point wasn't actually inside hull (todo fix)")
        end
    end
    
    return hull, data
end

# find a good initial simplex by random sampling
function goodsimplex(pts::V, ::Val{D}, I, K, nsamp=1000) where {V, D}
    T = eltype(eltype(V))
    N = length(pts)

    # oof can we end up with repeated indices?? 
    minidx   = MVector{D}(-1 for _ in 1:D)
    mincoord = MVector{D}(typemax(T) for _ in 1:D)
    maxidx   = MVector{D}(-1 for _ in 1:D)
    maxcoord = MVector{D}(typemin(T) for _ in 1:D)

    samples = min(N, nsamp)
    shouldsample = (samples == nsamp)
    @inbounds for i in 1:samples
        # either sample or just go consecutively
        idx = shouldsample ? rand(1:N) : i
        pt = pts[idx]

        for (d, c) in enumerate(pt)
            if c < mincoord[d]
                minidx[d], mincoord[d] = idx, c
            end
            if c > maxcoord[d]
                maxidx[d], maxcoord[d] = idx, c
            end
        end
    end

    minmaxidx = vcat(minidx, maxidx)
    maxvol = -one(T)
    maxsimp = zeros(SVector{D+1, I})
    for (i, comb) in enumerate(combinations(minmaxidx, D+1))
        (i > 1000) && break

        comb = SVector{D+1, I}(comb)
        vol = abs(simplexvolume(pointmatrix(pts, comb)))
        if vol > maxvol
            maxvol, maxsimp = vol, comb
        end
    end

    maxvol > 0 && return maxsimp

    # sampling only used repeated or coplanar points, so
    # exaustively search for a simplex with volume
    simplex = Hyperplane(SVector{1, I}(1), pts, K)
    for d = 2:(D+1)
        found = -1
        for (i, pt) in enumerate(pts)
            i ∈ simplex.point_indices && continue
            if hyperplane_dist(simplex, pt, pts) != 0
                found = i
                break
            end
        end
        if found == -1
            throw(ArgumentError("All the points are coplanar. Try projecting into a lower dimension."))
        end

        idxs = SVector{d, I}(simplex.point_indices..., found)
        simplex = Hyperplane(idxs, pts, K)
    end

    # loop above makes inference hard
    return SVector{D+1, I}(simplex.point_indices)::SVector{D+1, I}
end

function simplexvolume(pts)
    D, N = size(pts)
    mat, pt = droplast(pts, dims=(2,)), pts[:, end]

    if N == D+1 # mat is square
        return vol_exact(mat, pt) / factorial(N-1)
    else
        error()
        return sqrt(detn_exact(mat*mat')) / factorial(N-1)
    end
end

# Compute the arithmetic center of the hull's vertices
function center_average(hull)
    vs = vertices(hull)    
    return sum(hull.pts[i] for i in vs) / length(vs)
end

# Compute the area and volume of the hull.
# Area is found by summing the area of each simplex.
# Volume is found by summing the volume of pyramids,
# where the base is a facet and the apex is a fixed
# interior point. The volume of each pyramid is then
# (facet area)*(height to apex)/dimension.
function areavol(hull::Hull{D, T}) where {D, T}
    c = center_average(hull)
    
    areavol = sum(hull.facets) do f
        area = simplexvolume(pts[f.plane.point_indices])
        dist = hyperplane_dist(f.plane, c, hull.pts)

        SVector{2}(area, area * dist)
    end
    
    return areavol[1], areavol[2] / D
end

# assuming a and b don't have repeats
function count_intersection(a, b)
    c = 0
    for x in a
        if x in b
            c += 1
        end
    end
    return c
end

# replace a zero in v with x
function insertadj(v, x)
    for (i, y) in enumerate(v)
        (y == 0) && return setindex(v, x, i)
    end
    error("no more adj slots?")
    return v
end

# find the point indices that define the ridge shared by planes x and y.
# assumes x and y only differ by one point
@generated function ridgepoints(x::SVector{D, I}, y::SVector{D, I}) where {D, I}
    body = map(1:D) do i
        index = (1:i-1) ∪ (i+1:D)
        :( (x[$i] ∉ y) && return x[SVector{D-1}($(index...))] )
    end

    Expr(:block,
        body...,
        :(return zeros(SVector{D-1, I}))
    )
end

# So we don't reallocate these for every iteration
struct IterData{D, T, I, K}
    queue       ::Vector{Facet{D, I, K}}
    visible     ::Vector{Facet{D, I, K}}
    horizon     ::Vector{NTuple{2, Facet{D, I, K}}}
    newfacets   ::Vector{Facet{D, I, K}}
    cands       ::Vector{I}
    newplanes   ::Vector{Hyperplane{D, I, K}}
    maxidx      ::Vector{I}
    maxdist     ::Vector{T}

    function IterData{D, T, I, K}() where {D, T, I, K}
        new{D, T, I, K}((ft() for ft in fieldtypes(IterData{D, T, I, K}))...)
    end
end

function clear_iterdata(data)
    empty!(data.queue)
    empty!(data.visible)
    empty!(data.horizon)
    empty!(data.newfacets)
    empty!(data.cands)
    empty!(data.newplanes)
end

# partition the point indices in `candidates` to the above sets of facets in `fs`.
# `save` determines if unassigned candidate points should be saved in data.cands 
function mark_above(fs::AbstractVector{Facet{D, I, K}}, candidates, pts, data::IterData{D, T, I, K}, save) where {D, T, I, K}
    maxidx  = fill!(resize!(data.maxidx,  length(fs)), -1)
    maxdist = fill!(resize!(data.maxdist, length(fs)), -one(T))

    for idx in candidates
        pt = @inbounds pts[idx]
        marked = false
        for (fi, facet) in enumerate(fs)
            plane = facet.plane
            
            # don't check points that define this plane
            (idx in plane.point_indices) && continue
        
            dist = hyperplane_dist(plane, pt, pts)
            if dist > 0 # 'above'
                push!(facet.above, idx)
                
                if dist > (@inbounds maxdist[fi])
                    @inbounds maxidx[fi]  = idx
                    @inbounds maxdist[fi] = dist
                end
                marked = true
                break
            end
        end

        if save && !marked
            push!(data.cands, idx)
        end
    end

    for (fi, facet) in enumerate(fs)
        facet.furthest_above_point = @inbounds maxidx[fi]
    end
end

# Do an iteration of quickhull: insert the point furthest
# above `facet` into the hull.
function iter(hull::Hull{D, T, I, K, V}, facet, data) where {D, T, I, K, V}
    furthest_pt_idx = facet.furthest_above_point
    furthest_pt = hull.pts[furthest_pt_idx]

    # STEP 1: starting with `facet`, find facets are 'visible'
    # from `pt`, and create the horizon defined by the boundary
    # of these facets.
    queue, visible, horizon = data.queue, data.visible, data.horizon
    push!(queue, facet)
    push!(visible, facet)
    while !isempty(queue)
        f = pop!(queue)
        for h_i in f.adj
            h = hull.facets.arr[h_i]
            dist = hyperplane_dist(h.plane, furthest_pt, hull.pts)
            if dist < 0
                # the point was above f but below h, so
                # (f, h) defines a ridge on the horizon
                push!(horizon, (f, h))
            elseif dist >= 0
                if h ∉ visible
                    push!(visible, h)
                    push!(queue, h)
                end
            end
        end
    end

    # how many of the facets in `newfacets` are being
    # reused from existing facets
    reallocated = 0

    # STEP 2: build cone of facets from the horizon to `pt`
    newfacets, newplanes = data.newfacets, data.newplanes
    for (i, (f, h)) in enumerate(horizon)
        # the D-1 point ridge where f and h meet
        # ie: ridge = f.plane.point_indices ∩ h.plane.point_indices
        ridge = ridgepoints(f.plane.point_indices, h.plane.point_indices)
        newplane = Hyperplane(SVector{D}(ridge..., furthest_pt_idx), hull.pts, K)
        # make it face outwards
        newplane = hyperplane_awayfrom(newplane, hull.interior_pt, hull.pts)

        # We need to compute all the new planes before assigning
        # the planes to facets because we're reusing facets
        push!(newplanes, newplane)

        if i <= length(visible)
            newfacet = visible[i]
            reallocated += 1
        elseif has_unused_facet(hull.facets)
            newfacet = pop_unused_facet!(hull.facets)
            reallocated += 1
        else
            newfacet = Facet(newplane)
        end
        push!(newfacets, newfacet)
    end

    # for convenience split the new facets by
    # previously or newly allocated
    prev_allocated = @view newfacets[1:reallocated]
    new_allocated  = @view newfacets[reallocated+1:end]

    # now we can assign the planes after they've all been computed
    for (i, (f, h)) in enumerate(horizon)
        newfacets[i].adj = setindex(zero(SVector{D, I}), h.handle, 1)
        newfacets[i].plane = newplanes[i]
        
        # fix adjacency list of the facet not visible ("behind" the horizon)
        idx = findfirst(isequal(f.handle), h.adj)
        # since the new facets may not have a `handle` yet,
        # the index -i is taken to mean the ith facet in `newfacets`
        h.adj = setindex(h.adj, -i, idx)
    end

    # create the topology of the cone by linking adjacent facets
    for i = 1:length(newfacets)
        for j = 1:(i-1)
            f, h = newfacets[i], newfacets[j]
            c = count_intersection(f.plane.point_indices, h.plane.point_indices)
            # the planes interesect at D-1 points, so they share a ridge
            if c == D-1
                # see comment above about negative indices
                f.adj = insertadj(f.adj, -j)
                h.adj = insertadj(h.adj, -i)
            end
        end
    end

    # STEP 3: mark above point sets

    # first partition points above the visible facets to
    # the newly allocated facets. After this the remaining
    # points to be partitioned are in data.cand, so the
    # point sets of the previously allocated facets can
    # be clobbered
    itr = Iterators.flatten(Iterators.map(f -> f.above, visible))
    mark_above(new_allocated, itr, hull.pts, data, true)

    foreach(f -> empty!(f.above), prev_allocated)
    mark_above(prev_allocated, data.cands, hull.pts, data, false)
    
    # remove all the visible facets from the hull
    for (i, f) in enumerate(visible)
        remove_hull_facet!(hull.facets, f)
        if i > reallocated
            # if the facet was in the hull but not
            # reallocated, mark it as unused for use later
            push_unused_facet!(hull.facets, f)
        end
    end

    for f in newfacets
        push_hull_facet!(hull.facets, f)
    end

    # now that all facets have a handle, fix any negative
    # indices in the adjacency lists
    fixup_adj!(f) = f.adj = map(f.adj) do adj_idx
        adj_idx < 0 ? (newfacets[-adj_idx].handle % eltype(f.adj)) : adj_idx
    end

    foreach(fixup_adj!, newfacets)
    # we also need to fix facets that were below the horizon
    foreach(fixup_adj! ∘ last, horizon)

end


dimcheck(D, N) = (N < D+1) && throw(ArgumentError("Need at least $(D+1) points in $D dimensions (got $N points)"))

"""
    quickhull(points, options=Quickhull.Options())

Compute the convex hull of `points`. `points` can be a vector of
point-like objects (e.g. `Tuple` or `StaticVector`) or a (D, N) sized matrix
of numbers.

See documentation for `Quickhull.Options`.

`kernel` can optionally be specified to control how geometric predicates
are computed. For instance, `HyperplaneKernelInexact` is fast not numerically
robust, while `HyperplaneKernelExact_A` is slightly slower but robust.
"""
function quickhull(pts::AbstractMatrix{T}, opts::O=Options()) where {T, O <: Options}
    D, N = size(pts)
    dimcheck(D, N)

    v = dropdims(reinterpret(SVector{D, T}, pts), dims=1)

    return _quickhull(v, Val(D), opts)
end

function quickhull(pts::V, opts::O=Options()) where {V <: AbstractVector, O <: Options}
    D, N = length(eltype(pts)), length(pts)
    dimcheck(D, N)

    return _quickhull(pts, Val(D), opts)
end

function _quickhull(pts::V, ::Val{D}, opts) where {V, D}
    if opts.joggle
        pts = joggle(pts, opts.joggle_amount)
    end

    simplex = goodsimplex(pts, Val(D), opts.indextype, opts.kernel)

    hull, data = makesimplexhull(pts, simplex, Val(D), opts.indextype, opts.kernel)

    while true
        head = hull.facets.working_list_head
        if head < 1 # we're done!
            compact_unused_facets!(hull.facets)
            clear_iterdata(data)
            break
        end

        facet = hull.facets.arr[head]
        iter(hull, facet, data)
        clear_iterdata(data)
    end

    return hull
end

end # module Quickhull
