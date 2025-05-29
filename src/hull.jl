abstract type AbstractHull{D, T, I} end

"""
    points(hull)

The points the hull was constructed from. This includes points
inside the hull - see `vertices(hull)` and `vertexpoints(hull)`.
"""
points(hull::AbstractHull) = error("unimplemented")

"""
    vertices(hull)

The indices of points that are vertices of the hull.
"""
vertices(hull::AbstractHull) = error("unimplemented")

"""
    vertexpoints(hull)

The points that are vertices of the hull.
"""
vertexpoints(hull::AbstractHull) = @view points(hull)[vertices(hull)]

"""
    facets(hull)

The facets of the hull. A facet is defined by D vertices.
"""
facets(hull::AbstractHull) = error("unimplemented")

"""
    GeometryBasics.Mesh(hull::Quickhull.AbstractHull)

Create a `Mesh` from the points and facets of `hull`.
"""
function GeometryBasics.Mesh(hull::AbstractHull)
    GeometryBasics.Mesh(points(hull), facets(hull))
end


const NULL_INDEX        = -2
const FLAG_INDEX_UNUSED = -3
const FLAG_INDEX_DONE   = -4

@embed_SmallVec mutable struct Facet{D, I, K}
    plane::Hyperplane{D, I, K}
    adj::SVector{D, I}
    
    handle::Int # i.e. this facet's index in a FacetList

    # pointers for use as a linked list in FacetList
    prev_handle::Int
    next_handle::Int

    furthest_above_point::I

    # Vector of points above this facet
    above::@SmallVec{I, 4}

    function Facet(plane::Hyperplane{D, I, K}) where {D, I, K}
        new{D, I, K}(
            plane, zero(SVector{D, I}),
            NULL_INDEX, NULL_INDEX,
            0, -1)
    end
end

# two linked lists backed by an array:
#  - working list, containing facets that need to be
#    worked on, i.e. have > 0 points above them.
#    This list is doubly linked.
#  - unused list, containing facets that are no
#    longer in use and can be reused. Singly linked.
mutable struct FacetList{D, I, K}
    arr::Vector{Facet{D, I, K}}

    working_list_head::Int
    unused_list_head::Int

    valid_handles::Bool
    compacted::Bool

    FacetList{D, I, K}() where {D, I, K} = new{D, I, K}(Facet{D, K}[], NULL_INDEX, NULL_INDEX, true, false)
end

function list_length(fl::FacetList, handle)
    len = 0
    while handle != NULL_INDEX
        len += 1
        handle = fl.arr[handle].next_handle
    end
    return len
end

working_list_length(fl::FacetList) = list_length(fl, fl.working_list_head)
unused_list_length(fl::FacetList)  = list_length(fl, fl.unused_list_head)


Base.length(fl::FacetList) = length(fl.arr)
Base.eltype(fl::FacetList) = eltype(fl.arr)
Base.iterate(fl::FacetList) = iterate(fl.arr)
Base.iterate(fl::FacetList, state) = iterate(fl.arr, state)

# push a facet that is used in the hull currently
function push_hull_facet!(fl::FacetList, facet)
    if facet.handle == NULL_INDEX
        push!(fl.arr, facet)
        facet.handle = length(fl.arr)
    end
    facet.prev_handle = NULL_INDEX

    if isempty(facet.above)
        # The facet has no points in its above set, so mark
        # it as done. This doesn't necessarily mean this facet
        # will exist in the finished hull.
        facet.next_handle = FLAG_INDEX_DONE
    else
        facet.next_handle = fl.working_list_head
        fl.working_list_head = facet.handle

        if facet.next_handle != NULL_INDEX
            fl.arr[facet.next_handle].prev_handle = facet.handle
        end
    end
end

function remove_hull_facet!(fl::FacetList, facet)
    if facet.next_handle == FLAG_INDEX_DONE
        # Done facets are already removed from the working
        # list so there isn't anything to do
        return
    end

    if facet.next_handle != NULL_INDEX
        fl.arr[facet.next_handle].prev_handle = facet.prev_handle
    end

    if fl.working_list_head == facet.handle
        @assert facet.prev_handle < 1
        fl.working_list_head = facet.next_handle
    else
        fl.arr[facet.prev_handle].next_handle = facet.next_handle
    end
end

function push_unused_facet!(fl::FacetList, facet)
    facet.next_handle = fl.unused_list_head
    facet.prev_handle = FLAG_INDEX_UNUSED
    fl.unused_list_head = facet.handle
end

has_unused_facet(fl::FacetList) = fl.unused_list_head != NULL_INDEX

function pop_unused_facet!(fl::FacetList)
    if !has_unused_facet(fl)
        throw(ArgumentError("list empty"))
    end
    facet = fl.arr[fl.unused_list_head]
    fl.unused_list_head = facet.next_handle

    return facet
end

# Remove any facets marked as unused
# Wait this will mess up facet handles...
function compact_unused_facets!(fl::FacetList)
    fl.unused_list_head = NULL_INDEX
    i′ = 1
    for f in fl.arr
        if f.prev_handle != FLAG_INDEX_UNUSED
            fl.arr[i′] = f
            i′ += 1
        end
    end
    resize!(fl.arr, i′ - 1)

    fl.compacted = true
    fl.valid_handles = false
    return fl
end

function make_handles_valid!(fl::FacetList)
    fl.valid_handles && return fl

    handlemap = Dict{Int, Int}()
    for (i, f) in enumerate(fl.arr)
        handlemap[f.handle] = i
        f.handle = i
    end

    for f in fl.arr
        f.next_handle >= 1 && (f.next_handle = handlemap[f.next_handle])
        f.prev_handle >= 1 && (f.prev_handle = handlemap[f.prev_handle])
        
        f.adj = map(h -> handlemap[h], f.adj)
    end

    fl.valid_handles = true
    return fl
end

finished_facets(fl::FacetList) = filter(f -> isempty(f.above) && f.next_handle == FLAG_INDEX_DONE, fl.arr)

Base.@kwdef struct IterStat
    nvisible::Int
    nnew::Int
    ncands::Int
    duration_ns::Int
end

mutable struct Hull{D, T <: Number, I <: Integer, K, V <: AbstractVector} <: AbstractHull{D, T, I}
    pts::V
    facets::FacetList{D, I, K}
    vertices::Union{Vector{I}, Nothing}
    approx_interior_pt::SVector{D, T}

    statistics::Union{Nothing, Vector{IterStat}}

    Hull(pts::V, approx_interior_pt::SVector{D, T}, ::Type{I}, ::Type{K}) where {D, T, I, K, V} = new{D, T, I, K, V}(
        pts,
        FacetList{D, I, K}(),
        nothing,
        approx_interior_pt,
        nothing
    )
end

indextype(::Hull{D, T, I}) where {D, T, I} = I

function Base.show(io::IO, ::MIME"text/plain", hull::AbstractHull{D, T, I}) where {D, T, I}
    println(io, "Hull of $(length(points(hull))) points in $D dimensions")
    print(io, "  - $(length(vertices(hull))) Hull vertices: ")
    println(io, vertices(hull))
    print(io, "  - $(length(facets(hull))) Hull facets: ")
    println(io, facets(hull))
end

function vertices(hull::Hull)
    I = indextype(hull)
    if isnothing(hull.vertices)
        vset = Set{I}()
        for f in facets(hull)
            push!(vset, f...)
        end
        hull.vertices = collect(vset)
    end

    return hull.vertices::Vector{I}
end

points(hull::Hull) = mappedarray(pt -> Point(pt), hull.pts)

function facets(hull::Hull)
    if !hull.facets.compacted
        compact_unused_facets!(hull.facets)
    end

    return mappedarray(f -> NgonFace(f.plane.point_indices), hull.facets.arr)
end

function Base.insert!(hull::Hull{D, T, I, K, V}, pt) where {D, T, I, K, V}
    facet = find_facet_below(hull, pt)
    isnothing(facet) && return false

    push!(hull.pts, pt)
    pt_idx = lastindex(hull.pts)
    push!(facet.above, pt_idx)
    facet.furthest_above_point = pt_idx

    data = IterData{D, T, I, K}()
    iter(hull, facet, data)
    hull.facets.compacted = false
    return true
end

Base.in(pt, hull::Hull) = isnothing(find_facet_below(hull, pt))

function find_facet_below(hull, pt)
    make_handles_valid!(hull.facets)
    facet = first(hull.facets.arr)

    curr_dist = hyperplane_dist(facet.plane, pt, hull.pts)

    iters, maxiters = 0, 1000
    while curr_dist <= 0 && (iters += 1) <= maxiters
        maximum_dist, idx = findmax(facet.adj) do h_handle
            h = hull.facets.arr[h_handle]
            hyperplane_dist(h.plane, pt, hull.pts)
        end

        if maximum_dist > curr_dist
            curr_dist = maximum_dist
            facet = hull.facets.arr[facet.adj[idx]]
        else
            break
        end
    end

    if curr_dist <= 0
        handle = findfirst(hull.facets.arr) do f
            hyperplane_dist(f.plane, pt, hull.pts) > 0
        end

        return isnothing(handle) ? nothing : hull.facets.arr[handle]
    else
        return facet
    end
end
