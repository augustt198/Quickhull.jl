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
    marker::Int

    furthest_above_point::I

    # Vector of points above this facet
    above::@SmallVec{I, 4}

    function Facet(plane::Hyperplane{D, I, K}) where {D, I, K}
        new{D, I, K}(
            plane, zero(SVector{D, I}),
            NULL_INDEX, NULL_INDEX, NULL_INDEX,
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

    FacetList{D, I, K}() where {D, I, K} = new{D, I, K}(Facet{D, K}[], NULL_INDEX, NULL_INDEX)
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
end

finished_facets(fl::FacetList) = filter(f -> isempty(f.above) && f.next_handle == FLAG_INDEX_DONE, fl.arr)

@kwdef struct IterStat
    nseen::Int
    nnew::Int
    ncands::Int
    nmaxabovefacet::Int
    marks::Vector{Int}
    time::Int
end

mutable struct Hull{D, T <: Number, I <: Integer, K, V <: AbstractVector}
    pts::V
    facets::FacetList{D, I, K}
    vertices::Union{Vector{I}, Nothing}
    interior_pt::SVector{D, T}

    statistics::Union{Nothing, Vector{IterStat}}

    Hull(pts::V, interior_pt::SVector{D, T}, ::Type{I}, ::Type{K}) where {D, T, I, K, V} = new{D, T, I, K, V}(
        pts,
        FacetList{D, I, K}(),
        nothing,
        interior_pt,
        nothing
    )
end

indextype(::Hull{D, T, I}) where {D, T, I} = I

function Base.show(io::IO, ::MIME"text/plain", hull::Hull{D, T, I, K, V}) where {D, T, I, K, V}
    println(io, "Hull of $(length(hull.pts)) points in $D dimensions")
    println(io, "  - Point type: $(eltype(V))")
    println(io, "  - Kernel type: $K")
    print(io, "  - $(length(vertices(hull))) Hull vertices: ")
    println(io, vertices(hull))
    print(io, "  - $(length(facets(hull))) Hull facets: ")
    println(io, collect(facets(hull)))
end

facets(hull::Hull) = Iterators.map(f -> f.plane.point_indices, hull.facets)

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

# Filter out facets not on the bottom of the convex hull
function delaunay_facets(hull::Hull{D}) where {D}
    maxlift = maximum(last, hull.pts)

    return Iterators.filter(hull.facets) do facet
        above_pt = sum(SVector{D}(hull.pts[i]) for i in facet.plane.point_indices) / D
        above_pt = setindex(above_pt, above_pt[end] + 2maxlift, D)
        hyperplane_dist(facet.plane, above_pt, hull.pts) < 0
    end
end
