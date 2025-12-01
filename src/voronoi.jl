"""
    voronoi_centers(delaunay_hull)

The circumcenters of a Delaunay triangulation's faces.
"""
function voronoi_centers(hull)
    pts = points(hull)
    map(facets(hull)) do f
        c, r = circumsphere(pts[f])
        return c
    end
end

struct VoronoiRay
    pt::Int
    ridge::Vector{Int}
    opposite::Int
end

function _voronoi_edges_and_rays(hull, do_rays=false)
    fs = delaunay_facets(hull)

    edges = Tuple{Int, Int}[]

    handle_map = Dict{Int, Int}()
    for (i, f) in enumerate(fs)
        handle_map[f.handle] = i
    end

    rays = VoronoiRay[]

    # make this better once `adj` has proper ordering
    vert_count_map = Dict{Int, Int}()
    function compute_ray(i, f, exclude_adj)
        empty!(vert_count_map)
        for adj in f.adj
            adj == exclude_adj && continue
            h = hull.hull.facets.arr[adj]
            for v in ridgepoints(f.plane.point_indices, h.plane.point_indices)
                cnt = get(vert_count_map, v, 0)
                vert_count_map[v] = cnt + 1
            end
        end

        D = length(f.adj)
        ridge = Int[]
        for (v, cnt) in vert_count_map
            if cnt == D-2
                push!(ridge, v)
            end
        end
        
        opposite = -1
        for idx in f.plane.point_indices
            if vert_count_map[idx] == D-1
                opposite = idx
            end
        end

        push!(rays, VoronoiRay(i, ridge, opposite))
    end
    
    for (i, f) in enumerate(fs)
        for adj in f.adj
            j = get(handle_map, adj, nothing)
            if isnothing(j)
                do_rays && compute_ray(i, f, adj)
            elseif adj < f.handle # don't double count
                push!(edges, (i, j))
            end
        end
    end

    return (edges = edges, ray_ridges = rays)
end

"""
    voronoi_edges(delaunay_hull)

Return a vector of edges that form the Voronoi diagram
(excluding rays). Each edge is a tuple `(v1, v2)` of point
indices into `voronoi_centers(delaunay_hull)`. 
"""
voronoi_edges(hull) = _voronoi_edges_and_rays(hull, false).edges

"""
    voronoi_edges_points(delaunay_hull)

Return a vector of edges that form the Voronoi diagram
(excluding rays). Each edge is a tuple `(p1, p2)` of line segment
endpoints.
"""
function voronoi_edge_points(hull)
    pts = voronoi_centers(hull)
    edges = voronoi_edges(hull)

    [(pts[i], pts[j]) for (i, j) in edges]
end

"""
    voronoi_edge_points_homogeneous(delaunay_hull)

Return a vector of edges that form the Voronoi diagram
(including rays). Each edge is a tuple `(p1, p2)` of line
segment endpoints in homogeneous coordinates. This means
an extra coordinate has been added to the end of each point:
if the coordinate is 1 the point is a 'normal' point, and if
the coordinate is 0 the point is 'at infinity'.
"""
function voronoi_edge_points_homogeneous(hull)
    pts = voronoi_centers(hull)
    D = pointsdim(pts)
    edges, rays = _voronoi_edges_and_rays(hull, true)

    edge_pts = map(edges) do (i, j)
        (Point(pts[i]..., 1), Point(pts[j]..., 1))
    end

    for ray in rays
        pt = Point(pts[ray.pt]..., 1)
        ridge = SVector{D}(ray.ridge)
        plane = Hyperplane(ridge, points(hull), HyperplaneKernelInexact)
        opp_pt = points(hull)[ray.opposite]
        plane = hyperplane_awayfrom(plane, opp_pt, points(hull))
        ray_pt = Point(plane.kernel.normal..., 0)
        
        push!(edge_pts, (pt, ray_pt))
    end

    return edge_pts
end

"""
    voronoi_edge_points_projected(delaunay_hull, raylength=1)

Return a vector of edges that form the Voronoi diagram
(including rays). Each edge is a tuple `(p1, p2)` of line
segment endpoints. Rays are projected out to be length `raylength`.
"""
function voronoi_edge_points_projected(hull, raylength=1)
    edge_pts = voronoi_edge_points_homogeneous(hull)

    projected = map(edge_pts) do (p1, p2)
        p1′, p2′ = droplast(p1), droplast(p2)
        if p2[end] > 0
            return (p1′, p2′)
        else
            return (p1′, p1′ + raylength*p2′)
        end
    end

    return projected
end

function vertex_adjacency_map(hull)
    adj_map = Dict{Int, Vector{Int}}()
    for (i, f) in enumerate(facets(hull))
        for vert in f
            adj = get!(() -> Int[], adj_map, vert)
            push!(adj, i)
        end
    end

    adj_map
end

function voronoi_cell_hulls(hull)
    adj_map = vertex_adjacency_map(hull)
    pts = voronoi_centers(hull)
    D = pointsdim(pts)

    filter!(adj_map) do (vert, adj)
        return length(adj) >= D + 1
    end

    dictmap(f, d) = [f(k, v) for (k, v) in d]

    dictmap(adj_map) do vert, adj
        cell_hull = quickhull(pts[adj])
        vert, cell_hull
    end
end

function circumsphere(pts)
    n = length(pts)
    D = pointsdim(pts)
    @assert n == D + 1

    A = zero(SMatrix{D, D})
    b = zero(SVector{D})
    for i = 1:D
        for j = 1:D
            a_ij = pts[n][j] - pts[i][j]
            A = setindex(A, 2a_ij, i, j)
        end
        rhs = dot(pts[n], pts[n]) - dot(pts[i], pts[i])
        b = setindex(b, rhs, i)
    end

    center = A \ b
    radius = norm(center - pts[n])
    return Point(center...), radius
end
