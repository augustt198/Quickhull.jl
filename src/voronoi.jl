function voronoi_centers(hull)
    pts = points(hull)
    map(facets(hull)) do f
        c, r = circumsphere(pts[f])
        return c
    end
end

function voronoi_edges(hull)
    fs = delaunay_facets(hull)

    edges = Tuple{Int, Int}[]

    handle_map = Dict{Int, Int}()
    for (i, f) in enumerate(fs)
        handle_map[f.handle] = i
    end
    
    for (i, f) in enumerate(fs)
        for adj in f.adj
            # don't double count
            if adj < f.handle
                j = get(handle_map, adj, nothing)
                if !isnothing(j)
                    push!(edges, (i, j))
                end
            end
        end
    end

    return edges
end

function voronoi_edge_points(hull)
    pts = voronoi_centers(hull)
    edges = voronoi_edges(hull)

    [(pts[i], pts[j]) for (i, j) in edges]
end

function circumsphere(pts)
    n = length(pts)
    D = length(first(pts))
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
