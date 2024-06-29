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

    # build matrix for Cayley-Menger determinant
    M = zero(SMatrix{D+1, D+1})
    for i = 1:n
        for j = 1:(i-1)
            diff = pts[i] .- pts[j]
            d² = dot(diff, diff)
            M = setindex(M, d², i, j)
            M = setindex(M, d², j, i)
        end
    end

    x = M \ ones(SVector{D+1})
    inv_2R² = sum(x)
    # normalize to find barycentric coords
    bary = x ./ inv_2R²

    # convert to cartesian
    center = sum(i -> bary[i] * pts[i], 1:(D+1))
    radius = sqrt(inv(inv_2R²) / 2)
    return center, radius
end
