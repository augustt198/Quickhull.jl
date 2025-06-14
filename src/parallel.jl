function subhull_dc(pts, opts)
    hull, _ = subhull_dc_recur(pts, 1, length(pts), opts.subdivide.levels, opts)
    return hull
end

function subhull_dc_recur(pts, l, r, level, opts)
    chunk_len = cld(r - l + 1, opts.subdivide.chunks)
    chunk_starts = l:chunk_len:r
    tasks = map(chunk_starts) do base
        l′, r′ = base, min(base + chunk_len - 1, r)
        if level <= 1
            sub_pts = @view pts[l′:r′]
            if opts.subdivide isa ParallelSubdivide
                Threads.@spawn (_quickhull_main(sub_pts, opts), l′:r′)
            else
                (_quickhull_main(sub_pts, opts), l′:r′)
            end
        else
            if opts.subdivide isa ParallelSubdivide
                Threads.@spawn subhull_dc_recur(pts, l′, r′, level - 1, opts)
            else
                subhull_dc_recur(pts, l′, r′, level - 1, opts)
            end
        end
    end

    hulls = fetch.(tasks)

    vs = merge_vertices(hulls)
    hull_input = @view pts[vs]
    return _quickhull_main(hull_input, opts), vs
end

function merge_vertices(hulls)
    vs = Int[]
    vcount = sum(h -> length(vertices(h[1])), hulls)
    sizehint!(vs, vcount)

    for (h, idxs) in hulls
        for v in vertices(h)
            push!(vs, idxs[v])
        end
    end

    return vs
end
