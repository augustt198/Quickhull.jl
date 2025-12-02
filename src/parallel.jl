subhull_dc(pts, opts) = subhull_dc_recur(pts, 1, length(pts), opts.subdivide.levels, opts)

function subhull_dc_recur(pts, l, r, level, opts)
    D = pointsdim(pts)
    n = r - l + 1

    # a bit of extra complexity because we need to ensure that
    # we don't subdivide the input so that any chunk has less
    # than D + 1 points (in which case we can't even form a simplex
    # hull). So we ensure the calculated chunk lengths are not too
    # short and the last chunk isn't truncated too short either.
    chunk_len = max(D + 1, cld(n, opts.subdivide.chunks))
    chunk_starts = l:chunk_len:r
    last_start = last(chunk_starts)
    if r - last_start + 1 < D + 1
        chunk_starts = l:chunk_len:(last_start - 1)
        last_start = last(chunk_starts)
    end

    tasks = map(chunk_starts) do base
        l′ = base
        r′ = (base == last_start) ? r : base + chunk_len - 1

        if level <= 1
            sub_pts = @view pts[l′:r′]
            if opts.subdivide isa ParallelSubdivide
                Threads.@spawn _quickhull_main(sub_pts, opts)
            else
                _quickhull_main(sub_pts, opts)
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
    return _quickhull_main(hull_input, opts)
end

function merge_vertices(hulls)
    vs = Int[]
    vcount = sum(h -> length(vertices(h)), hulls)
    sizehint!(vs, vcount)

    for h in hulls
        idxs = parentindices(h.pts)[1]
        for v in vertices(h)
            push!(vs, idxs[v])
        end
    end

    return vs
end
