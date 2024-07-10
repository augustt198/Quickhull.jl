function merge_vertices(hulls)
    vs_set = Set{Int}()
    for (h, idxs) in hulls
        for v in vertices(h)
            push!(vs_set, idxs[v])
        end
    end
    return vs_set
end

function subhull_dc(pts, opts)
    _subhull_dc(pts, 1, length(pts), opts.subdivide.levels, opts)[1]
end

function _subhull_dc(pts, l, r, level, opts)
    chunk_len = cld(r - l, opts.subdivide.chunks)
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
                Threads.@spawn _subhull_dc(pts, l′, r′, level - 1, subdivopt)
            else
                _subhull_dc(pts, l′, r′, level - 1, subdivopt)
            end
        end
    end

    hulls = fetch.(tasks)

    vs = collect(merge_vertices(hulls))
    hull_input = @view pts[vs]
    return _quickhull_main(hull_input, opts), vs
end
