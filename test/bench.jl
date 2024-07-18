using Quickhull
using BenchmarkTools
using Profile
using Random
using QHull

using GLMakie

include("utils.jl")

SPHERE_POLY = GLMakie.GeometryBasics.Polygon(
    decompose(Point2f, Circle(Point2f(0), .5)),
    [decompose(Point2f, Circle(Point2f(0), 0.3))]
)

MARKER_DICT = Dict(
    samplebox => :rect,
    sampleball => :circle,
    samplesphere => SPHERE_POLY
)

function bench()
    sample_ranges = Dict(
        samplebox => [100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
        sampleball => [100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
        samplesphere => [100, 1000, 10_000, 100_000, 1_000_000] # much harder
    )

    fig = Figure()
    ax1 = Axis(fig[1, 1], title="2 Dimensions", xscale=log10)

    ax2 = Axis(fig[2, 1], title="3 Dimensions", xscale=log10)

    for (D, ax) in zip((2, 3), (ax1, ax2))

        for (sample_func, Ns) in sample_ranges
            println("Running $(D)D $sample_func with $Ns points")
            data_mine_inexact = map(Ns) do N
                opts = Quickhull.Options(kernel=Quickhull.HyperplaneKernelInexact)
                res = @benchmark quickhull(pts, $opts) setup=(Random.seed!(1738); pts=$sample_func($N, $D))
                median(res).time / 1e9
            end
            data_mine_exact = map(Ns) do N
                opts = Quickhull.Options(kernel=Quickhull.HyperplaneKernelExactSIMD)
                res = @benchmark quickhull(pts, $opts) setup=(Random.seed!(1738); pts=$sample_func($N, $D))
                median(res).time / 1e9
            end
            data_qhull = map(Ns) do N
                res = @benchmark QHull.chull(pts) setup=(Random.seed!(1738); pts=Matrix($sample_func($N, $D)') )
                median(res).time / 1e9
            end

            marker = MARKER_DICT[sample_func]
            line = scatterlines!(ax, Ns, data_qhull ./ data_mine_inexact, label="$(string(sample_func)[7:end]) (inexact)", marker=marker)
            scatterlines!(ax, Ns, data_qhull ./ data_mine_exact, label="$(string(sample_func)[7:end]) (exact)", color=line.color, linestyle=:dot, marker=marker)
        end
        hlines!(ax, [1], color=:gray)
    end

    Label(fig[1:2,0], "Speedup relative to QHull", rotation=pi/2)
    Label(fig[3,:], "Number of points")
    fig[1:2,2] = Legend(fig, ax1, "Test case")

    fig
end

bench()
