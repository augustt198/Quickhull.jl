using Test
import Random, DataStructures

using Quickhull
import QHull

include("utils.jl")

function vertexset(facets)
    foldl((s, f) -> union!(s, f.plane.point_indices), facets; init=Set{Int}())
end
vertexset(hull::QHull.Chull) = Set(hull.vertices)

for D in (2, 3)
    @testset verbose=true "$D-D hull" begin

        sample_ranges = Dict(
            gridbox => [5, 10, 20, 30],
            samplebox => [50, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
            sampleball => [50, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
            samplesphere => [50, 100, 1000, 10_000] # much harder
        )
        for (sample_func, Ns) in sample_ranges
            @testset "$sample_func" begin
                for N in Ns
                    Random.seed!(1738 + 2)
                    pts = sample_func(N, D)
                    pts_T = Matrix(pts')

                    t1 = @elapsed begin
                        my_hull = quickhull(pts, Quickhull.Options(kernel=Quickhull.HyperplaneKernelExact_A))
                    end
                    t2 = @elapsed begin
                        their_hull = QHull.chull(pts_T)
                    end
                    my_facets = Quickhull.finished_facets(my_hull.facets)

                    println("Elapsed ($sample_func $D / $N): $t1 -- $t2     [$(length(my_facets)) facets] [$(size(their_hull.simplices))]")
                    @test vertices(my_hull) == vertexset(their_hull)
                end
            end
        end
    end
end
