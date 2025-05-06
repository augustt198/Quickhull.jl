using Test
using Quickhull

import GeometryBasics
import Random

import Conda
Conda.add("scipy")

import QHull

include("utils.jl")

for D in (2, 3)
    @testset "$D-D hull" begin

        sample_ranges = Dict(
            gridbox => [5, 10, 20, 30],
            samplebox => [50, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
            sampleball => [50, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
            samplesphere => [50, 100, 1000, 10_000] # much harder
        )
        for (sample_func, Ns) in sample_ranges
            @testset "$sample_func" begin
                for N in Ns
                    Random.seed!(1234)
                    pts = sample_func(N, D)
                    pts_T = Matrix(pts')

                    t1 = @elapsed begin
                        my_hull = quickhull(pts, Quickhull.Options(kernel=Quickhull.HyperplaneKernelExactSIMD))
                    end
                    t2 = @elapsed begin
                        their_hull = QHull.chull(pts_T)
                    end
                    my_facets = Quickhull.finished_facets(my_hull.facets)

                    #println("Elapsed ($sample_func $D / $N): $t1 -- $t2     [$(length(my_facets)) facets] [$(size(their_hull.simplices))]")
                    @test Set(vertices(my_hull)) == Set(their_hull.vertices)
                end
            end
        end
    end
end

@testset "kernels" begin
    Random.seed!(1234)
    pts = sampleball(100_000, 3)
    their_hull = QHull.chull(Matrix(pts'))

    Ks = (Quickhull.HyperplaneKernelInexact,
        Quickhull.HyperplaneKernelExact_A,
        Quickhull.HyperplaneKernelExactSIMD)

    for K ∈ Ks
        my_hull = quickhull(pts, Quickhull.Options(kernel=K))
        @test Set(vertices(my_hull)) == Set(their_hull.vertices)
    end
end

@testset "inference" begin
    @inferred quickhull([(rand(), rand()) for i = 1:100])
    @inferred quickhull([rand(SVector{3}) for i = 1:100])
    @inferred quickhull([rand(GeometryBasics.Point3f) for i = 1:100])

    hull = quickhull([rand(GeometryBasics.Point3f) for i = 1:100])
    @inferred points(hull)
    @inferred vertices(hull)
    @inferred facets(hull)

    @inferred delaunay([(rand(), rand()) for i = 1:100])
    @inferred delaunay([rand(SVector{3}) for i = 1:100])
    @inferred delaunay([rand(GeometryBasics.Point3f) for i = 1:100])

    tri = delaunay([rand(GeometryBasics.Point2f) for i = 1:100])
    @inferred points(tri)
    @inferred vertices(tri)
    @inferred facets(tri)
end
