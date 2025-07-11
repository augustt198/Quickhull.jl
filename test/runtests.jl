using Test
using Quickhull

using GeometryBasics
import Random

import DirectQhull

include("utils.jl")

function hullcompare(my_hull, qhull_hull, only_test_verts=false)
    verts_match = Set(vertices(my_hull)) == Set(qhull_hull.vertices)

    if only_test_verts
        @test verts_match
        return
    end

    my_facets = Set(sort.(facets(my_hull)))
    qhull_facets = reinterpret(eltype(my_facets), qhull_hull.simplices)[:]
    qhull_facets = Set(sort.(qhull_facets))

    facets_match = my_facets == qhull_facets

    @test verts_match && facets_match
end

@testset "predicates" begin
    Random.seed!(999)
    for i = 1:100
        mat, pt = rand(SMatrix{3, 3}), rand(SVector{3})
        @test Quickhull.vol_exact_slow(mat, pt) == Quickhull.vol_exact_adaptive_multifloat(mat, pt)
    end

    for i = 1:100
        mat = rand(SMatrix{3, 3})
        pt = sum(mat, dims=2)[:, 1] / 3
        @test Quickhull.vol_exact_slow(mat, pt) == Quickhull.vol_exact_adaptive_multifloat(mat, pt)
    end
end

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

                    my_hull = quickhull(pts, Quickhull.Options(kernel=Quickhull.HyperplaneKernelExactSIMD))
                    qhull_hull = DirectQhull.ConvexHull(pts)

                    only_test_verts = sample_func == gridbox
                    hullcompare(my_hull, qhull_hull, only_test_verts)
                end
            end
        end

    end
end

@testset "inserts" begin
    Random.seed!(1738)

    pts = [rand(Point3d) for i = 1:1000]
    hull = quickhull(pts)

    for i = 1:1000
        insert!(hull, rand(Point3d))
    end

    hull_rebuilt = quickhull(pts)

    @test Set(vertices(hull)) == Set(vertices(hull_rebuilt))
    @test length(facets(hull)) == length(facets(hull_rebuilt))
end

@testset "kernels" begin
    Random.seed!(1234)
    pts = sampleball(100_000, 3)
    their_hull = DirectQhull.ConvexHull(pts)

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
    @inferred quickhull([rand(Point3f) for i = 1:100])

    hull = quickhull([rand(Point3f) for i = 1:100])
    @inferred points(hull)
    @inferred vertices(hull)
    @inferred facets(hull)

    @inferred delaunay([(rand(), rand()) for i = 1:100])
    @inferred delaunay([rand(SVector{3}) for i = 1:100])
    @inferred delaunay([rand(Point3f) for i = 1:100])

    tri = delaunay([rand(Point2f) for i = 1:100])
    @inferred points(tri)
    @inferred vertices(tri)
    @inferred facets(tri)
end

@testset "robustness" begin
    nf0 = nextfloat(0.0)
    ϵ = eps()

    pts_2d = [Point(-10.0, 0.0), Point(0.0, nf0), Point(10.0, 0.0)]
    @test Set(vertices(quickhull(pts_2d))) == Set([1, 2, 3])

    pts_2d = [Point(-10.0, 0.0), Point(0.0, nf0), Point(0.0, nextfloat(0.0, 2)), Point(10.0, 0.0)]
    @test Set(vertices(quickhull(pts_2d))) == Set([1, 3, 4])
    
    pts_2d = [Point(1.0, 1.0), Point(1.0+ϵ, 1.0), Point(1.0, 1.0+ϵ), Point(1.0+ϵ, 1.0+ϵ)]
    @test Set(vertices(quickhull(pts_2d))) == Set([1, 2, 3, 4])

    h = 1e100
    pts_2d = [Point(h, 1.0), Point(h+eps(h), 1.0), Point(h, 1.0+ϵ), Point(h+eps(h), 1.0+ϵ)]
    @test Set(vertices(quickhull(pts_2d))) == Set([1, 2, 3, 4])

    pts_3d = [Point(1.0, 1.0, 1.0), Point(1.0+ϵ, 1.0, 1.0), Point(1.0, 1.0+ϵ, 1.0), Point(1.0, 1.0, 1.0+ϵ)]
    @test Set(vertices(quickhull(pts_3d))) == Set([1, 2, 3, 4])
end

@testset "parallel" begin
    Random.seed!(1234)

    for D = 2:3
        pts = randn(D, 100_000)

        hull_reg = quickhull(pts)

        subdiv = Quickhull.ParallelSubdivide(chunks=4, levels=2)
        hull_par = quickhull(pts, Quickhull.Options(subdivide=subdiv))

        @test Set(vertexpoints(hull_reg)) == Set(vertexpoints(hull_par))
    end
end

using Aqua
@testset "Aqua.jl" begin
    Aqua.test_all(Quickhull)
end
