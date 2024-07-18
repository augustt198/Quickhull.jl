using LinearAlgebra
using StaticArrays

function samplebox(N, D)
    return rand(D, N)
end

function samplesphere(N, D)
    pts = randn(D, N)
    return pts ./ norm.(eachcol(pts))'
end

function sampleball(N, D)
    r = rand(N) .^ (1/D)
    return samplesphere(N, D) .* r'
end

function gridbox(N, D)
    ax = range(-1, 1, length=N)
    SV = SVector{D, eltype(ax)}
    pts = vec(SV.(Iterators.product(ntuple(_ -> ax, D)...)))
    return reduce(hcat, pts)
end

hypercube(D) = [float((idx >> i) & 1) for idx=0:(2^D-1), i=0:(D-1)]
