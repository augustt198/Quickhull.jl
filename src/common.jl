const PointIndex = Int32
const FacetIndex = Int32

presized_array(::Type{T}, n) where {T} = empty!(Vector{T}(undef, n))
constructorof(::Type{T}) where T = Base.typename(T).wrapper

# reinterpret a D x N matrix into a length N vector
# of length D points (static vectors)
function matrix2points(pts::Matrix{T}) where {T}
    D, _ = size(pts)
    return dropdims(reinterpret(SVector{D, T}, pts), dims=1)
end

function pointmatrix(pts, indices::SVector{N, T}) where {N, T}
    D = length(first(pts))
    mat = hcat(SVector{D}.(pts[Int.(indices)])...)
    @assert mat isa SMatrix
    return mat
end

function joggle(points, amt)
    D = length(points[1])

    map(points) do pt
        s = 2rand() - 1
        SVector{D}(x + amt * s * eps(x) for x in pt)
    end
end

@inline exactify(x) = big.(rationalize.(x, tol=0))
