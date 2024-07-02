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

joggle(points, amt) = map(points) do pt
    pt + amt * (2rand() - 1) * eps.(x)
end

@inline exactify(x) = big.(rationalize.(x, tol=0))

droplast(x; dims=(ndims(x),)) = _droplast(x, dims)

Base.@constprop :aggressive function _droplast(x, dims)
    indices = ntuple(ndims(x)) do i
        if i ∈ dims
            if x isa StaticArray
                N = size(x, i)
                return SVector{N-1}(i for i = 1:(N-1))
            else
                return axes(x, i)[begin:end-1]
            end
        else
            return (:)
        end
    end
    view(x, indices...)
end
