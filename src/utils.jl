presized_array(::Type{T}, n) where {T} = empty!(Vector{T}(undef, n))
constructorof(::Type{T}) where T = Base.typename(T).wrapper

# reinterpret a D x N matrix into a length N vector
# of length D points (static vectors)
function matrix2points(pts::Matrix{T}) where {T}
    D, _ = size(pts)
    return dropdims(reinterpret(SVector{D, T}, pts), dims=1)
end

function pointmatrix(pts, indices::SVector{N, T}) where {N, T}
    D = pointsdim(pts)
    mat = hcat(SVector{D}.(pts[Int.(indices)])...)
    @assert mat isa SMatrix
    return mat
end

joggle(points, amt) = map(points) do pt
    T = eltype(pt)
    pt + T(amt) * (2rand(T) - 1) * eps.(pt)
end

@inline exactify(x) = big(rationalize(x, tol=0))

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
    return x isa Array ? view(x, indices...) : x[indices...]
end

# this is used for types, so it should be inferrable as a constant
pointsdim(pts::AbstractVector) = length(first(pts))
