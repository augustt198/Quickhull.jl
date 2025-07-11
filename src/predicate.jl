using Polynomials

function symbolic_det(rows, perm=false)
    n = length(rows)
    n == 1 && return rows[1][1]

    expr = nothing
    for (i, sym) in enumerate(rows[1])
        newrows = [row[(1:i-1) ∪ (i+1:n)] for row in rows[2:n]]

        ex = symbolic_det(newrows, perm)
        if isnothing(expr)
            expr = :($sym * $ex)
        else
            expr = ifelse(isodd(i) || perm, :(muladd($sym, $ex, $expr)), :(muladd(-$sym, $ex, $expr)))
        end
    end

    return expr
end

# is this why they called it floating point "permanent" ... 
symbolic_permanent(rows) = symbolic_det(rows, true)

# i.e. assume the matrix starts error free
function symbolic_permanent_exactinit(rows)
    n = length(rows)
    if n == 2
        a, b, c, d = rows[1][1], rows[1][2], rows[2][1], rows[2][2]
        return :(abs($a * $d) + abs($b * $c))
    end

    expr = nothing
    for (i, sym) in enumerate(rows[1])
        newrows = [row[(1:i-1) ∪ (i+1:n)] for row in rows[2:n]]

        ex = symbolic_permanent_exactinit(newrows)
        if isnothing(expr)
            expr = :($sym * $ex)
        else
            expr = :(muladd($sym, $ex, $expr))
        end
    end

    return expr
end

const ε = Polynomial([0, 1])

function symbolic_det_relative_error(N, base=ε, ε_val=exactify(eps()))
    N == 1 && return base

    expr = nothing
    for i = 1:N
        δ =  symbolic_det_relative_error(N - 1)
        mul = ε + (1 + ε)*(base + δ + base*δ)
        if isnothing(expr)
            expr = mul
        else
            if expr(ε_val) > mul(ε_val)
                expr = ε + (1 + ε)*expr
            else
                expr = ε + (1 + ε)*mul
            end
        end
    end

    return expr
end

@generated function _det_expanded(mat::SMatrix{D, D, T, D2}) where {D, T, D2}
    mat2 = [[:(mat[$i, $j]) for j = 1:D] for i = 1:D]
    return symbolic_det(mat2)
end

@generated function _perm_expanded(mat::SMatrix{D, D, T, D2}) where {D, T, D2}
    mat2 = [[:(abs(mat[$i, $j])) for j = 1:D] for i = 1:D]
    return  symbolic_permanent(mat2)
end

function detn_exact_slow(mat::SMatrix{N, N, T}) where {N, T}
    d = det(exactify.(mat))
    return copysign(T(d), d)
end

function vol_exact_slow(mat::SMatrix{N, N, T}, pt::SVector{N, T}) where {N, T}
    d = _det_expanded(exactify.(mat) .- exactify.(pt))
    return copysign(T(d), d)
end

@generated function vol_exact(mat::SMatrix{N, N, T}, pt::SVector{N, T}) where {N, T}
    assignments = [
        :($(Symbol(:a, i, j)) = M[$i, $j]) for i = 1:N for j = 1:N
    ]

    mat1 = [[Symbol(:a, i, j) for j = 1:N] for i = 1:N]
    mat2 = [[:(abs( $(Symbol(:a, i, j)) )) for j = 1:N] for i = 1:N]
    res_det = symbolic_det(mat1)
    res_perm = symbolic_permanent(mat2)
    ε_val = exactify(eps(T))
    poly = symbolic_det_relative_error(N, ε, ε_val)
    rel = (poly * (1 + ε))(ε_val)
    rel_f = T(rel)

    ex = Expr(:block,
        :(M = mat .- pt),
        assignments...,
        quote
            _det = $res_det
            _perm = $res_perm
            if abs(_det) > $rel_f * _perm
                return _det
            else
                vol_exact_adaptive_multifloat(mat, pt)
            end
        end)
    ex
end


# A MultiFloat that automatically extends its size so arithmetic is always exact.
# Could still be optimized more
struct ExactMultiFloat{T <: AbstractFloat, D} <: AbstractFloat
    x::MultiFloat{T, D}

    ExactMultiFloat(x::MultiFloat{T, D}) where {T <: AbstractFloat, D} = new{T, D}(x)
end

for op in (:+, :-)
    @eval function Base.$op(a::ExactMultiFloat{T, D1}, b::ExactMultiFloat{T, D2}) where {T, D1, D2}
        D = D1 + D2
        a′ = MultiFloat{T, D}(a.x)
        b′ = MultiFloat{T, D}(b.x)
        ExactMultiFloat($op(a′, b′))
    end
end

Base.:-(a::ExactMultiFloat) = ExactMultiFloat(-a.x)

function Base.:*(a::ExactMultiFloat{T, D1}, b::ExactMultiFloat{T, D2}) where {T, D1, D2}
    D = 2 * D1 * D2
    a′ = MultiFloat{T, D}(a.x)
    b′ = MultiFloat{T, D}(b.x)
    ExactMultiFloat(a′ * b′)
end

Base.muladd(a::ExactMultiFloat, b::ExactMultiFloat, c::ExactMultiFloat) = a * b + c

function vol_exact_multifloat(mat::SMatrix{N, N, T}, pt::SVector{N, T}) where {N, T}
    mat = ExactMultiFloat.(MultiFloat{T, 1}.(mat))
    pt =  ExactMultiFloat.(MultiFloat{T, 1}.(pt))

    M = mat .- pt
    d = _det_expanded(M).x
    return copysign(T(d), d)
end

using MacroTools

macro bodyinbounds(expr)
    MacroTools.postwalk(expr) do ex
        if ex isa Expr && ex.head == :function
            newbody = :(@inbounds $(ex.args[2]))
            return Expr(:function, ex.args[1], newbody)
        end
        return ex
    end
end

# Implements multi-float arithmetic in an exact way. There are
# basically three lengths in play here, which may be confusing.
# `N` is the underlying length of the static vector of terms. `length`
# is the number of terms that are actually occupied (nonzero), stored
# so expansion_sum, compress_expansion, etc. don't do more work than
# necessary. Nmax limits the maximum N value that results from doing
# arithmetic on AdaptiveMultiFloats, i.e. N adjusts automatically up
# to Nmax. The purpose of this is to avoid blow-up in the number of terms.
# For instance, the result of the orient3d test requires 192 terms in the
# absolute worst case - this leads to excruciating compile times. Usually,
# only a fraction of those 192 terms are needed, so we can limit Nmax to
# something like 8 or 16. If the result of an arithmetic operation can't
# be represented exactly for the given Nmax, LossOfPrecisionException
# will be raised. This should be extremely rare, so it can be handled by
# using slow arbitrary-precision arithmetic.
struct AdaptiveMultiFloat{T <: AbstractFloat, N, Nmax}
    terms::SVector{N, T}
    length::Int
end

expansion_sum(e::SVector{N1, T}, f::SVector{N2, T}, elen, flen) where {N1, N2, T} = expansion_sum(e, f, elen, flen, Val{N1 + N2}())

@bodyinbounds function expansion_sum(e::SVector{N1, T}, f::SVector{N2, T}, elen, flen, ::Val{Nnew}) where {N1, N2, T, Nnew}
    res = zero(SVector{Nnew, T})
    for j = 1:elen
        res = setindex(res, e[j], j)
    end

    for i = 1:flen
        hi = f[i]
        lo = zero(hi)
        for j = 1:elen
            hi, lo = MultiFloats.two_sum(res[i + j - 1], hi)
            res = setindex(res, lo, i + j - 1)
        end
        res = setindex(res, hi, i + elen)
    end

    return res
end

@bodyinbounds function compress_expansion(e::SVector, elen)
    g = zero(e)
    Q = e[elen]
    bottom = elen

    for i = (elen-1):-1:1
        Q, q = MultiFloats.fast_two_sum(Q, e[i])
        if !iszero(q)
            g = setindex(g, Q, bottom)
            bottom -= 1
            Q = q
        end
    end
    g = setindex(g, Q, bottom)

    h = zero(e)
    top = 1
    for i = (bottom+1):elen
        Q, q = MultiFloats.fast_two_sum(g[i], Q)
        if !iszero(q)
            h = setindex(h, q, top)
            top += 1
        end
    end
    h = setindex(h, Q, top)
    return h, top
end

@bodyinbounds function zero_elim(e::SVector)
    res = zero(e)
    idx = 1
    for (i, x) in enumerate(e)
        if !iszero(x)
            res = setindex(res, x, idx)
            idx += 1
        end
    end

    length = idx - 1
    return res, length
end

expansion_sum(rand(SVector{4}), rand(SVector{4}), 3, 3)

struct LossOfPrecisionException <: Exception
end

@bodyinbounds function Base.:+(a::AdaptiveMultiFloat{T, N1, Nmax}, b::AdaptiveMultiFloat{T, N2, Nmax}) where {T, N1, N2, Nmax}
    summed = expansion_sum(a.terms, b.terms, a.length, b.length)
    terms, len = compress_expansion(summed, a.length + b.length)
    if len <= Nmax
        nterms = min(length(summed), Nmax)
        terms_trunc = SVector{nterms}(terms[i] for i = 1:nterms)

        AdaptiveMultiFloat{T, nterms, Nmax}(terms_trunc, len)
    else
        throw(LossOfPrecisionException())
    end
end

@bodyinbounds function scale_expansion(e::SVector{N, T}, elen, b) where {N, T}
    h = zero(SVector{2*N, T})
    Q, hnew = MultiFloats.two_prod(e[1], b)
    h = setindex(h, hnew, 1)

    for i = 2:elen
        Ti, ti = MultiFloats.two_prod(e[i], b)
        Q, hnew = MultiFloats.two_sum(Q, ti)
        h = setindex(h, hnew, 2*i - 2)
        Q, hnew = MultiFloats.fast_two_sum(Ti, Q)
        h = setindex(h, hnew, 2*i - 1)
    end
    h = setindex(h, Q, 2*elen)

    return h
end

@bodyinbounds function Base.:*(a::AdaptiveMultiFloat{T, N1, Nmax}, b::AdaptiveMultiFloat{T, N2, Nmax}) where {T, N1, N2, Nmax}
    prod_sum = zero(SVector{2 * N1 * N2, T})

    prod_sum_length = 0
    for a_idx = 1:a.length
        a_term = a.terms[a_idx]
        
        prod = scale_expansion(b.terms, b.length, a_term)
        prod_sum = expansion_sum(prod_sum, prod, prod_sum_length, length(prod), Val{length(prod_sum)}())
        prod_sum_length += 2 * N2
    end

    terms, len = compress_expansion(prod_sum, length(prod_sum))
    if len <= Nmax
        nterms = min(length(prod_sum), Nmax)
        terms_trunc = SVector{nterms}(terms[i] for i = 1:nterms)

        AdaptiveMultiFloat{T, nterms, Nmax}(terms_trunc, len)
    else
        throw(LossOfPrecisionException())
    end
end

function Base.:-(a::AdaptiveMultiFloat{T, N, Nmax}) where {T, N, Nmax}
    AdaptiveMultiFloat{T, N, Nmax}(-a.terms, a.length)
end

Base.:-(a::AdaptiveMultiFloat, b::AdaptiveMultiFloat) = a + -b

Base.muladd(a::AdaptiveMultiFloat, b::AdaptiveMultiFloat, c::AdaptiveMultiFloat) = a*b + c

function approximate(a::AdaptiveMultiFloat)
    if a.length == 0
        return zero(eltype(a.terms))
    else
        return @inbounds a.terms[a.length]
    end
end

function AdaptiveMultiFloat(mf::MultiFloat{T, N}) where {T, N}
    AdaptiveMultiFloat{T, N, N}(reverse(mf._limbs), N)
end

function AdaptiveMultiFloat(mf::MultiFloat{T, N}, ::Val{Nmax}) where {T, N, Nmax}
    AdaptiveMultiFloat{T, N, Nmax}(reverse(mf._limbs), N)
end

function AdaptiveMultiFloat(x::AbstractFloat, ::Val{Nmax}) where {Nmax}
    AdaptiveMultiFloat{typeof(x), 1, Nmax}(SVector{1}(x), 1)
end

function vol_exact_adaptive_multifloat(mat::SMatrix{N, N, T}, pt::SVector{N, T}) where {N, T}
    Nmax = N <= 3 ? Val{8}() : Val{16}()
    mat = AdaptiveMultiFloat.(mat, Nmax)
    pt =  AdaptiveMultiFloat.(pt, Nmax)

    try
        M = mat .- pt
        d = _det_expanded(M)
        return approximate(d)
    catch e
        if e isa LossOfPrecisionException
            return vol_exact_slow(mat, pt)
        else
            rethrow(e)
        end
    end
end
