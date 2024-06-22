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
            #expr = ifelse(isodd(i) || perm, :($expr + $sym * $ex), :($expr - $sym * $ex))
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
        return :(abs($a * $c) + abs($b * $d))
    end

    expr = nothing
    for (i, sym) in enumerate(rows[1])
        newrows = [row[(1:i-1) ∪ (i+1:n)] for row in rows[2:n]]

        ex = symbolic_permanent_exactinit(newrows)
        if isnothing(expr)
            expr = :($sym * $ex)
        else
            expr = :($expr + $sym * $ex)
        end
    end

    return expr
end

const ε = Polynomial([0, 1])

function symbolic_det_relative_error(N, base=ε, ε_val=big(rationalize(eps())))
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

function detn_exact_slow(mat::SMatrix{N, N, T}) where {N, T}
    d = det(exactify(mat))
    return copysign(T(d), d)
end

function vol_exact_slow(mat::SMatrix{N, N, T}, pt::SVector{N, T}) where {N, T}
    d = det(exactify(mat) .- exactify(pt))
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
    ε_val = big(rationalize(eps(T)))
    poly = symbolic_det_relative_error(N, ε, ε_val)
    rel = (poly * (1 + ε))(big(rationalize(eps())))
    rel_f = Float64(rel)

    ex = Expr(:block,
        :(M = mat .- pt),
        assignments...,
        quote
            _det = $res_det
            _perm = $res_perm
            if abs(_det) > $rel_f * _perm
                return _det
            else
                vol_exact_slow(mat, pt)
            end
        end)
    ex
end
