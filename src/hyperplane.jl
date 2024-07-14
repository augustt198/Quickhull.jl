abstract type HyperplaneKernel end

struct Hyperplane{D, I <: Integer, K <: HyperplaneKernel}
    point_indices::SVector{D, I}
    kernel::K
end

function Hyperplane(pis::SVector{D, I}, points, ::Type{K}) where {D, I, K}
    Hyperplane(pis, make_kernel(K, pis, points))
end

function hyperplane_towards(plane, pt, pts)
    return hyperplane_dist(plane, pt, pts) > 0 ? plane : hyperplane_invert(plane)
end

function hyperplane_awayfrom(plane, pt, pts)
    return hyperplane_dist(plane, pt, pts) < 0 ? plane : hyperplane_invert(plane)
end


# Different kernels for hyperplane calculations:

# An inexact kernel that uses the typical dot product representation
# of a plane.
struct HyperplaneKernelInexact{D, T <: Number} <: HyperplaneKernel
    normal::SVector{D, T}
    offset::T
end

function make_kernel(::Type{K}, point_indices::SVector{D, I}, pts::AbstractVector) where {D, I, K <: HyperplaneKernelInexact}
    pts_T = pointmatrix(pts, point_indices)
    mat = vcat(hcat(pts_T, @SVector zeros(D)), (@SVector ones(D+1))')

    ns = qr(mat).Q[:,end] # nullspace

    normal, offset = droplast(ns), ns[end]
    mag = norm(normal)

    HyperplaneKernelInexact(normal / mag, offset / mag)
end

function hyperplane_dist(plane::Hyperplane{D, I, K}, pt, pts) where {D, I, K <: HyperplaneKernelInexact}
    k = plane.kernel
    dot(k.normal, pt) + k.offset
end

function hyperplane_invert(plane::Hyperplane{D, I, K}) where {D, I, K <: HyperplaneKernelInexact}
    k = plane.kernel
    Hyperplane(plane.point_indices, K(-k.normal, -k.offset))
end

# An exact kernel that stores the coordinates of the points
# defining the plane, so an exact determinant can be taken
struct HyperplaneKernelExact_A{D, T, D2, L} <: HyperplaneKernel
    mat::SMatrix{D, D2, T, L}
    sign::T
end

function make_kernel(::Type{K}, point_indices::SVector{D, I}, pts::AbstractVector) where {D, I, K <: HyperplaneKernelExact_A}
    mat = pointmatrix(pts, point_indices)
    HyperplaneKernelExact_A(mat, one(eltype(mat)))
end

function hyperplane_dist(plane::Hyperplane{D, I, K}, pt, pts) where {D, I, K <: HyperplaneKernelExact_A}
    k = plane.kernel

    pt′ = SVector{length(pt)}(pt)

    if size(k.mat, 1) == size(k.mat, 2)
        return vol_exact(k.mat, pt′) * k.sign
    else
        M_big = exactify(k.mat) .- exactify(pt′)
        det_big = det(M_big' * M_big) # this is always nonnegative
        return sqrt(eltype(k.mat)(det_big)) * k.sign
    end
end

function hyperplane_invert(plane::Hyperplane{D, I, K}) where {D, I, K <: HyperplaneKernelExact_A}
    k = plane.kernel
    Hyperplane(plane.point_indices, K(k.mat, -k.sign))
end



@generated function _det_expanded(mat::SMatrix{D, D, T, D2}) where {D, T, D2}
    mat2 = [[:(mat[$i, $j]) for j = 1:D] for i = 1:D]
    return symbolic_det(mat2)
end

@generated function _perm_expanded(mat::SMatrix{D, D, T, D2}) where {D, T, D2}
    mat2 = [[:(abs(mat[$i, $j])) for j = 1:D] for i = 1:D]
    return  symbolic_permanent(mat2)
end

# this is broken right now
struct HyperplaneKernelExact_C{T, Dp1} <: HyperplaneKernel
    minors::SVector{Dp1, T}
    minor_perms::SVector{Dp1, T}
    sign::T
end

function make_kernel(::Type{K}, point_indices::SVector{D, I}, pts::AbstractVector) where {D, I, K <: HyperplaneKernelExact_C}
    mat = pointmatrix(pts, point_indices)'

    T = eltype(mat)

    mat′ = hcat(mat, ones(SVector{D, T}))
    function minor(i)
        # drop ith index
        ind = SVector{D}(j >= i ? j+1 : j for j = 1:D)
        cut = mat′[:, ind]
        _det_expanded(cut)
    end

    function minor_perm(i)
        ind = SVector{D}(j >= i ? j+1 : j for j = 1:D)
        cut = mat′[:, ind]
        _perm_expanded(cut)
    end

    minors = SVector{D+1, T}(minor(i) for i = 1:D+1)
    minor_perms = SVector{D+1, T}(minor_perm(i) for i = 1:D+1)

    HyperplaneKernelExact_C(minors, minor_perms, one(T))
end

@inline function det_using_minors(minors, minor_perms, pt)
    D = length(pt)
    Dp1 = length(minors)

    d = p = zero(eltype(minors))
    for i = 1:D
        d_mul = pt[i] * minors[i]
        p_mul = abs(pt[i]) * minor_perms[i]

        if i == 1
            d = d_mul
            p = p_mul
        else
            d += isodd(i) ? d_mul : -d_mul
            p += p_mul
        end
    end

    d += isodd(Dp1) ? minors[Dp1] : -minors[Dp1]
    p += abs(minor_perms[Dp1])

    return d, p
end

@generated function det_using_minors_relerror2(::Val{D}, ::Type{T}) where {D, T}
    ε = Polynomial([0, 1])
    ε_val = big(rationalize(eps(T)))
    mul_term_err = ε + (1 + ε)*symbolic_det_relative_error(D-1, 0, ε_val)

    cum_err = mul_term_err
    for i = 2:D
        if cum_err(ε_val) > mul_term_err(ε_val)
            cum_err = ε + (1 + ε)*cum_err
        else
            cum_err = ε + (1 + ε)*mul_term_err
        end
    end
    # cum_err = ε + (1 + ε)*(cum_err + ε) # huhhhh

    rel_err_arb = cum_err(ε_val)
    final = (1 + ε_val) * rel_err_arb

    return T(final)
end

function hyperplane_dist(plane::Hyperplane{D, I, HyperplaneKernelExact_C{T, Dp1}}, pt, pts) where {D, I, Dp1, T}
    k = plane.kernel

    minors, minor_perms = k.minors, k.minor_perms
    d = p = zero(eltype(minors))
    for i = 1:D
        if i == 1
            d = pt[i] * minors[i]
            p = abs(pt[i]) * minor_perms[i]
        else
            sgn = isodd(i) ? 1 : -1
            d = muladd(sgn * pt[i], minors[i], d)
            p = muladd(abs(pt[i]), minor_perms[i], p)
        end
    end

    d += isodd(Dp1) ? minors[Dp1] : -minors[Dp1]
    p += minor_perms[Dp1]

    if abs(d) > det_using_minors_relerror2(Val(Dp1), T) * p
        return d * k.sign
    else
        mat = pointmatrix(pts, plane.point_indices)
        pt′ = SVector{D}(pt)
        return vol_exact_slow(mat, pt′) * -k.sign
    end
end

function hyperplane_invert(plane::Hyperplane{D, I, K}) where {D, I, K <: HyperplaneKernelExact_C}
    k = plane.kernel
    Hyperplane(plane.point_indices, K(k.minors, k.minor_perms, -k.sign))
end

struct HyperplaneKernelExactSIMD{T, Dp1} <: HyperplaneKernel
    cofactor_pairs::SVector{Dp1, SIMD.Vec{2, T}}
    sign::T
end

function make_kernel(::Type{K}, point_indices::SVector{D, I}, pts::AbstractVector) where {D, I, K <: HyperplaneKernelExactSIMD}
    pt_D = length(first(pts))
    if D != pt_D
        # fallback when not simplex
        return make_kernel(HyperplaneKernelExact_A, point_indices, pts)
    end

    mat = pointmatrix(pts, point_indices)'
    T = eltype(mat)

    mat′ = hcat(mat, ones(SVector{D, T}))
    function cofactor(i)
        # drop ith index
        ind = SVector{D}(j >= i ? j+1 : j for j = 1:D)
        cut = mat′[:, ind]
        (-1)^(i+1) * _det_expanded(cut)
    end

    function cofactor_perm(i)
        ind = SVector{D}(j >= i ? j+1 : j for j = 1:D)
        cut = mat′[:, ind]
        _perm_expanded(cut)
    end

    cofactor_pairs = SVector{D+1, SIMD.Vec{2, T}}(SIMD.Vec{2, T}((cofactor(i), cofactor_perm(i))) for i = 1:D+1)
    HyperplaneKernelExactSIMD(cofactor_pairs, one(T))
end

function hyperplane_dist(plane::Hyperplane{D, I, HyperplaneKernelExactSIMD{T, Dp1}}, pt, pts) where {D, I, Dp1, T}
    k = plane.kernel
    cps = k.cofactor_pairs

    pt_simd = SIMD.Vec{D, T}((pt...,))
    pt_abs_simd = abs(pt_simd)
    dp = SIMD.Vec{2, T}(0)
    for i = 1:D
        coord = SIMD.Vec{2, T}((pt_simd[i], pt_abs_simd[i]))
        if i == 1
            dp = muladd(coord, cps[1], cps[Dp1])
        else
            dp = muladd(coord, cps[i], dp)
        end
    end

    d, p = dp[1], dp[2]

    if abs(d) > det_using_minors_relerror2(Val(Dp1), T) * p
        return d * k.sign
    else
        mat = pointmatrix(pts, plane.point_indices)
        pt′ = SVector{D}(pt)
        return vol_exact_multifloat(mat, pt′) * k.sign * (-1)^D # not sure why
    end
end

function hyperplane_invert(plane::Hyperplane{D, I, K}) where {D, I, K <: HyperplaneKernelExactSIMD}
    k = plane.kernel
    Hyperplane(plane.point_indices, K(k.cofactor_pairs, -k.sign))
end
