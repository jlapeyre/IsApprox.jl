"""
    AbstractApprox

Supertype of types specifying how equality or approximate equality should
be computed. See `Equal`, `Approx`, and `EachApprox`.
"""
abstract type AbstractApprox end

"""
    Equal <: AbstractApprox

Approximate equality test that actually demands strict equality.
"""
struct Equal <: AbstractApprox
end

"""
    EachApprox <: AbstractApprox

Demands that each pair of elements are approximately equal.
`kw` are keyword pairs that are forwarded to `isapprox`.
"""
struct EachApprox{T} <: AbstractApprox
    kw::T
end
EachApprox(; kws...) = EachApprox(kws)

"""
   Approx(; kws...) <: AbstractApprox

Specifies using the legacy `isapprox` interface. For example,
for `AbstractMatrix`, matrix norms are used to test closeness.
`kw` are keyword pairs that are forwarded to `isapprox`.
For example, `Approx(atol=1e-9)`.
"""
struct Approx{T} <: AbstractApprox
    kw::T
end
Approx(; kws...) = Approx(kws)

"""
    isapprox(::Equal, x, y)

Return `true` if `x == y`.
"""
Base.isapprox(::Equal, x, y) = (x == y)

"""
    isapprox(a::Union{EachApprox, Approx}, x, y)

Use the definition of approximate equality specified by `a` to determine
if `x` is approximately `y`.
"""
Base.isapprox(a::Union{EachApprox, Approx}, x, y) = isapprox(x, y; pairs(a.kw)...)

"""
    isapprox(a::EachApprox, A::AbstractArray, B::AbstractArray)

Compute element-wise approximate equality.
"""
function Base.isapprox(a::EachApprox, A::AbstractArray, B::AbstractArray)
    n1 = length(A)
    n2 = length(B)
    if n1 != n2
        return false
    end
    for (x, y) in zip(A, B)
        if ! isapprox(a, x, y)
            return false
        end
    end
    return true
end
