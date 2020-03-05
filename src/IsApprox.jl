module IsApprox

export AbstractApprox, Equal, EachApprox, Approx

abstract type AbstractApprox end

"""
    Equal <: AbstractApprox

Approximate equality test that actually demands strict equality
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
    Approx <: AbstractApprox

Specifies using the legacy `isapprox` interface. For example,
for `AbstractMatrix`, matrix norms are used to test closeness.
`kw` are keyword pairs that are forwarded to `isapprox`.
"""
struct Approx{T} <: AbstractApprox
    kw::T
end
Approx(; kws...) = Approx(kws)

Base.isapprox(::Equal, x, y) = (x == y)

Base.isapprox(a::Union{EachApprox, Approx}, x, y) = isapprox(x, y; pairs(a.kw)...)

# Elementwise approximate equality
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

include("applications.jl")

end # module
