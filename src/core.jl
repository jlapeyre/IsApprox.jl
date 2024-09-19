"""
    AbstractApprox

Supertype of types specifying how equality or approximate equality should
be computed. See `Equal`, `Approx`, and `EachApprox`.
"""
abstract type AbstractApprox end

if VERSION < v"1.7"
    const Pairs = Base.Iterators.Pairs
else
    const Pairs = Base.Pairs
end

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
    function EachApprox(kw::Pairs)
        new{typeof(kw)}(kw)
    end
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
    function Approx(kw::Pairs)
        new{typeof(kw)}(kw)
    end
end

Approx(; kws...) = Approx(kws)

"""
    isapprox(x, y, ::Equal)

Return `true` if `x == y`.
"""
Base.isapprox(x, y, ::Equal) = (x == y)

"""
    isapprox(x, y, a::Union{EachApprox, Approx})

Use the definition of approximate equality specified by `a` to determine
if `x` is approximately `y`.
"""
Base.isapprox(x, y, a::Union{EachApprox, Approx}) = isapprox(x, y; a.kw...)

"""
    isapprox(A::AbstractArray, B::AbstractArray, a::EachApprox)

Compute element-wise approximate equality.
"""
function Base.isapprox(A::AbstractArray, B::AbstractArray, approx::EachApprox)
    n1 = length(A)
    n2 = length(B)
    if n1 != n2
        return false
    end
    for (x, y) in zip(A, B)
        if ! isapprox(x, y, approx)
            return false
        end
    end
    return true
end

"""
    UpToPhase <: AbstractApprox

Demands that each pair of elements are approximately equal up to a phase,
that is a number whose absolute value is one.
`kw` are keyword pairs that are forwarded to `isapprox`.
"""
struct UpToPhase{T} <: AbstractApprox
    kw::T
    function UpToPhase(kw::Pairs)
        new{typeof(kw)}(kw)
    end
end
UpToPhase(; kws...) = UpToPhase(kws)

function Base.isapprox(x::Number, y::Number, a::UpToPhase)
    aa = Approx(;a.kw...)
    if isapprox(x, zero(x), aa)
        return isapprox(y, zero(y), aa)
    elseif isapprox(y, zero(y), aa)
        return isapprox(x, zero(x), aa)
    end
    return isunitary(x / y, aa)
end

function Base.isapprox(A::AbstractArray, B::AbstractArray, _app::UpToPhase)
    n1 = length(A)
    n2 = length(B)
    if n1 != n2
        return false
    end
    app = Approx(;_app.kw...)
    seen_non_zero_flag = false
    z = zero(eltype(A)) # TODO use promotion
    for (a, b) in zip(A, B)
        if iszero(a, app)
            !iszero(b, app) && return false
        elseif iszero(b, app)
            !iszero(a, app) && return false
        else
            if ! seen_non_zero_flag
                z = a/b
                isunitary(z, app) || return false
                seen_non_zero_flag = true
            else
                isapprox(a/b, z, app) || return false
            end
        end
    end
    return true
end
