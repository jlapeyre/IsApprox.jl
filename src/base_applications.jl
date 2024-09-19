### The functions in this file extend functions of the same names in `Base` and `LinearAlgebra`.
### It may be necessary to instead replace the `Base` and `LinearAlgebra` functions
### if `IsApprox` is widely adopted.

using LinearAlgebra: Hermitian, Symmetric, HermOrSym

import Base: isone, iszero, isreal, isinteger
import LinearAlgebra: ishermitian, issymmetric, istriu, istril, isbanded, isdiag

### isone, iszero

Base.isone(x, approx_test::AbstractApprox) = isapprox(x, one(x), approx_test)
Base.iszero(x, approx_test::AbstractApprox) = isapprox(x, zero(x), approx_test)
iszero(approx_test::AbstractApprox) = x -> iszero(x, approx_test)
iszero(x::AbstractArray, approx_test::AbstractApprox) = all(iszero(approx_test), x)
isone(x::BigInt, ::Equal) = Base.isone(x)
iszero(x::BigInt, ::Equal) = Base.iszero(x)

isone(A::StridedMatrix, approx_test::Approx) = isapprox(A, one(A), approx_test)
iszero(A::StridedMatrix, approx_test::Approx) = isapprox(A, zero(A), approx_test)

# dense.jl
const ISONE_CUTOFF = 2^21 # 2M

function isone(A::StridedMatrix, approx_test::AbstractApprox)
    m, n = size(A)
    m != n && return false # only square matrices can satisfy x == one(x)
    if sizeof(A) < ISONE_CUTOFF
        _isone_triacheck(A, m, approx_test)
    else
        _isone_cachefriendly(A, m, approx_test)
    end
end

@inline function _isone_triacheck(A::StridedMatrix, m::Int, approx_test::AbstractApprox)
    @inbounds for i in 1:m, j in i:m
        if i == j
            isone(A[i,i], approx_test) || return false
        else
            iszero(A[i,j], approx_test) && iszero(A[j,i], approx_test) || return false
        end
    end
    return true
end

# Inner loop over rows to be friendly to the CPU cache
@inline function _isone_cachefriendly(A::StridedMatrix, m::Int, approx_test::AbstractApprox)
    @inbounds for i in 1:m, j in 1:m
        if i == j
            isone(A[i,i], approx_test) || return false
        else
            iszero(A[j,i], approx_test) || return false
        end
    end
    return true
end

### ishermitian, issymmetric

# The call `ishermitian(A::AbstractMatrix, B::AbstractMatrix) lowers
# to exactly the same code as that in LinearAlgebra
function ishermitian(A::AbstractMatrix, approx_test::AbstractApprox)
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    for i = indsn, j = i:last(indsn)
        if ! isapprox(A[i,j], adjoint(A[j,i]), approx_test)
            return false
        end
    end
    return true
end

# This method uses the isapprox interface, which compares using a norm.
ishermitian(A::AbstractMatrix, approx_test::Approx) = isapprox(A, adjoint(A), approx_test)

ishermitian(x::Number, approx_test::AbstractApprox) = isapprox(x, conj(x), approx_test)
ishermitian(A::Hermitian, ::AbstractApprox) = true
ishermitian(A::Hermitian, ::Approx) = true
ishermitian(A::Symmetric{<:Real}, ::AbstractApprox) = true
ishermitian(A::Symmetric{<:Real}, ::Approx) = true
ishermitian(A::Symmetric{<:Complex}, approx_test::AbstractApprox) = isreal(A, approx_test)
ishermitian(A::Symmetric{<:Complex}, approx_test::Approx) = isreal(A, approx_test)
issymmetric(A::Hermitian{<:Real}, ::AbstractApprox) = true
issymmetric(A::Hermitian{<:Complex}, approx_test::AbstractApprox) = isreal(A, approx_test)
issymmetric(A::Symmetric, ::AbstractApprox) = true

issymmetric(A::AbstractMatrix{<:Real}, approx_test::AbstractApprox) =
    ishermitian(A, approx_test)

# Copied from LinearAlgebra. Why does the iteration over i differ slightly from ishermitian ?
function issymmetric(A::AbstractMatrix, approx_test::AbstractApprox)
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    for i = first(indsn):last(indsn), j = (i):last(indsn)
        if ! isapprox(A[i,j], transpose(A[j,i]), approx_test)
            return false
        end
    end
    return true
end

issymmetric(x::Number, approx_test::AbstractApprox) = isapprox(x, x, approx_test)

### isreal

# complex.jl
isreal(x::Real, approx_test::AbstractApprox) = true
isreal(z::Complex, approx_test::AbstractApprox) = isapprox(real(z), z, approx_test)
# Old way
#isreal(z::Complex, approx_test::AbstractApprox) = iszero(imag(z), approx_test)

isreal(x::AbstractArray{<:Real}, approx_test::AbstractApprox) = true

isreal(approx_test::AbstractApprox) = x -> isreal(x, approx_test)
isreal(x::AbstractArray, approx_test::AbstractApprox) = all(isreal(approx_test),x)

isreal(A::HermOrSym{<:Real}, approx_test::AbstractApprox) = true
function isreal(A::HermOrSym, approx_test::AbstractApprox)
    n = size(A, 1)
    @inbounds if A.uplo == 'U'
        for j in 1:n
            for i in 1:(j - (A isa Hermitian))
                if !isreal(A.data[i,j], approx_test)
                    return false
                end
            end
        end
    else
        for j in 1:n
            for i in (j + (A isa Hermitian)):n
                if !isreal(A.data[i,j], approx_test)
                    return false
                end
            end
        end
    end
    return true
end

### isinteger

# This must be changed if this is integrated into Base
# isinteger(x) = isinteger(x, Equal())
# We use union or explicit types to avoid method ambiguity. Is there a way around this ?
isinteger(x::BigFloat, ::Equal) = Base.isinteger(x)
isinteger(x::Rational, ::Equal) = Base.isinteger(x)

# number.jl
isinteger(x::Integer, ::AbstractApprox) = true

# floatfuncs.jl
## The original is x - trunc(x) == 0. So this implementation might differ
## from the stock version of isinteger for some user's subtype of AbstractFloat.
## We choose to this implementation because the default relative tolerance is
## reasonable. That is, `isapprox(approx_test, x - trunc(x), 0)` requires
## specifying `atol`.
isinteger(x::AbstractFloat, approx_test::AbstractApprox) = isapprox(x, trunc(x), approx_test)

# complex.jl
isinteger(z::Complex, approx_test::AbstractApprox) =
    isreal(z, approx_test) && isinteger(real(z), approx_test)

# TODO: Need to think about difference between real(z) and abs(z) regarding tolerance in all
# methods for complex numbers (not just isinteger)
isinteger(z::Complex, approx_test::UpToPhase) = isinteger(abs(z), approx_test)

isinteger(x::Rational, approx_test::AbstractApprox) = isinteger(float(x), approx_test)

### istriu, istril, isbanded, isdiag

# For compatibility
_require_one_based_indexing(A...) = !Base.has_offset_axes(A...) || throw(ArgumentError("offset arrays are not supported but got an array with index other than 1"))

# TODO: Why is approx a kw arg here and below?
function istriu(A::AbstractMatrix, k::Integer, approx::AbstractApprox)
    _require_one_based_indexing(A)
    m, n = size(A)
    for j in 1:min(n, m + k - 1)
        for i in max(1, j - k + 1):m
            iszero(A[i, j], approx) || return false
        end
    end
    return true
end
istriu(::Number, ::AbstractApprox) = true

function istril(A::AbstractMatrix, k::Integer, approx::AbstractApprox)
    _require_one_based_indexing(A)
    m, n = size(A)
    for j in max(1, k + 2):n
        for i in 1:min(j - k - 1, m)
            iszero(A[i, j], approx) || return false
        end
    end
    return true
end
istril(::Number, ::AbstractApprox) = true

isbanded(A::AbstractMatrix, kl::Integer, ku::Integer, approx::AbstractApprox) =
    istriu(A, kl, approx) && istril(A, ku, approx)

isdiag(A::AbstractMatrix, approx::AbstractApprox) = isbanded(A, 0, 0, approx)
isdiag(x::Number, approx::AbstractApprox) = true
isdiag(A::HermOrSym, approx::AbstractApprox) =
    isdiag(A.uplo == 'U' ? LinearAlgebra.UpperTriangular(A.data) :
    LinearAlgebra.LowerTriangular(A.data), approx)
