### The functions in this file extend functions of the same names in `Base` and `LinearAlgebra`.
### It may be necessary to instead replace the `Base` and `LinearAlgebra` functions
### if `IsApprox` is widely adopted.

using LinearAlgebra: Hermitian, Symmetric, HermOrSym

import Base: isone, iszero, isreal, isinteger
import LinearAlgebra: ishermitian, issymmetric, istriu, istril, isbanded, isdiag

### isone, iszero

"""
    isone(x, approx::AbstractApprox)

Return true if x is approximately equal to one, according to approx.
- Equal: exact equality.
- Approx: forwards to isapprox(x, one(x); approx.kw...), typically norm-based for arrays.
- EachApprox: element-wise approximate equality for arrays.
Optimized paths exist for StridedMatrix to avoid excess allocation.
"""
Base.isone(x, approx_test::AbstractApprox) = isapprox(x, one(x), approx_test)

"""
    iszero(x, approx::AbstractApprox)

Return true if x is approximately equal to zero, according to approx.
- Equal: exact equality.
- Approx: forwards to isapprox(x, zero(x); approx.kw...), typically norm-based for arrays.
- EachApprox: element-wise approximate equality for arrays.
Use iszero(approx)(x) to obtain a predicate for broadcasting/folds.
"""
Base.iszero(x, approx_test::AbstractApprox) = isapprox(x, zero(x), approx_test)

iszero(approx_test::AbstractApprox) = x -> iszero(x, approx_test)
function iszero(x::AbstractArray, approx::AbstractApprox)
    @inbounds for xi in x
        iszero(xi, approx) || return false
    end
    return true
end
# iszero(x::AbstractArray, approx_test::AbstractApprox) = all(iszero(approx_test), x)

isone(x::BigInt, ::Equal) = Base.isone(x)
iszero(x::BigInt, ::Equal) = Base.iszero(x)

isone(A::StridedMatrix, approx_test::Approx) = isapprox(A, one(A), approx_test)
iszero(A::StridedMatrix, approx_test::Approx) = isapprox(A, zero(A), approx_test)

# dense.jl
const ISONE_CUTOFF = 2^21 # 2M

function isone(A::StridedMatrix, approx_test::AbstractApprox)
    m, n = size(A)
    m != n && return false # only square matrices can satisfy x == one(x)
    return if sizeof(A) < ISONE_CUTOFF
        _isone_triacheck(A, m, approx_test)
    else
        _isone_cachefriendly(A, m, approx_test)
    end
end

@inline function _isone_triacheck(A::StridedMatrix, m::Int, approx_test::AbstractApprox)
    @inbounds for i in 1:m, j in i:m
        if i == j
            isone(A[i, i], approx_test) || return false
        else
            iszero(A[i, j], approx_test) && iszero(A[j, i], approx_test) || return false
        end
    end
    return true
end

# Inner loop over rows to be friendly to the CPU cache
@inline function _isone_cachefriendly(A::StridedMatrix, m::Int, approx)
    @inbounds for i in 1:m
        isone(A[i, i], approx) || return false
        for j in 1:(i - 1)
            iszero(A[j, i], approx) || return false
            iszero(A[i, j], approx) || return false
        end
    end
    return true
end

### ishermitian, issymmetric

# The call `ishermitian(A::AbstractMatrix, B::AbstractMatrix) lowers
# to exactly the same code as that in LinearAlgebra

"""
    ishermitian(A::AbstractMatrix, approx::AbstractApprox)

Return true if A is Hermitian under the notion of approximate equality given by approx.
- Equal: exact Hermitian test.
- Approx: compares A to adjoint(A) via isapprox with norm-based semantics.
- EachApprox: checks symmetry pairwise element-wise against adjoint entries.
Non-square arrays return false.
"""
function ishermitian(A::AbstractMatrix, approx_test::AbstractApprox)
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    @inbounds for i in indsn, j in i:last(indsn)
        isapprox(A[i, j], adjoint(A[j, i]), approx_test) || return false
    end
    return true
end

# This method uses the isapprox interface, which compares using the Frobenius norm.
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

"""
    issymmetric(A::AbstractMatrix, approx::AbstractApprox)

Return true if A is symmetric under approx.
- Equal: exact symmetry.
- Approx: compares A to transpose(A) via isapprox (norm-based for arrays).
- EachApprox: element-wise comparison against transpose(A).
Non-square arrays return false.
"""
issymmetric(A::AbstractMatrix{<:Real}, approx_test::AbstractApprox) =
    ishermitian(A, approx_test)

# Copied from LinearAlgebra. Why does the iteration over i differ slightly from ishermitian ?
function issymmetric(A::AbstractMatrix, approx_test::AbstractApprox)
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    @inbounds for i in first(indsn):last(indsn), j in i:last(indsn)
        isapprox(A[i, j], transpose(A[j, i]), approx_test) || return false
    end
    # for i = first(indsn):last(indsn), j = (i):last(indsn)
    #     if ! isapprox(A[i,j], transpose(A[j,i]), approx_test)
    #         return false
    #     end
    # end
    return true
end

issymmetric(x::Number, approx_test::AbstractApprox) = isapprox(x, x, approx_test)

### isreal

# complex.jl

"""
    isreal(x, approx::AbstractApprox)

Return true if x is real under approx.
- Real numbers: always true.
- Complex numbers: compares real(z) to z under approx (equivalently tests small imaginary part).
- Arrays: true if all entries are real under approx.
- Hermitian/Symmetric wrappers: true if stored data are (approximately) real.
"""
isreal(x::Real, approx_test::AbstractApprox) = true
isreal(z::Complex, approx_test::AbstractApprox) = isapprox(real(z), z, approx_test)
# Old way
#isreal(z::Complex, approx_test::AbstractApprox) = iszero(imag(z), approx_test)

isreal(x::AbstractArray{<:Real}, approx_test::AbstractApprox) = true

isreal(approx_test::AbstractApprox) = x -> isreal(x, approx_test)
isreal(x::AbstractArray, approx_test::AbstractApprox) = all(isreal(approx_test), x)

isreal(A::HermOrSym{<:Real}, approx_test::AbstractApprox) = true
function isreal(A::HermOrSym, approx_test::AbstractApprox)
    n = size(A, 1)
    data = A.data
    isherm = A isa Hermitian
    isupper = A.uplo == 'U'
    @inbounds for j in 1:n
        isreal(data[j, j], approx_test) || return false
    end
    if isupper
        @inbounds for j in 1:n
            iend = j - (isherm ? 1 : 0)
            for i in 1:iend
                isreal(data[i, j], approx_test) || return false
            end
        end
    else
        @inbounds for j in 1:n
            istart = j + (isherm ? 1 : 0)
            for i in istart:n
                isreal(data[i, j], approx_test) || return false
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
"""
    isinteger(x, approx::AbstractApprox)

Return true if x is (approximately) an integer under approx.
- Integer types: always true.
- AbstractFloat: compares x to trunc(x) under approx.
- Complex: requires (approximately) real and integer real part.
- Rational: coerces to float and applies the float rule.
- UpToPhase: for Complex, treats phase-insensitive magnitude for the integer check.
"""
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

"""
    istriu(A::AbstractMatrix, k::Integer, approx::AbstractApprox)

Return true if A is upper triangular with offset k under approx.
Elements with i > j + k must be (approximately) zero.
Larger k relaxes the condition (allows more subdiagonals), negative k tightens it.
Offset arrays are not supported.
"""
function istriu(A::AbstractMatrix, k::Integer, approx::AbstractApprox)
    _require_one_based_indexing(A)
    m, n = size(A)
    @inbounds for j in 1:n
        i0 = j + k + 1
        i0 <= m || continue
        for i in i0:m
            iszero(A[i, j], approx) || return false
        end
    end
    return true
end
istriu(::Number, ::AbstractApprox) = true


"""
    istril(A::AbstractMatrix, k::Integer, approx::AbstractApprox)

Return true if A is lower triangular with offset k under approx.
Elements with j > i + k must be (approximately) zero.
Larger k relaxes the condition (allows more superdiagonals), negative k tightens it.
Offset arrays are not supported.
"""
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

"""
    isbanded(A::AbstractMatrix, kl::Integer, ku::Integer, approx::AbstractApprox)

Return true if A is banded with lower bandwidth kl and upper bandwidth ku under approx.
That is, A[i, j] ≈ 0 when i - j > kl or j - i > ku.
"""
isbanded(A::AbstractMatrix, kl::Integer, ku::Integer, approx::AbstractApprox) =
    istriu(A, kl, approx) && istril(A, ku, approx)

"""
    isdiag(A::AbstractMatrix, approx::AbstractApprox)

Return true if A is diagonal under approx (i.e., banded with kl = ku = 0).
Supports Hermitian and Symmetric wrappers by delegating to triangular views.
"""
isdiag(A::AbstractMatrix, approx::AbstractApprox) = isbanded(A, 0, 0, approx)
isdiag(x::Number, approx::AbstractApprox) = true
isdiag(A::HermOrSym, approx::AbstractApprox) =
    isdiag(
    A.uplo == 'U' ? LinearAlgebra.UpperTriangular(A.data) :
        LinearAlgebra.LowerTriangular(A.data), approx
)
