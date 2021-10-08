### The functions in this file extend functions of the same names in `Base` and `LinearAlgebra`.
### It may be necessary to instead replace the `Base` and `LinearAlgebra` functions
### if `IsApprox` is widely adopted.

### isone, iszero

isone(x, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, x, one(x))
iszero(x, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, x, zero(x))
iszero(x::AbstractArray, approx_test::AbstractApprox=Equal()) =
    all(y -> iszero(y, approx_test), x)
isone(x::BigInt, ::Equal) = Base.isone(x)
iszero(x::BigInt, ::Equal) = Base.iszero(x)

isone(A::StridedMatrix, approx_test::Approx) = isapprox(approx_test, A, one(A))
iszero(A::StridedMatrix, approx_test::Approx) = isapprox(approx_test, A, zero(A))

# dense.jl
const ISONE_CUTOFF = 2^21 # 2M

function isone(A::StridedMatrix, approx_test::AbstractApprox=Equal())
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
function ishermitian(A::AbstractMatrix, approx_test::AbstractApprox=Equal())
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    for i = indsn, j = i:last(indsn)
        if ! isapprox(approx_test, A[i,j], adjoint(A[j,i]))
            return false
        end
    end
    return true
end

# This uses the isapprox interface, which compares using a norm
ishermitian(A::AbstractMatrix, approx_test::Approx) = isapprox(approx_test, A, adjoint(A))

ishermitian(x::Number, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, x, conj(x))

issymmetric(A::AbstractMatrix{<:Real}, approx_test::AbstractApprox=Equal()) =
    ishermitian(A, approx_test)

# Copied from LinearAlgebra. Why does the iteration over i differ slightly from ishermitian ?
function issymmetric(A::AbstractMatrix, approx_test::AbstractApprox=Equal())
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    for i = first(indsn):last(indsn), j = (i):last(indsn)
        if ! isapprox(approx_test, A[i,j], transpose(A[j,i]))
            return false
        end
    end
    return true
end

issymmetric(x::Number, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, x, x)

### isreal

# complex.jl
isreal(x::Real, approx_test::AbstractApprox=Equal()) = true
isreal(z::Complex, approx_test::AbstractApprox=Equal()) = iszero(imag(z), approx_test)

### isinteger

# This must be changed if this is integrated into Base
isinteger(x) = isinteger(x, Equal())
# We use union or explicit types to avoid method ambiguity. Is there a way around this ?
isinteger(x::BigFloat, ::Equal) = Base.isinteger(x)
isinteger(x::Rational, ::Equal) = Base.isinteger(x)

# number.jl
isinteger(x::Integer, ::AbstractApprox=Equal()) = true

# floatfuncs.jl
## The original is x - trunc(x) == 0. So this implementation might differ
## from the stock version of isinteger for some user's subtype of AbstractFloat.
## We choose to this implementation because the default relative tolerance is
## reasonable. That is, `isapprox(approx_test, x - trunc(x), 0)` requires
## specifying `atol`.
isinteger(x::AbstractFloat, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, x, trunc(x))

# complex.jl
isinteger(z::Complex, approx_test::AbstractApprox=Equal()) =
    isreal(z, approx_test) & isinteger(real(z), approx_test)

isinteger(x::Rational, approx_test::AbstractApprox) = isinteger(float(x), approx_test)

### istriu, istril, isbanded, isdiag

function istriu(A::AbstractMatrix, k::Integer = 0; approx::AbstractApprox=Equal())
    Base.require_one_based_indexing(A)
    m, n = size(A)
    for j in 1:min(n, m + k - 1)
        for i in max(1, j - k + 1):m
            iszero(A[i, j], approx) || return false
        end
    end
    return true
end
istriu(x::Number, ::AbstractApprox) = true

function istril(A::AbstractMatrix, k::Integer = 0; approx::AbstractApprox=Equal())
    Base.require_one_based_indexing(A)
    m, n = size(A)
    for j in max(1, k + 2):n
        for i in 1:min(j - k - 1, m)
            iszero(A[i, j], approx) || return false
        end
    end
    return true
end
istril(x::Number, ::AbstractApprox) = true

isbanded(A::AbstractMatrix, kl::Integer, ku::Integer; approx::AbstractApprox=Equal()) =
    istriu(A, kl; approx=approx) && istril(A, ku; approx=approx)

isdiag(A::AbstractMatrix; approx::AbstractApprox=Equal()) = isbanded(A, 0, 0; approx=approx)
isdiag(x::Number; approx::AbstractApprox=Equal()) = true
