import LinearAlgebra

## We use ispossemidef for numbers as well as matrices.
## This follows `ishermitian`, etc.

"""
    ispossemidef(x::Number, approx_test::AbstractApprox=Equal())

Return `true` if `x` is (approximately) non-negative.
"""
function ispossemidef(x::Real, approx_test::AbstractApprox=Equal())
    return x > zero(x) || iszero(x, approx_test)
end

## Used for both pos def and pos semidef
function _isposdef(z::Complex, approx_test::AbstractApprox, posdeffunc)
    return posdeffunc(real(z), approx_test) && iszero(imag(z), approx_test)
end

function ispossemidef(z::Complex, approx_test::AbstractApprox=Equal())
    return _isposdef(z, approx_test, ispossemidef)
end

function isposdef(z::Complex, approx_test::AbstractApprox=Equal())
    return _isposdef(z, approx_test, isposdef)
end

"""
    isposdef(x::Number, approx_test::AbstractApprox=Equal())

Return `true` if `x` is approximately greater than zero.
"""
isposdef(x::Number) = isposdef(x, Equal())
isposdef(x::Real, ::Equal) = x > zero(x)
## For methods other than `Equal`, x can be zero or negative
isposdef(x::Number, approx_test::AbstractApprox) = ispossemidef(x, approx_test)

## For both positive definite and positive semidefinite
function _isposdef(m::AbstractMatrix, approx_test::AbstractApprox, posdeffunc)
    ! ishermitian(m, approx_test) && return false
    evs = LinearAlgebra.eigvals(LinearAlgebra.Hermitian(m))
    return all(x -> posdeffunc(x, approx_test), evs)
end

"""
    ispossemidef(m::AbstractMatrix, approx_test::AbstractApprox=Equal())

Return `true` if `m` is positive semidefinite.
"""
function ispossemidef(m::AbstractMatrix, approx_test::AbstractApprox=Equal())
    return _isposdef(m, approx_test, ispossemidef)
end

"""
    isposdef(m::AbstractMatrix, approx_test::AbstractApprox=Equal())

Return `true` if `m` is positive definite.
"""
function isposdef(m::AbstractMatrix, approx_test::AbstractApprox)
    return _isposdef(m, approx_test, isposdef)
end

isposdef(A::AbstractMatrix) = isposdef(A, Equal())
# copied from dense.jl
isposdef(A::AbstractMatrix, ::Equal) =
    ishermitian(A, Equal()) && LinearAlgebra.isposdef(LinearAlgebra.cholesky(LinearAlgebra.Hermitian(A); check = false))

## Compared two methods:
## a) Allocate, ie m' * m. b) iterate over columns
## Found:
## 1. iterating order, ie rows vs cols is correct
## 2. Avoiding allocation improves efficiency even for 2x2 matrices, but...
## 3. Doing the allocation, eg m' * m can be slightly faster. Eg for 100x100 dense identity matrix.
## 4. For rand(100, 100), iterating over columns is 1000 times faster. Fails on first column.
## `approx_test` is `Equal` or `EachApprox`.
function _isunitary(m::AbstractMatrix, approx_test::AbstractApprox, dotf, transposef)
    rowinds = axes(m)[2]
    for i in rowinds
        isapprox(approx_test, dotf(view(m, :, i), view(transposef(m), :, i)), 1) || return false
        for j in i+1:last(rowinds)
            isapprox(approx_test, dotf(view(m, :, i), view(transposef(m), :, j)), 0) || return false
        end
    end
    return true
end

## This is slower even for small matrices.
## Then why am I using it ? (May 2021)
function isunitary(m::AbstractMatrix, approx_test::Approx)
    return  isapprox(approx_test, m' * m, LinearAlgebra.I)
end

_identity(x) = x
"""
    isunitary(m::AbstractMatrix, approx_test::AbstractApprox=Equal())

Return `true` if `m` is unitary. If `m` is real, this tests orthogonality.
"""
isunitary(m::AbstractMatrix, approx_test::AbstractApprox=Equal()) =
    _isunitary(m, approx_test, LinearAlgebra.dot, _identity)

isunitary(x::Number, approx_test::AbstractApprox=Equal()) = isone(abs(x), approx_test)
isunitary(J::LinearAlgebra.UniformScaling, approx_test::AbstractApprox=Equal()) = isunitary(J.λ, approx_test)

"""
    _dotu(x::AbstractVector, y::AbstractVector)

Dot product with no complex conjugation of `x`. This dispatches to `LinearAlgebra.dot`
if `x` and `y` are `Real`.
"""
_dot(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = LinearAlgebra.dot(x, y)

_dot(x::AbstractVector{Complex{T}}, y::AbstractVector{Complex{T}}) where {T <: Union{Float64, Float32}} =
    LinearAlgebra.BLAS.dotu(x, y)

## There has got to be a better way than all of these _dot methods
function _dot(x::AbstractVector{<:Complex}, y::AbstractVector{<:Complex})
    return sum(*(z...) for z in zip(x, y))
end

"""
    isinvolution(m::AbstractMatrix, approx_test::AbstractApprox=Equal())

Return `true` if `m * m == I`
"""
isinvolution(m::AbstractMatrix, approx_test::AbstractApprox=Equal()) = _isunitary(m, approx_test, _dot, transpose)

function isinvolution(m::AbstractMatrix, approx_test::Approx)
    return  isapprox(approx_test, m * m, LinearAlgebra.I)
end

function isidempotent(m::AbstractMatrix, approx_test::AbstractApprox=Equal())
    return isapprox(approx_test, m * m, m)
end

function isnormal(m::AbstractMatrix, approx_test::AbstractApprox=Equal())
    return isapprox(approx_test, m * m', m' * m)
end

"""
    commutes(X, Y, approx_test::AbstractApprox=Equal())

Return `true` if `X` and `Y` commute.
"""
commutes(X, Y, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, X * Y, Y * X)

"""
    anticommutes(X, Y, approx_test::AbstractApprox=Equal())

Return `true` if `X` and `Y` anticommute.
"""
anticommutes(X, Y, approx_test::AbstractApprox=Equal()) = isapprox(approx_test, X * Y, -(Y * X))

# Some of the following from QuantumInfo.jl
# might be implemented here:
# ischannel
# iscp
# isliouvillian
# istp
# isunital

# Yao.jl
#
# iscommute
# isnormalized
# isreflexive  <-- isinvolution is a better name, See above

# From other libraries
# isidempotent
# isnormal
