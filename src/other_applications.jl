import LinearAlgebra

## We use ispossemidef for numbers as well as matrices.
## This follows `ishermitian`, etc.

import LinearAlgebra: isposdef

ispossemidef(x) = ispossemidef(x, Equal())

"""
    ispossemidef(x::Number, approx_test::AbstractApprox)

Return `true` if `x` is (approximately) non-negative.
"""
function ispossemidef(x::Real, approx_test::AbstractApprox)
    return x > zero(x) || iszero(x, approx_test)
end

## Used for both pos def and pos semidef
function _isposdef(z::Complex, approx_test::AbstractApprox, posdeffunc)
    return posdeffunc(real(z), approx_test) && iszero(imag(z), approx_test)
end

function ispossemidef(z::Complex, approx_test::AbstractApprox)
    return _isposdef(z, approx_test, ispossemidef)
end

function isposdef(z::Complex, approx_test::AbstractApprox)
    return _isposdef(z, approx_test, isposdef)
end

"""
    isposdef(x::Number, approx_test::AbstractApprox)

Return `true` if `x` is approximately greater than zero.
"""
isposdef(x::Real, ::Equal) = x > zero(x)
#isposdef(x::Number) = isposdef(x, Equal())
## For methods other than `Equal`, x can be zero or negative
isposdef(x::Number, approx_test::AbstractApprox) = ispossemidef(x, approx_test)

## For both positive definite and positive semidefinite
function _isposdef(m::AbstractMatrix, approx_test::AbstractApprox, posdeffunc)
    ! ishermitian(m, approx_test) && return false
    evs = LinearAlgebra.eigvals(LinearAlgebra.Hermitian(m))
    return all(x -> posdeffunc(x, approx_test), evs)
end

"""
    ispossemidef(m::AbstractMatrix, approx_test::AbstractApprox)

Return `true` if `m` is positive semidefinite.
"""
function ispossemidef(m::AbstractMatrix, approx_test::AbstractApprox)
    return _isposdef(m, approx_test, ispossemidef)
end

"""
    isposdef(m::AbstractMatrix, approx_test::AbstractApprox)

Return `true` if `m` is positive definite.
"""
function isposdef(m::AbstractMatrix, approx_test::AbstractApprox)
    return _isposdef(m, approx_test, isposdef)
end

# TODO: I think I no longer need this, after chaning arg order,
# and importing all predicate functions
#isposdef(A::AbstractMatrix) = isposdef(A, Equal())
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
function _isunitary(m::AbstractMatrix, approx_test::AbstractApprox, dotf::F1, transposef::F2) where {F1, F2}
    _one = one(eltype(m))
    rowinds = axes(m)[2]
    for i in rowinds
        isapprox(dotf(view(m, :, i), view(transposef(m), :, i)), _one, approx_test) || return false
        for j in i+1:last(rowinds)
            isapprox(dotf(view(m, :, i), view(transposef(m), :, j)) + _one, _one, approx_test) || return false
        end
    end
    return true
end

isunitary(x) = isunitary(x, Equal())

## Use vector norm.
## Slower, but more generally useful.
function isunitary(m::AbstractMatrix, approx_test::Approx)
    return  isapprox(m' * m, LinearAlgebra.I, approx_test)
end

_identity(x) = x
"""
    isunitary(m::AbstractMatrix, approx_test::AbstractApprox)

Return `true` if `m` is unitary. If `m` is real, this tests orthogonality.
"""
isunitary(m::AbstractMatrix, approx_test::AbstractApprox) =
    _isunitary(m, approx_test, LinearAlgebra.dot, _identity)

# abs2 is much faster, but we would need to use sqrt to adjust the tolerance, thus losing any advantage.
isunitary(x::Number, approx_test::AbstractApprox) = isone(abs(x), approx_test)
isunitary(J::LinearAlgebra.UniformScaling, approx_test::AbstractApprox) = isunitary(J.λ, approx_test)

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

isinvolution(x) = isinvolution(x, Equal())

"""
    isinvolution(m::AbstractMatrix, approx_test::AbstractApprox)

Return `true` if `m * m == I`
"""
isinvolution(m::AbstractMatrix, approx_test::AbstractApprox) = _isunitary(m, approx_test, _dot, transpose)

function isinvolution(m::AbstractMatrix, approx_test::Approx)
    return  isapprox(m * m, LinearAlgebra.I, approx_test)
end

isidempotent(x) = isidempotent(x, Equal())

function isidempotent(m::AbstractMatrix, approx_test::AbstractApprox)
    return isapprox(m * m, m, approx_test)
end

isnormal(x) = isnormal(x, Equal())

function isnormal(m::AbstractMatrix, approx_test::AbstractApprox)
    return isapprox(m * m', m' * m, approx_test)
end

"""
    commutes(X, Y, approx_test::AbstractApprox)

Return `true` if `X` and `Y` commute.
"""
commutes(X, Y, approx_test::AbstractApprox) = isapprox(X * Y, Y * X, approx_test)

commutes(X, Y) = commutes(X, Y, Equal())

"""
    anticommutes(X, Y, approx_test::AbstractApprox)

Return `true` if `X` and `Y` anticommute.
"""
anticommutes(X, Y, approx_test::AbstractApprox) = isapprox(X * Y, -(Y * X), approx_test)

anticommutes(X, Y) = anticommutes(X, Y, Equal())

"""
    isnormalized(itr, approx_test::AbstractApprox)

Return `true` if `itr` is normalized, that is, if the items sum to one.
If `itr` is a dictionary, the items are the values.
"""
isnormalized(itr, approx_test::AbstractApprox) = isnormalized(Base.IteratorEltype(itr), itr, approx_test)
isnormalized(::Base.EltypeUnknown, itr, approx_test::AbstractApprox) = isapprox(sum(itr), 1, approx_test)
isnormalized(::Base.HasEltype, itr, approx_test::AbstractApprox) = isapprox(sum(itr), one(eltype(itr)), approx_test)
isnormalized(d::_AbstractDict{<:Any, V}, approx_test::AbstractApprox) where V = isapprox(sum(values(d)), one(V), approx_test)
# isnormalized(x::Base.HasEltype, approx_test::AbstractApprox) = throw(MethodError(isnormalized, (x, approx_test)))
# isnormalized(x::Base.EltypeUnknown, approx_test::AbstractApprox) = throw(MethodError(isnormalized, (x, approx_test)))
isnormalized(x) = isnormalized(x, Equal())

# A Vector has no algebraic interpretation wrt isposdef, ispossemidef. So we don't define
# a method. To help isprobdist, we do the following.
_all_possemidef(itr, approx_test) = all(x -> ispossemidef(x, approx_test), itr)
_all_possemidef(d::AbstractDict, approx_test) = _all_possemidef(values(d), approx_test)

"""
    isprobdist(itr, approx_test=Equal())

Return `true` if the items in `itr` form a probability distribution, that is
sum to one and are non-negative.
If `itr` is a dictionary, `values(itr)` is tested.
"""
isprobdist(itr, approx_test=Equal()) = isnormalized(itr, approx_test) && _all_possemidef(itr, approx_test)

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
# [done] isnormalized
# isreflexive  <-- isinvolution is a better name, See above

# From other libraries
# isidempotent
# isnormal
