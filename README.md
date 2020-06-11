# IsApprox

## Introduction

`IsApprox` implements an interface for specifying different tests for approximate equality.
It is also fun and hip.

Requirements of `IsApprox` are:
* It should provide a drop in replacement for (and extend) `isapprox` as well as several
application functions, such as `isone` and `issymmetric`.
* Replacements of existing methods must incur no run-time penalty.

## Examples

See this [Jupyter notebook](https://github.com/jlapeyre/IsApprox.jl/blob/master/Notebooks/IsApprox.ipynb)
for examples. See also the [test suite](https://github.com/jlapeyre/IsApprox.jl/blob/master/test/runtests.jl).

## Motivation

For computing eigenvalues, `LinearAlgebra` wants to know if a matrix is exactly
Hermitian. Quantum information packages, on the other hand, might want to know if a matrix
is approximately Hermitian. Furthermore, many functions that check whether a property
(approximately) holds are interdependent. For example `isdiag` calls functions that
eventually call `iszero`. And `isposdef` calls `ishermitian`. Furthermore again, one might
want to check approximate equality in norm; or elementwise for speed. One might want to specify
different tolerances and have them propagate. In practice, packages
(including `QuantumInformation`, `QuantumInfo`, and `Yao`)
tend to reimplement tests in ways that fail to meet these criteria and fail to be composable.
Clearly, a general interface for approximate equality is needed.

## Description

Examples of test types provided are exact equality, elementwise approximate equality, and
approximate equality in norm.  `IsApprox` allows users to specify different definitions of
closeness, via a zero-cost abstraction. That is, specifying the definition of closeness
need not incur a run-time cost.  The code that implements tests for properties such as
symmetry or positivity may then be somewhat decoupled from the specification of
closeness. Furthermore, a simple, small, collection of closeness measures should be
adequate for the vast majority of use cases.

Three subtypes of `AbstractApprox` are included, `Equal`, `Approx` and `EachApprox`.

`IsApprox` implements the interface at least partially for each of: `isone`, `iszero`, `ishermitian`, `issymmetric`,
`isreal`, `isinteger`, `istriu`, `istril`, `isbanded`, `isdiag`, `isposdef`,
`ispossemidef`, `isunitary`, `isinvolution`.

Consider `ishermitian`.

* `ishermitian(A)` or equivalently `ishermitian(A, Equal())` demands exact equality.
This implementation and the function of the same name in `LinearAlgebra` lower to the same code.
That is, the `IsApprox` interface adds no performance penalty.

* `ishermitian(A, Approx(kws...))` has the same semantics as `Base.isapprox`. In this
case, we test that `A` is close to Hermitian in some norm. In this case, a separate code
path is required, namely

```julia
ishermitian(A::AbstractMatrix, approx::Approx) = isapprox(approx, A, adjoint(A))
```

* `ishermitian(A, EachApprox(kws...))`. `EachApprox` specifies element-wise closeness.
If `A` is not close to Hermitian, this test is much faster than `Approx` because
only order `1` elements must be tested. This implementation shares a code path
with that for `Equal`.

## Interfaces

### AbstractApprox

`AbstractApprox`, `Equal`, `Approx`, and `EachApprox` are exported.

### `isapprox`

This extends `Base.isapprox` with methods that take an initial argument of type `AbstractApprox`.
The application functions below take an optional argument of type `AbstractApprox` in the final
position and (may) forward this argument to `isapprox`.

### `isone`, `iszero`, `ishermitian`, etc.

These are not exported, and do not extend the `Base` and `LinearAlgebra` functions of the same names.
They take an optional final argument of type `AbstractApprox`. They are not exported because they
would overwrite existing definitions. However, the `AbstractApprox` interface could be moved into
`Base`.

There are also functions, which *are* exported, that are in neither `Base` nor the standard library, such
as `IsApprox.isunitary`. These follow the parameter ordering and calling conventions
as `IsApprox.isone`, etc.


<!--  LocalWords:  IsApprox isapprox isone issymmetric LinearAlgebra isdiag iszero hoc
 -->
<!--  LocalWords:  isposdef ishermitian elementwise reimplement ishermitan QuantumInfo
 -->
<!--  LocalWords:  QuantumInformation Yao positivity subtypes AbstractApprox EachApprox
 -->
<!--  LocalWords:  kws julia adjoint isunitary
 -->
