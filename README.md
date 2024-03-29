[![Build Status](https://github.com/jlapeyre/IsApprox.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jlapeyre/IsApprox.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jlapeyre/IsApprox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlapeyre/IsApprox.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET QA](https://img.shields.io/badge/JET.jl-%E2%9C%88%EF%B8%8F-%23aa4444)](https://github.com/aviatesk/JET.jl)

# IsApprox

## Introduction

`IsApprox` implements an interface for applying different definitions of "approximate" in tests for approximate (or exact) equality.
It is also fun and hip.

Design requirements of `IsApprox` are:

* It should provide a drop-in replacement for (and extend) `isapprox` as well as
several application functions, such as `isone` and `issymmetric`.  In
particular, many functions that currently check for a property exactly (to
machine precision) will instead use `IsApprox` to implement both exact and
approximate comparison.

* Replacements of existing methods (eg. `isone(::Float64)`) must incur no run-time penalty.
In practice, this means specifying the notion of "approximate" via types, eg `Equal` and `Approx`
so that the compiler inlines the comparison code.

## Examples

See this [Jupyter notebook](https://github.com/jlapeyre/IsApprox.jl/blob/master/Notebooks/IsApprox.ipynb)
for examples. See also the [test suite](https://github.com/jlapeyre/IsApprox.jl/blob/master/test/runtests.jl).

## Motivation

For some applications, `LinearAlgebra` wants to know if a matrix is exactly
Hermitian. Quantum information packages, on the other hand, might want to know
if a matrix is approximately (or exactly) Hermitian. Furthermore, many functions that check
whether a property (approximately) holds are interdependent. For example
`isdiag` calls functions that eventually call `iszero`. And `isposdef` calls
`ishermitian`. Furthermore again, one might want to check approximate equality
in norm; or elementwise. One might want to specify a tolerance and have it
propagate. In practice, packages
tend to reimplement tests in ways that do not satisfy all these criteria,
and fail to be composable.
Such packages include [`QuantumInformation`](https://github.com/iitis/QuantumInformation.jl)(
[code example](https://github.com/iitis/QuantumInformation.jl/blob/b47400ebb09d10cc1eba5f7bf06badeb6cfe5429/src/utils.jl#L93-L113))
,
[`QuantumInfo`](https://github.com/BBN-Q/QuantumInfo.jl)(
[code example](https://github.com/BBN-Q/QuantumInfo.jl/blob/cbafdc7f295e7d56e41116c5ef4eca9500d45909/src/basics.jl#L248-L255)),
and
[`Yao`](https://github.com/QuantumBFS/Yao.jl)(
[code example](https://github.com/QuantumBFS/YaoBase.jl/blob/master/src/inspect.jl)).
Clearly, a general interface for approximate
equality is needed.

## Description

`IsApprox` allows users to specify different definitions of
closeness, via a zero-cost abstraction. That is, specifying the definition of closeness
need not incur a run-time cost.  The code that implements tests for properties such as
symmetry or positivity may then be somewhat decoupled from the specification of
closeness. Furthermore, a simple, small, collection of closeness measures should be
adequate for the vast majority of use cases.

Four subtypes of `AbstractApprox` are included, `Equal`, `Approx`, `EachApprox`, and `UpToPhase`.


`IsApprox` implements the interface at least partially for each of: `isone`, `iszero`, `ishermitian`, `issymmetric`,
`isreal`, `isinteger`, `istriu`, `istril`, `isbanded`, `isdiag`, `isposdef`,
`ispossemidef`, `isunitary`, `isinvolution`, `isnormalized`, `isprobdist`.

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

## API

### AbstractApprox

`AbstractApprox`, `Equal`, `Approx`, `UpToPhase`, and `EachApprox` are exported.

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


## Style

This package will probably try to follow the [Blue Style Guide](https://github.com/invenia/BlueStyle).
An important rule is broken immediately: predicates are written `isprop` rather than `is_prop`.
And `ispropmod1mod2` rather than `is_prop_mod1_mod2`. The main reason is that some of these
functions exist by the same name in `Base`. And some are very closely related.

<!--  LocalWords:  IsApprox isapprox isone issymmetric LinearAlgebra isdiag iszero hoc
 -->
<!--  LocalWords:  isposdef ishermitian elementwise reimplement ishermitan QuantumInfo
 -->
<!--  LocalWords:  QuantumInformation Yao positivity subtypes AbstractApprox EachApprox
 -->
<!--  LocalWords:  kws julia adjoint isunitary
 -->
