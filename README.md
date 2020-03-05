# IsApprox

`IsApprox` implements an interface for specifying different tests for approximate equality.
Examples are exact equality, elementwise approximate equality, and approximate equality in norm.
Currently, packages often implement their own version of functions such as `ishermitan` because
the required notion of close depends on the use case. `IsApprox` is an attempt to allow
users to instead specify different notions of closeness. The code that implements
tests for properties such as symmetry or positivity, may then be somewhat decoupled from the
specification of closeness. Futhermore, a simple, small, collection of closeness measures
should be adequate for the vast majority of use cases.

Three subtypes of `AbstractApprox` are included, `Equal`, `Approx` and `EachApprox`.

An example application, `IsApprox.ishermitian` is included.

* `ishermitian(A)` or equivalently `ishermitian(A, Equal())` demands exact equality.
This implemenation and the function of the same name in `LinearAlgebra` lower to the same code.
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
