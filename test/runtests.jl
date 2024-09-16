using IsApprox
using IsApprox: IsApprox, isone, iszero, isreal, isinteger, ispossemidef, isposdef, isdiag
using Dictionaries: Dictionary
using Test
import LinearAlgebra

if VERSION >= v"1.7" && VERSION <= v"1.11"
    @testset "JET" begin
        include("jet_test.jl")
    end
end

include("aqua_test.jl")

@testset "isnormalized, isprobdist" begin
    v1 = collect(1:10) ./ sum(1:10)
    g1 = (x for x in v1)
    d1 = Dict(i => v1[i] for i in eachindex(v1))
    for c in (v1, g1, d1)
        @test isprobdist(c)
        @test isnormalized(c)
    end
    _v2 = collect(1:10)
    _v2[5] = -5
    v2 = _v2 ./ sum(_v2)
    g2 = (x for x in v2)
    d2 = Dict(i => v2[i] for i in eachindex(v2))
    for c in (v2, g2, d2)
        @test ! isprobdist(c)
        @test isnormalized(c)
    end
    v3 = copy(v2)
    v3[end] *= 2
    @test !isprobdist(v3)

    v4 = [x + 1e-10 * (rand() - 0.5) for x in v1]
    @test ! isprobdist(v4)
    @test isprobdist(v4, Approx())
    @test ! isprobdist(v4, Approx(atol=1e-16))
end


@testset "isone, iszero" begin
    @test isone(1.0)
    @test ! isone(0.0)
    @test iszero(0.0)
    @test ! iszero(1.0)
    @test isone(1.0, Equal())
    @test ! isone(1.0 + 1e-10, Equal())
    @test isone(1.0 + 1e-10, Approx())
    @test isone(1.0 + 1e-7, Approx(atol=1e-6))
    @test ! isone(1.0 + 1e-7, Approx())

    # maximum element is 9.9e-11
    # norm(m) == 1.8892148837710255e-10
    m = [9.914830625828892e-11 4.473996666408231e-11 3.066184129948095e-11;
         9.914231937308404e-12 7.07270051593276e-11 6.487445728681865e-11;
         9.587109787642197e-11 6.520784660158512e-11 1.2903560723094866e-11]

    @test iszero(m, Approx(atol=2e-10)) # in norm
    @test ! iszero(m, Approx(atol=1e-10))
    @test iszero(m, EachApprox(atol=2e-10)) # element-wise
    @test iszero(m, EachApprox(atol=1e-10))
    @test ! iszero(m) # exact
    @test ! iszero(m, Equal()) # exact
end

@testset "isreal" begin
    isreal = IsApprox.isreal

    @test isreal(1)
    @test isreal(1.0)
    @test isreal(1.0, Approx())
    @test ! isreal(1.0 + 1e-10im)
    @test ! isreal(1.0 + 1e-7im, Approx())
    @test isreal(1.0 + 1e-10im, Approx(atol=1e-9))

    m = real.(rand(ComplexF64, 3, 3))
    @test isreal(m)
    m_noisy = m + 1e-14 * randn(ComplexF64, 3, 3)
    @test !isreal(m_noisy)

    for approx in (Approx, EachApprox)
        @test isreal(m_noisy, approx())
        @test isreal(m_noisy, approx(atol=1e-8))
        @test !isreal(m_noisy, approx(atol=1e-16))
        @test isreal(m_noisy, approx())
        @test isreal(m_noisy, approx(atol=1e-8))
        @test !isreal(m_noisy, approx(atol=1e-16))
    end
end

@testset "isinteger" begin
    @test isinteger(1)
    @test isinteger(1, Approx())
    @test isinteger(1 + 1e-8, Approx())
    @test ! isinteger(1 + 1e-5, Approx())
    @test isinteger(1000 + 1e-5, Approx())
    @test isinteger(1 + im * 1e-8, Approx())
    @test ! isinteger(1 + im * 1e-5, Approx())
    @test isinteger(exp(im*2), UpToPhase())
    @test ! isinteger(1.1 * exp(im*2), UpToPhase())
end

@testset "isdiag" begin
    m0 = [1 0; 0 1]
    m = [1.0 6.108385298833888e-20; 7.691926633708195e-20 1.0]
    @test isdiag(m0)
    @test ! isdiag(m)
    @test isdiag(m; approx=Approx(atol=1e-10))
end

@testset "isposdef" begin
    n = 3
    m = one(rand(n, n))
    m2 = m + rand(n, n) * 1e-10
    @test ispossemidef(m)
    @test ! ispossemidef(m2)
    @test ispossemidef(m2, Approx(atol=1e-9))
    @test ispossemidef(m2, EachApprox(atol=1e-9))
    @test ispossemidef(0)
    @test ispossemidef(1)
    @test ! ispossemidef(-1)
    @test isposdef(1 + 0im)
    @test ispossemidef(-1e-5, Approx(atol=1e-3))

    @test isposdef(1)
    @test ! isposdef(0)
    @test ! isposdef(-1e5)
    @test isposdef(-1e-10, Approx(atol=1e-8))
    @test isposdef(1 + 1e-10im, Approx(atol=1e-8))
    @test ! isposdef(1 + 1e-10im, Equal())

    m3 = [1 0; 0 0]
    @test ispossemidef(m3)
    @test ! isposdef(m3)
    @test isposdef(m3, Approx())
end

@testset "isunitary isinvolution" begin
    s0 = [1 0; 0 1]
    s1 = [0 1; 1 0]
    s2 = [0 -im; im 0]
    s3 = [1 0; 0 -1]

    @test isunitary(s3)
    @test isunitary(s2)
    @test isinvolution(s3)
    @test isinvolution(s2)

    # merr = rand(0:100, 8,8) * 1e-12
    errmat = [5.3e-11 9.4e-11 3.5e-11 6.6e-11 8.6e-11 4.9e-11 5.2e-11 1.6e-11;
              2.8e-11 3.9e-11 4.0e-11 6.0e-12 3.1e-11 2.9e-11 9.7e-11 7.1e-11;
              8.0e-11 2.5e-11 8.5e-11 7.2e-11 5.7e-11 4.8e-11 5.8e-11 4.1e-11;
              5.5e-11 9.8e-11 4.7e-11 8.7e-11 3.6e-11 9.0e-11 4.3e-11 8.0e-11;
              7.0e-11 1.5e-11 2.1e-11 9.8e-11 9.7e-11 2.1e-11 7.6e-11 4.0e-11;
              4.8e-11 3.9e-11 9.8e-11 0.0 4.0e-11 4.0e-12 4.0e-11 5.3e-11;
              6.8e-11 8.3e-11 4.0e-12 1.5e-11 6.7e-11 3.4e-11 1.0e-11 9.0e-12;
              7.0e-11 5.2e-11 2.2e-11 7.2e-11 9.0e-12 5.4e-11 4.9e-11 5.3e-11]

    m = kron(s1, s2, s3)
    merr = m + errmat

    @test isunitary(m)
    @test isunitary(m, Equal())
    @test ! isunitary(merr)
    @test ! isunitary(merr, Equal())
    @test isunitary(merr, Approx())
    @test ! isunitary(merr, EachApprox())
    @test ! isunitary(merr, EachApprox(atol=1e-10))
    @test isunitary(merr, EachApprox(atol=5e-10))
    @test ! isunitary(merr, Approx(atol=1e-10))

    m1 = [1 0; 1 -1]
    _m1err = (rand(2, 2) .- 0.5) * 1e-10
    m1err = m1 + _m1err

    @test isinvolution(m1)
    @test ! isunitary(m1)
    @test isinvolution(m1, Equal())
    @test ! isinvolution(m1err)
    @test isinvolution(m1err, EachApprox(atol=1e-9))
    @test isinvolution(m1err, Approx(atol=1e-9))
    @test ! isunitary(m1err, Approx(atol=1e-9))
    @test ! isinvolution(m1err, Approx(atol=1e-12))

    @test isunitary(1)
    @test ! isunitary(2)
    @test isunitary((1e-10 + exp(im * 2.5)) * LinearAlgebra.I, Approx())

end

@testset "UpToPhase" begin
    @test isapprox(UpToPhase(), 2, 2)
    @test isapprox(UpToPhase(), 2.1, 2.1 + 1e-10)
    @test ! isapprox(UpToPhase(atol=1e-12), 2.1, 2.1 + 1e-10)
    @test isapprox(UpToPhase(), 2.1, 2.1 * cis(2*pi*3.1))
    m = rand(5, 5)
    @test isapprox(UpToPhase(), m, m .* cis(2*pi*1.3))
    @test isapprox(UpToPhase(), m, m .* (1+1e-10)*cis(2*pi*1.3))
    @test ! isapprox(UpToPhase(), m, m .* (1+1e-7)*cis(2*pi*1.3))
end

@testset "Dictionary" begin
    # bug at commit ebfe206de7d
    n = 5
    _v = rand(n)
    v = _v ./ sum(_v)
    d = Dictionary(1:n, v)
    @test IsApprox._all_possemidef(d, Approx())
    @test isnormalized(d, Approx())
    @test isprobdist(d, Approx())
end

@testset "Bad construction" begin
    @test_throws MethodError Approx(3)
    @test_throws MethodError Approx(1e-10)
    @test_throws MethodError EachApprox(1e-10)
    @test_throws MethodError EachApprox("dog")
    @test_throws MethodError UpToPhase(1e-10)
end

# These methods only exist to avoid method ambiguity
@testset "BigInt BigFloat Rational" begin
    isone = IsApprox.isone
    iszero = IsApprox.iszero
    isinteger = IsApprox.isinteger
    @test isone(big(1))
    @test !isone(big(0))
    @test !iszero(big(1))
    @test iszero(big(0))
    @test isone(big(1), Equal())
    @test !isone(big(0), Equal())
    @test !iszero(big(1), Equal())
    @test iszero(big(0), Equal())
    @test isinteger(rationalize(42.0))
    @test ! isinteger(1//2)
end

@testset "Hermitian Symmetric" begin
    ishermitian = IsApprox.ishermitian
    isreal = IsApprox.isreal
    issymmetric = IsApprox.issymmetric
    isdiag = IsApprox.isdiag

    m = rand(ComplexF64, 4, 4);
    m_herm = m + m'
    @test !ishermitian(m)
    @test ishermitian(m_herm)
    noise = 1e-14 * randn(ComplexF64, 4, 4);
    m_herm_noisy = m_herm + noise
    @test ! ishermitian(m_herm_noisy)
    @test ! ishermitian(m_herm_noisy, Equal())
    @test ishermitian(m_herm_noisy, Approx())
    @test ishermitian(m_herm_noisy, Approx(atol=1e-8))
    @test ! ishermitian(m_herm_noisy, Approx(atol=1e-16))

    m_Herm = LinearAlgebra.Hermitian(m)
    @test ishermitian(m_Herm)
    @test ishermitian(m_Herm, Approx())
    @test !isreal(m_herm)
    @test !isreal(m_herm, Approx())
    @test !isreal(m_Herm)
    @test !isreal(m_Herm, Approx())
    @test !isdiag(m_herm)
    @test !isdiag(m_Herm)
    @test !isdiag(m_Herm; approx=Approx())

    m = rand(Float64, 4, 4);
    m_symm = m + m'
    @test !ishermitian(m)
    @test ishermitian(m_symm)
    noise = 1e-14 * randn(Float64, 4, 4);
    m_symm_noisy = m_symm + noise
    @test ! ishermitian(m_symm_noisy)
    @test ! ishermitian(m_symm_noisy, Equal())
    @test ishermitian(m_symm_noisy, Approx())
    @test ishermitian(m_symm_noisy, Approx(atol=1e-8))
    @test ! ishermitian(m_symm_noisy, Approx(atol=1e-16))

    m_Symm = LinearAlgebra.Symmetric(m)
    @test ishermitian(m_Symm)
    @test ishermitian(m_Symm, Approx())
    @test isreal(m_symm)
    @test isreal(m_symm, Approx())
    @test isreal(m_Symm)
    @test isreal(m_Symm, Approx())
    @test !isdiag(m_symm)
    @test !isdiag(m_Symm)
    @test !isdiag(m_Symm, approx=Approx())
end

@testset "istriu, istril" begin
    istriu = IsApprox.istriu
    istril = IsApprox.istril
    for istri in (istril, istril)
        @test istri(1)
        @test istri(1, Approx())
        @test istri(1, EachApprox())
    end
end
