using IsApprox
using IsApprox: isone, iszero, isreal, isinteger, ispossemidef, isposdef, isdiag
using Test

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

@testset "isreal, isinteger" begin
    @test isreal(1)
    @test isreal(1.0)
    @test isreal(1.0, Approx())
    @test ! isreal(1.0 + 1e-10im)
    @test ! isreal(1.0 + 1e-10im, Approx())
    @test isreal(1.0 + 1e-10im, Approx(atol=1e-9))
    @test isinteger(1)
    @test isinteger(1, Approx())
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

@testset "isunitary" begin
    s0 = [1 0; 0 1]
    s1 = [0 1; 1 0]
    s2 = [0 -im; im 0]
    s3 = [1 0; 0 -1]

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
end
