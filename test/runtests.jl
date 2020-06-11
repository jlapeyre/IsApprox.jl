using IsApprox
using Test

@testset "isone, iszero" begin
    @test IsApprox.isone(1.0)
    @test ! IsApprox.isone(0.0)
    @test IsApprox.iszero(0.0)
    @test ! IsApprox.iszero(1.0)
    @test IsApprox.isone(1.0, Equal())
    @test ! IsApprox.isone(1.0 + 1e-10, Equal())
    @test IsApprox.isone(1.0 + 1e-10, Approx())
    @test IsApprox.isone(1.0 + 1e-7, Approx(atol=1e-6))
    @test ! IsApprox.isone(1.0 + 1e-7, Approx())

    # maximum element is 9.9e-11
    # norm(m) == 1.8892148837710255e-10
    m = [9.914830625828892e-11 4.473996666408231e-11 3.066184129948095e-11;
         9.914231937308404e-12 7.07270051593276e-11 6.487445728681865e-11;
         9.587109787642197e-11 6.520784660158512e-11 1.2903560723094866e-11]

    @test IsApprox.iszero(m, Approx(atol=2e-10)) # in norm
    @test ! IsApprox.iszero(m, Approx(atol=1e-10))
    @test IsApprox.iszero(m, EachApprox(atol=2e-10)) # element-wise
    @test IsApprox.iszero(m, EachApprox(atol=1e-10))
    @test ! IsApprox.iszero(m) # exact
    @test ! IsApprox.iszero(m, Equal()) # exact
end

@testset "isreal, isinteger" begin
    @test IsApprox.isreal(1)
    @test IsApprox.isreal(1.0)
    @test IsApprox.isreal(1.0, Approx())
    @test ! IsApprox.isreal(1.0 + 1e-10im)
    @test ! IsApprox.isreal(1.0 + 1e-10im, Approx())
    @test IsApprox.isreal(1.0 + 1e-10im, Approx(atol=1e-9))
    @test IsApprox.isinteger(1)
    @test IsApprox.isinteger(1, Approx())
end

@testset "isdiag" begin
    m0 = [1 0; 0 1]
    m = [1.0 6.108385298833888e-20; 7.691926633708195e-20 1.0]
    @test IsApprox.isdiag(m0)
    @test ! IsApprox.isdiag(m)
    @test IsApprox.isdiag(m; approx=Approx(atol=1e-10))
end

@testset "isposdef" begin
    n = 3
    m = one(rand(n, n))
    m2 = m + rand(n, n) * 1e-10
    @test IsApprox.ispossemidef(m)
    @test ! IsApprox.ispossemidef(m2)
    @test IsApprox.ispossemidef(m2, Approx(atol=1e-9))
    @test IsApprox.ispossemidef(m2, EachApprox(atol=1e-9))
    @test IsApprox.ispossemidef(0)
    @test IsApprox.ispossemidef(1)
    @test ! IsApprox.ispossemidef(-1)
    @test IsApprox.ispossemidef(-1e-5, Approx(atol=1e-3))

    @test IsApprox.isposdef(1)
    @test ! IsApprox.isposdef(0)
    @test ! IsApprox.isposdef(-1e5)
    @test IsApprox.isposdef(-1e-10, Approx(atol=1e-8))
    @test IsApprox.isposdef(1 + 1e-10im, Approx(atol=1e-8))
    @test ! IsApprox.isposdef(1 + 1e-10im, Equal())

    m3 = [1 0; 0 0]
    @test IsApprox.ispossemidef(m3)
    @test ! IsApprox.isposdef(m3)
    @test IsApprox.isposdef(m3, Approx())
end
