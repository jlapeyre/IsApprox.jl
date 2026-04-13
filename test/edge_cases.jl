using IsApprox: isnormal, isdiag, isunitary, ishermitian, issymmetric, isposdef

@testset "edge cases" begin
    @testset "isinvolution edge-case (upper-only error)" begin
        using LinearAlgebra: I
        N = zeros(3, 3)
        N[1, 2] = 1.0
        N[2, 3] = 1.0
        m = Matrix{Float64}(I, 3, 3) + N
        @test !isinvolution(m, Equal())
    end

    @testset "triangular and banded k offsets" begin
        A = [1 0 0; 1e-12 1 0; 0 1e-12 1];
        @test IsApprox.istriu(A, 0, EachApprox(atol=1e-9));
        @test !IsApprox.istriu(A, -1, Equal());
        @test IsApprox.isbanded(A, 1, 0, EachApprox(atol=1e-9));
    end

    @testset "UpToPhase zero-handling (arrays)" begin
        A = zeros(ComplexF64, 3);
        B = zeros(ComplexF64, 3);
        @test isapprox(A, B, UpToPhase());
        @test !isapprox(A, [0, 0, 1], UpToPhase());
    end

    @testset "UpToPhase leading zeros" begin
        A = [0.0, 1.0, 2.0];
        B = cis(0.7) .* A;
        B[1] = 0.0;
        @test isapprox(A, B, UpToPhase());
    end

    @testset "isnormal on scalars (current behavior)" begin
        @test_throws MethodError isnormal(1.0);
        @test_throws MethodError isnormal(0.0);
        @test isnormal([1 0; 0 1]);
    end

    @testset "isdiag for Hermitian uplo paths" begin
        M = [1 0 0; 0 2 1e-12im; 0 -1e-12im 3];
        @test isdiag(Hermitian(M, :U), EachApprox(atol=1e-9));
        @test isdiag(Hermitian(M, :L), EachApprox(atol=1e-9));
    end

    @testset "isunitary column dot checks" begin
        M = Matrix{Float64}(I, 3, 3);
        M[1, 2] = 1e-10;
        @test isunitary(M, EachApprox(atol=1e-9));
        @test !isunitary(M, Equal());
    end

    @testset "istril with positive k and rectangular" begin
        A = [1 0 0 0; 1e-12 1 0 0; 0 1e-12 1 0];
        @test IsApprox.istril(A, 0, EachApprox(atol=1e-9));
        @test !IsApprox.istril(A, -2, Equal());
    end

    @testset "issymmetric/ishermitian wrong-shape combos" begin
        A = rand(ComplexF64, 2, 3);
        @test !issymmetric(A, Approx());
        @test !ishermitian(A, EachApprox());
    end

    @testset "isposdef/semidef complex small imag" begin
        z = 1 + 1e-12im;
        @test isposdef(z, Approx(atol=1e-10));
        @test !isposdef(z, Equal());
    end
end # @testset "edge cases"
