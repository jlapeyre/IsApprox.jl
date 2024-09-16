@setup_workload begin
    nothing
    using IsApprox
    import LinearAlgebra
    @compile_workload begin
        for T in (Float64, Float32)
            for m in (rand(T, 4, 4), rand(Complex{T}, 2, 2))
                for approx in (Equal(), Approx(), EachApprox())
                    IsApprox.isposdef(m, approx)
                    IsApprox.ispossemidef(m, approx)
                    IsApprox.ishermitian(m, approx)
                    IsApprox.issymmetric(m, approx)
                    IsApprox.iszero(m, approx)
                    IsApprox.isone(m, approx)
                    IsApprox.isreal(m, approx)
                    IsApprox.isdiag(m, approx)
                    IsApprox.isunitary(m, approx)
                    IsApprox.isidempotent(m, approx)
                    IsApprox.isnormalized(m, approx)
                end
            end
        end
    end
end
