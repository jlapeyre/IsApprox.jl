# The call `ishermitian(A::AbstractMatrix, B::AbstractMatrix) lowers
# to exactly the same code as that in LinearAlgebra
function ishermitian(A::AbstractMatrix, approx::AbstractApprox=Equal())
    indsm, indsn = axes(A)
    if indsm != indsn
        return false
    end
    for i = indsn, j = i:last(indsn)
        if ! isapprox(approx, A[i,j], adjoint(A[j,i]))
            return false
        end
    end
    return true
end

# This uses the isapprox interface, which compares using a norm
ishermitian(A::AbstractMatrix, approx::Approx) = isapprox(approx, A, adjoint(A))

ishermitian(x::Number, approx::AbstractApprox=Equal()) = isapprox(approx, x, conj(x))
