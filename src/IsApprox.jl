module IsApprox

using PrecompileTools: @setup_workload, @compile_workload

using Dictionaries: Dictionaries
export AbstractApprox, Equal, EachApprox, Approx, UpToPhase
export isposdef, ispossemidef, isunitary, isinvolution, isidempotent, isnormal, commutes, anticommutes
export isnormalized, isprobdist
# From LinearAlgebra
export ishermitian, issymmetric, istriu, istril, isbanded, isdiag

# This is from DictTools.jl which is not yet registered
"""
    _AbstractDict{T, V}

Either an `AbstractDict` or an `AbstractDictionary`. A union type
"""
const _AbstractDict{T, V} = Union{AbstractDict{T,V}, Dictionaries.AbstractDictionary{T,V}}

include("core.jl")
include("base_applications.jl")
include("other_applications.jl")
include("precompile.jl")

end # module
