module IsApprox

using Dictionaries: Dictionaries
export AbstractApprox, Equal, EachApprox, Approx, UpToPhase
export isposdef, ispossemidef, isunitary, isinvolution, isidempotent, isnormal, commutes, anticommutes
export isnormalized, isprobdist

# This is from DictTools.jl which is not yet registered
"""
    _AbstractDict{T, V}

Either an `AbstractDict` or an `AbstractDictionary`. A union type
"""
const _AbstractDict{T, V} = Union{AbstractDict{T,V}, Dictionaries.AbstractDictionary{T,V}}

include("core.jl")
include("base_applications.jl")
include("other_applications.jl")

end # module
