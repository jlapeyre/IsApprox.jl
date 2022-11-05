module IsApprox

import DictTools

export AbstractApprox, Equal, EachApprox, Approx, UpToPhase
export ispossemidef, isunitary, isinvolution, isidempotent, isnormal, commutes, anticommutes
export isnormalized, isprobdist

include("core.jl")
include("base_applications.jl")
include("other_applications.jl")

end # module
