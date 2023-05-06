using IsApprox
using Aqua: Aqua

const ThePackage = IsApprox

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(ThePackage)
end

# This often gives false positive
@testset "aqua project toml formatting" begin
    Aqua.test_project_toml_formatting(ThePackage)
end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(ThePackage)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(ThePackage)
end

# TODO: Not sure exactly which versions are ok.
# For <= v1.6, there are ambiguities in Dictionaries
if VERSION >= v"1.7"
    @testset "aqua test ambiguities" begin
        Aqua.test_ambiguities([ThePackage, Core, Base])
    end
end

@testset "aqua piracy" begin
    Aqua.test_piracy(ThePackage)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(ThePackage)
end

@testset "aqua state deps" begin
    Aqua.test_stale_deps(ThePackage)
end
