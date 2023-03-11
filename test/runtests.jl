using Test

using LabelledArrays
using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.ModelComparison
using AlgebraicPetri.OpenTransitions
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Programs
using AlgebraicRewriting

@testset "Core" begin
    include("core.jl")
end

@testset "Types" begin
    include("types.jl")
end

@testset "Petri" begin
    include("petri.jl")
end

@testset "Epidemiology" begin
    include("epidemiology.jl")
end

@testset "TypedPetris" begin
  include("typed_petri.jl")
end

@testset "BilayerNetworks" begin
  include("bilayernetworks.jl")
end

@testset "ModelComparison" begin
  include("modelcomparison.jl")
end

@testset "Catalyst Interop" begin
  include("CatalystInterop.jl")
end

@testset "ModelingToolkit Interop" begin
  include("ModelingToolkitInterop.jl")
end

@testset "TypedPetris" begin
  include("typed_petri.jl")
end

@testset "OpenTransitions" begin
  include("opentransitions.jl")
end

@testset "Rewriting" begin
  include("rewriting.jl")
end