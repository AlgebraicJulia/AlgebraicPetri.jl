using Test

@testset "Catalyst Package Extension" begin
  include("AlgebraicPetriCatalystExt.jl")
end

@testset "ModelingToolkit Package Extension" begin
  include("AlgebraicPetriModelingToolkitExt.jl")
end
