using Test

using Petri
using LabelledArrays
using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using Catlab.Theories
using Catlab.Present
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets

@testset "Core" begin
    include("core.jl")
end

@testset "Types" begin
    include("types.jl")
end

@testset "Visualization" begin
    include("visualization.jl")
end

@testset "Petri" begin
    include("petri.jl")
end

@testset "Epidemiology" begin
    include("epidemiology.jl")
end
