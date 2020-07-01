using Test

using Petri
using AlgebraicPetri
using Catlab.Theories
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets

dom_list(p::PetriCospan) = left(p).(1:dom(left(p)).n)
codom_list(p::PetriCospan) = right(p).(1:dom(right(p)).n)

function compare_petricospan(p₁::PetriCospan, p₂::PetriCospan)
    @test dom(p₁) == dom(p₂)
    @test codom(p₁) == codom(p₂)
    @test dom_list(p₁) == dom_list(p₂)
    @test codom_list(p₁) == codom_list(p₂)
    @test decoration(p₁) == decoration(p₂)
end

@testset "Core" begin
    include("core.jl")
end

@testset "Petri" begin
    include("petri.jl")
end