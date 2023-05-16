using Test

@testset "AlgebraicPetri" begin
    include("algebraicpetri/AlgebraicPetri.jl")
end

@testset "BilayerNetworks" begin
  include("BilayerNetworks.jl")
end

@testset "Epidemiology" begin
    include("Epidemiology.jl")
end

@testset "ModelComparison" begin
  include("ModelComparison.jl")
end

@testset "OpenTransitions" begin
  include("OpenTransitions.jl")
end

@testset "TypedPetris" begin
  include("TypedPetri.jl")
end

@testset "Package Extensions" begin
  include("ext/extensions.jl")
end
