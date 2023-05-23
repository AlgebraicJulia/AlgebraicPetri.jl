using Test

@testset "Core" begin
  include("core.jl")
end

@testset "Types" begin
    include("types.jl")
end

@testset "Petri" begin
    include("petri.jl")
end

# TODO: Uncomment once AlgebraicRewriting has been updated
# @testset "Rewriting" begin
#   include("rewriting.jl")
# end
