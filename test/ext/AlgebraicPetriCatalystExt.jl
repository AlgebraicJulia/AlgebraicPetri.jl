module TestAlgebraicPetriCatalystExt

using Test
using AlgebraicPetri
using Catalyst

# Test with SIR model
sir = PetriNet(3, (1,2)=>(2,2), (2)=>(3))
sir_rxn = ReactionSystem(sir)
@test sir_rxn isa ReactionSystem

end
