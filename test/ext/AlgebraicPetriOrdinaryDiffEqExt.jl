module TestAlgebraicPetriOrdinaryDiffEqExt

using Test
using AlgebraicPetri
using OrdinaryDiffEq

sir_petri = PetriNet(3, ((1, 2), (2, 2)), (2, 3))
sir_rxn = ReactionNet{Number,Int}([990, 10, 0], (0.001, ((1, 2) => (2, 2))), (0.25, (2 => 3)))

sol_petri = solve(ODEProblem(sir_petri, [990, 10, 0], (0,100.0), [0.001, 0.25]), Tsit5())
sol_rxn = solve(ODEProblem(sir_rxn, (0,100.0)), Tsit5())

@test sol_petri.u == sol_rxn.u

end
