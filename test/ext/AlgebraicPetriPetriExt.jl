module TestAlgebraicPetriPetriExt

using Test
using AlgebraicPetri
import Petri

sir_petri = PetriNet(3, ((1, 2), (2, 2)), (2, 3))
sir_lpetri = LabelledPetriNet([:S, :I, :R], :inf => ((:S, :I), (:I, :I)), :rec => (:I, :R))
sir_rxn = ReactionNet{Number,Int}([990, 10, 0], (0.001, ((1, 2) => (2, 2))), (0.25, (2 => 3)))
sir_lrxn = LabelledReactionNet{Number,Int}((:S => 990, :I => 10, :R => 0), (:inf, 0.001) => ((:S, :I) => (:I, :I)), (:rec, 0.25) => (:I => :R))

@test Petri.Model(sir_petri) == Petri.Model(sir_rxn)
@test Petri.Model(sir_lpetri) == Petri.Model(sir_lrxn)

end
