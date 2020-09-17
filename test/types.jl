sir_petri = PetriNet(3, ((1, 2), (2, 2)), (2, 3))
sir_lpetri = LabelledPetriNet([:S, :I, :R], :inf=>((:S, :I), (:I, :I)), :rec=>(:I, :R))
sir_rxn = ReactionNet{Number, Int}([990, 10, 0], (.0001)=>((1, 2)=>(2,2)), (.25)=>(2=>3))
sir_lrxn = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .0001)=>((:S, :I)=>(:I,:I)), (:rec, .25)=>(:I=>:R))

sir_tpetri= PetriNet(TransitionMatrices(sir_petri))

@test sir_tpetri == sir_petri
@test Petri.Model(sir_petri) == Petri.Model(sir_rxn)
@test Petri.Model(sir_lpetri) == Petri.Model(sir_lrxn)

@test concentrations(sir_rxn) == [990, 10, 0]
@test rates(sir_rxn) == [.0001, .25]

@test concentrations(sir_lrxn) == LVector(S=990, I=10, R=0)
@test rates(sir_lrxn) == LVector(inf=.0001, rec=.25)

@test ns(sir_petri) == 3
add_species!(sir_petri)
@test ns(sir_petri) == 4

@test nt(sir_petri) == 2
add_transitions!(sir_petri, 2)
@test nt(sir_petri) == 4

@test ni(sir_petri) == 3
@test no(sir_petri) == 3
add_input!(sir_petri, 3, 1)
add_output!(sir_petri, 3, 4)
add_input!(sir_petri, 4, 4)
add_output!(sir_petri, 4, 3)
@test ni(sir_petri) == 5
@test no(sir_petri) == 5
@test sir_petri == PetriNet(4, ((1, 2), (2, 2)), (2, 3), (1, 4), (4, 3))
