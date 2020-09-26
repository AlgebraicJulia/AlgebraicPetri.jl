sir_petri = PetriNet(3, ((1, 2), (2, 2)), (2, 3))
sir_lpetri = LabelledPetriNet([:S, :I, :R], :inf=>((:S, :I), (:I, :I)), :rec=>(:I, :R))
β(u,t) = 1 / sum(u)
γ = .25
sir_rxn = ReactionNet{Function, Int}([990, 10, 0], (β)=>((1, 2)=>(2,2)), (t->γ)=>(2=>3))
sir_lrxn = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .001)=>((:S, :I)=>(:I,:I)), (:rec, .25)=>(:I=>:R))

sir_tpetri= PetriNet(TransitionMatrices(sir_petri))

@test sir_tpetri == sir_petri
@test Petri.Model(sir_petri) == Petri.Model(sir_rxn)
@test Petri.Model(sir_lpetri) == Petri.Model(sir_lrxn)

@test concentration(sir_rxn, 1) == 990
@test rate(sir_rxn, 1) == β

@test concentrations(sir_rxn) == [990, 10, 0]
@test typeof(rates(sir_rxn)) <: Array{Function}

@test concentrations(sir_lrxn) == LVector(S=990, I=10, R=0)
@test rates(sir_lrxn) == LVector(inf=.001, rec=.25)

du = [0.0, 0.0, 0.0]
out = vectorfield(sir_rxn)(du, concentrations(sir_rxn), rates(sir_rxn), 0.01)
@test out[1] ≈ -9.9
@test out[2] ≈ 7.4
@test out[3] ≈ 2.5

du = LVector(S=0.0, I=0.0, R=0.0)
out = vectorfield(sir_lrxn)(du, concentrations(sir_lrxn), rates(sir_lrxn), 0.01)
@test out.S ≈ -9.9
@test out.I ≈ 7.4
@test out.R ≈ 2.5

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