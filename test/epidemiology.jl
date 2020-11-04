sir_petri = LabelledPetriNet([:S,:I,:R], :inf=>((:S,:I)=>(:I,:I)), :rec=>(:I=>:R))

sir = infection ⋅ recovery

@test apex(sir) == sir_petri
