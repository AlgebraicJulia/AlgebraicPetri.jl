sir_petri = LabelledPetriNet([:S,:I,:R], :inf=>((:S,:I)=>(:I,:I)), :rec=>(:I=>:R))

sir = transmission ⋅ recovery

@test apex(F_epi(sir)) == sir_petri
