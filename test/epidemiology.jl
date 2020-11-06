sir_petri = LabelledPetriNet([:S,:I,:R], :inf=>((:S,:I)=>(:I,:I)), :rec=>(:I=>:R))

sir = infection ⋅ recovery

@test apex(sir) == sir_petri

sir_relation = @relation (s,i,r) begin
    infection(s,i)
    recovery(i,r)
end

@test apex(oapply_epi(sir_relation)) == sir_petri