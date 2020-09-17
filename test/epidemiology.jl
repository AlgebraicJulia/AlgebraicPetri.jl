sir_petri = Petri.Model(3, Dict(1=>(Dict(1=>1,2=>1), Dict(2=>2)), 2=>(Dict(2=>1), Dict(3=>1))))

sir = transmission ⋅ recovery

@test Petri.Model(decoration(F_epi(sir))) == sir_petri
