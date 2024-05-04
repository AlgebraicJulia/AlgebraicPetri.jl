""" Specific generators and useful tools for constructing epidemiological
systems
"""
module Epidemiology

using AlgebraicPetri
using Catlab

export spontaneous_petri, exposure_petri, oapply_epi, 
    infection, exposure, illness, recovery, death

"""
Generates an OpenLabelledPetriNet containing a single transition `z` that moves a
token from place `x` to place `y`
"""
spontaneous_petri(x::Symbol, y::Symbol, z::Symbol) = Open(LabelledPetriNet(unique([x,y]), z=>(x, y)))

"""
Generates an OpenLabelledPetriNet containing a single transition `transition`
which takes as input a single token from place `x` and `y` and produces a single
token in place `z` and `y`.
"""
exposure_petri(x::Symbol, y::Symbol, z::Symbol, transition::Symbol) =
    Open(LabelledPetriNet(unique([x,y,z]), transition=>((x,y)=>(z,y))))

""" LabelledPetriNet which describes the infection process of tokens in state S
by tokens in state I
"""
infection = exposure_petri(:S, :I, :I, :inf)

""" LabelledPetriNet which describes the exposure process where tokens in I
"expose" tokens in S, changing them from S to E
"""
exposure = exposure_petri(:S, :I, :E, :exp)

""" LabelledPetriNet which describes the illness process which moves tokens from E to I.
"""
illness = spontaneous_petri(:E,:I,:ill)

""" LabelledPetriNet which describes the recovery process which moves tokens from I to R
"""
recovery = spontaneous_petri(:I,:R,:rec)

""" LabelledPetriNet which describes the death process which moves tokens from I to D
"""
death = spontaneous_petri(:I,:D,:death)

epi_dict = Dict(:infection=>infection, :exposure=>exposure, :illness=>illness, :recovery=>recovery, :death=>death)

""" oapply_epi(ex, args...)

Generates a LabelledPetriNet under a composition pattern described by the
undirected wiring diagram `ex`. This requires that the boxes in `ex` are only
labelled with labels from the following set:
```
[:infection, :exposure, :illness, :recovery, :death]
```
"""
oapply_epi(ex, args...) = oapply(ex, epi_dict, args...)

end
