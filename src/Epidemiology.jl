module Epidemiology

using AlgebraicPetri
using Catlab
using Catlab.Theories
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra.FinSets

export oapply_epi, infection, exposure, illness, recovery, death

spontaneous_petri(x::Symbol, y::Symbol, z::Symbol) = Open(LabelledPetriNet(unique([x,y]), z=>(x, y)))
exposure_petri(x::Symbol, y::Symbol, z::Symbol, transition::Symbol) =
    Open(LabelledPetriNet(unique([x,y,z]), transition=>((x,y)=>(z,y))))

infection = exposure_petri(:S, :I, :I, :inf)
exposure = exposure_petri(:S, :I, :E, :exp)
illness = spontaneous_petri(:E,:I,:ill)
recovery = spontaneous_petri(:I,:R,:rec)
death = spontaneous_petri(:I,:D,:death)

epi_dict = Dict(:infection=>infection, :exposure=>exposure, :illness=>illness, :recovery=>recovery, :death=>death)

""" oapply_epi
"""
oapply_epi(ex, args...) = oapply(ex, epi_dict, args...)

end
