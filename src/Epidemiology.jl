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
oapply_epi(ex) = oapply(ex, epi_dict)

end

# sir = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .3/1000)=>((:S, :I)=>(:I,:I)), (:rec, .2)=>(:I=>:R))
# siir = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :A=>0, :R=>0), (:inf_i, .3/1000)=>((:S, :I)=>(:I,:I)), (:inf_a, .0003)=>((:S,:A)=>(:I,:A)), (:asym_i, .0003)=>((:S,:I)=>(:I,:A)), (:asym_a, .0003)=>((:S,:A)=>(:A,:A)),(:rec_i, .2)=>(:I=>:R),(:rec_a, .2)=>(:A=>:R))
# seair = LabelledReactionNet{Number, Int}((:S=>990, :E=>0, :I=>10, :A=>0, :R=>0), (:exp_i, .3/1000)=>((:S, :I)=>(:I,:E)), (:exp_a, .3/1000)=>((:S, :A)=>(:A, :E)), (:inf, .0003)=>(:E=>:I), (:asym, .0003)=>(:E=>:A), (:rec_a, .2)=>(:A=>:R), (:rec_i, .2)=>(:I=>:R))