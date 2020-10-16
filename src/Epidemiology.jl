module Epidemiology

using AlgebraicPetri
using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra.FinSets

export InfectiousDiseases, FunctorGenerators, F_epi, S, E, I, R, D, transmission, exposure, illness, recovery, death

ob(x::Symbol) = codom(Open([x], LabelledPetriNet(x), [x]))
spontaneous_petri(x::Symbol, y::Symbol, z::Symbol) = Open([x], LabelledPetriNet([x,y], z=>(x, y)), [y])
transmission_petri = Open([:S], LabelledPetriNet([:S,:I], :inf=>((:S,:I)=>(:I,:I))), [:I])
exposure_petri = Open([:S, :I], LabelledPetriNet([:S,:I,:E], :exp=>((:S,:I)=>(:E,:I))), [:E])

""" InfectiousDiseases
"""
@present InfectiousDiseases(FreeBiproductCategory) begin
    S::Ob
    E::Ob
    I::Ob
    R::Ob
    D::Ob
    transmission::Hom(S⊗I, I)
    exposure::Hom(S⊗I, E)
    illness::Hom(E,I)
    recovery::Hom(I,R)
    death::Hom(I,D)
end

""" generators
"""
S,E,I,R,D,transmission,exposure,illness,recovery,death = generators(InfectiousDiseases);

""" FunctorGenerators
"""
const FunctorGenerators = Dict(S=>ob(:S), E=>ob(:E), I=>ob(:I), R=>ob(:R), D=>ob(:D),
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri(:E,:I,:ill), recovery=>spontaneous_petri(:I,:R,:rec), death=>spontaneous_petri(:I,:D,:death))

""" F_epi
"""
F_epi(ex) = functor((OpenLabelledPetriNetOb, OpenLabelledPetriNet), ex, generators=FunctorGenerators)

end

# sir = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .3/1000)=>((:S, :I)=>(:I,:I)), (:rec, .2)=>(:I=>:R))
# siir = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :A=>0, :R=>0), (:inf_i, .3/1000)=>((:S, :I)=>(:I,:I)), (:inf_a, .0003)=>((:S,:A)=>(:I,:A)), (:asym_i, .0003)=>((:S,:I)=>(:I,:A)), (:asym_a, .0003)=>((:S,:A)=>(:A,:A)),(:rec_i, .2)=>(:I=>:R),(:rec_a, .2)=>(:A=>:R))
# seair = LabelledReactionNet{Number, Int}((:S=>990, :E=>0, :I=>10, :A=>0, :R=>0), (:exp_i, .3/1000)=>((:S, :I)=>(:I,:E)), (:exp_a, .3/1000)=>((:S, :A)=>(:A, :E)), (:inf, .0003)=>(:E=>:I), (:asym, .0003)=>(:E=>:A), (:rec_a, .2)=>(:A=>:R), (:rec_i, .2)=>(:I=>:R))