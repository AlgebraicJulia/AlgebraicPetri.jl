module Epidemiology

using AlgebraicPetri
using Petri
using Catlab
using Catlab.Theories

export InfectiousDiseases, FunctorGenerators, F_epi, S, E, I, R, D, transmission, exposure, illness, recovery, death

ob = PetriCospanOb(1)
spontaneous_petri = PetriCospan([1], PetriWithRates(1:2, [(Dict(1=>1), Dict(2=>1))], [0], [0,0]), [2])
transmission_petri = PetriCospan([1,2], PetriWithRates(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))], [0], [0,0]), [2])
exposure_petri = PetriCospan([1, 2], PetriWithRates(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))], [0], [0,0,0]), [3])

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
const FunctorGenerators = Dict(S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri)

""" F_epi
"""
F_epi(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=FunctorGenerators)

end