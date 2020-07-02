using Petri
using AlgebraicPetri
using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.Graphics.Graphviz: Graph

ob = PetriCospanOb(1)

spontaneous_petri = PetriCospan(
        Cospan(FinOrdFunction([1], 1, 2),
               FinOrdFunction([2], 1, 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1), Dict(2=>1))]))

transmission_petri = PetriCospan(
        Cospan(FinOrdFunction([1], 1, 2),
               FinOrdFunction([2], 1, 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1, 2=>1), Dict(2=>2))]))

exposure_petri = PetriCospan(
        Cospan(FinOrdFunction([1, 2], 2, 3),
               FinOrdFunction([3, 2], 2, 3)
        ), id(PetriFunctor), Petri.Model([1, 2, 3], [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]))

@present Epidemiology(FreeBiproductCategory) begin
    S::Ob
    E::Ob
    I::Ob
    R::Ob
    D::Ob
    transmission::Hom(otimes(S,I), I)
    exposure::Hom(otimes(S,I), otimes(E,I))
    illness::Hom(E,I)
    recovery::Hom(I,R)
    death::Hom(I,D)
end

S,E,I,R,D,transmission,exposure,illness,recovery,death = generators(Epidemiology)

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
        S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri))

sir = transmission ⋅ recovery
f_sir = F(sir)

to_graphviz(sir, orientation=LeftToRight, labels=true)
Graph(decoration(f_sir))

sei = exposure ⋅ (illness ⊗ id(I)) ⋅ ∇(I)
seir = sei ⋅ recovery
f_seir = F(seir)

to_graphviz(seir, orientation=LeftToRight, labels=true)
Graph(decoration(f_seir))

seird = sei ⋅ Δ(I) ⋅ (death ⊗ recovery)
f_seird = F(seird)

to_graphviz(seird, orientation=LeftToRight, labels=true)
Graph(decoration(f_seird))

# TODO: Add support for types so we can simplify
# seir = exposure ⋅ (illness ⊗ recovery)
# seird = seir ⋅ (death ⊗ id(R))