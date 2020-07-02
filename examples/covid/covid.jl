using Petri
using AlgebraicPetri
using Catlab
using Catlab.Programs
using Catlab.Theories
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.Graphics.Graphviz: Graph

# spontaneous(A::FinOrd, B::FinOrd)           # A -> B
# Petri.Model([1, 2], [([1], [2])])


# transmission(A::FinOrd, B::FinOrd)          # A ⊗ B -> B ⊗ B
# Petri.Model([1, 2], [([1, 2], [2, 2])])


# exposure(A::FinOrd, B::FinOrd, C::FinOrd)   # A ⊗ B -> C ⊗ B
# Petri.Model([1, 2, 3], [([1, 2], [3, 2])])

# SIR   = transmission(S, I) ⋅ ∇(I) ⋅ (spontaneous(I, R)
# SEIR  = exposure(S, I, E) ⋅ (spontaneous(E, I) ⊗ spontaneous(I, R))
# SEIRD = SEIR ⋅ (id(I) ⊗ spontaneous(R, D))

spontaneous = PetriCospan(
        Cospan(FinOrdFunction([1], 1, 2),
               FinOrdFunction([2], 1, 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1), Dict(2=>1))]))

transmission = PetriCospan(
        Cospan(FinOrdFunction([1], 1, 2),
               FinOrdFunction([2], 1, 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1, 2=>1), Dict(2=>2))]))

exposure = PetriCospan(
        Cospan(FinOrdFunction([1, 2], 2, 3),
               FinOrdFunction([3, 2], 2, 3)
        ), id(PetriFunctor), Petri.Model([1, 2, 3], [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]))

SIR = transmission ⋅ spontaneous

Graph(decoration(SIR))

# doesn't work
SEIR = exposure ⋅ (spontaneous ⊗ spontaneous)

Graph(decoration(SEIR))

SEIRD = SEIR ⋅ (id(PetriCospanOb(1)) ⊗ spontaneous)

decoration(SEIRD)

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

sir = transmission ⋅ ∇(I) ⋅ recovery

to_graphviz(sir, orientation=LeftToRight, labels=true)

seir = exposure ⋅ (illness ⊗ recovery)

to_graphviz(seir, orientation=LeftToRight, labels=true)

seird = seir ⋅ (death ⊗ id(R))

to_graphviz(seird, orientation=LeftToRight, labels=true)