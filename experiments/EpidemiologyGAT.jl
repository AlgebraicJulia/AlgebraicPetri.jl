using AlgebraicPetri
using Petri
using Catlab
using Catlab.GAT
using Catlab.Theories
using Catlab.Programs
using Catlab.WiringDiagrams

using Catlab.CategoricalAlgebra.ShapeDiagrams # remove later
using Catlab.Graphics # remove later

import Catlab.Theories: dom, codom, id, compose, ⋅, ∘, otimes, ⊗, munit,
                        braid, σ, mcopy, Δ, mmerge, ∇, create, □, delete, ◊,
                        pair, copair, proj1, proj2, coproj1, coproj2
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

@theory BiproductCategory(Ob,Hom) => EpidemiologyC(Ob,Hom) begin
  spontaneous(A::Ob, B::Ob)::(A → B)
  transmission(A::Ob, B::Ob)::(A ⊗ B → B)
  exposure(A::Ob, B::Ob, C::Ob)::(A ⊗ B → C ⊗ B)
end

@instance EpidemiologyC(PetriCospanOb, PetriCospan) begin
  @import dom, codom, compose, id, otimes, munit, braid, mcopy, mmerge, create, delete, pair, copair, proj1, proj2, coproj1, coproj2
  spontaneous(A::PetriCospanOb, B::PetriCospanOb) = PetriCospan([1], Petri.Model(1:2, [(Dict(1=>1), Dict(2=>1))]), [2])
  transmission(A::PetriCospanOb, B::PetriCospanOb) = PetriCospan([1,2], Petri.Model(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))]), [2])
  exposure(A::PetriCospanOb, B::PetriCospanOb, C::PetriCospanOb) = PetriCospan([1, 2], Petri.Model(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]), [3, 2])
end

@present BasicEpi(FreeEpidemiologyC) begin
  S::Ob
  E::Ob
  I::Ob
  R::Ob
  D::Ob
end

S,E,I,R,D = generators(BasicEpi)

(t::Presentation)(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(map(x->x=>PetriCospanOb(1), generators(t))))

sir = transmission(S,I) ⋅ spontaneous(I,R)

display_wd(sir)

BasicEpi(sir)

sei = exposure(S,I,E) ⋅ (spontaneous(E,I) ⊗ id(I)) ⋅ ∇(I)

seir = sei ⋅ spontaneous(I,R)

BasicEpi(seir)

display_wd(seir)

seird = sei ⋅ Δ(I) ⋅ (spontaneous(I,D) ⊗ spontaneous(I,R))

BasicEpi(seird)

display_wd(seird)