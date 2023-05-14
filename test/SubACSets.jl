module TestSubACSets

using Test

using AlgebraicPetri, AlgebraicPetri.SubACSets
using AlgebraicPetri.SubACSets: strip_attributes

m1 = LabelledPetriNet(
  [:X1, :Y1, :W1, :Z1],
  :f1 => (:X1 => (:W1, :Y1, :Z1))
)

m2 = LabelledPetriNet(
  [:X2, :W2, :Y2, :Z2],
  :f2 => ((:W2, :X2, :Y2) => :Z2)
)

mca1, mca1_morphs = mca(m1, m2)

@test is_isomorphic(mca1,strip_attributes(LabelledPetriNet(
  [:X1, :Y1, :W1, :Z1],
  :f1 => (:X1 => :Z1)
)))
@test length(mca1_morphs) == 2
@test length(mca1_morphs[1]) == 3
@test length(mca1_morphs[2]) == 3


@test is_isomorphic(PetriNet(mca1),mca(PetriNet(m1), PetriNet(m2))[1])

m3 = LabelledPetriNet(
  [:W3, :X3, :Y3, :Z3],
  :f3 => ((:W3, :X3) => (:Y3, :Z3))
)

mca3, mca3_morphs = mca([m3, m2, m1])
@test is_isomorphic(mca3,strip_attributes(LabelledPetriNet(
  [:W3, :X3, :Y3, :Z3],
  :f3 => (:X3 => :Z3)
)))
@test length(mca3_morphs) == 3
@test length(mca3_morphs[1]) == 4
@test length(mca3_morphs[2]) == 3
@test length(mca3_morphs[3]) == 3

m4 = LabelledPetriNet(
  [:X4, :Y4],
  :f41 => (:X4 => :Y4),
  :f42 => (:Y4 => :Y4)  
)
mca4, mca4_morphs = mca([m3, m2, m1, m4])
@test is_isomorphic(mca4,strip_attributes(LabelledPetriNet(
  [:X4, :Y4],
  :f41 => (:X4 => :Y4)
)))
@test length(mca4_morphs) == 4
@test length(mca4_morphs[1]) == 4
@test length(mca4_morphs[2]) == 3
@test length(mca4_morphs[3]) == 3
@test length(mca4_morphs[4]) == 2

end
