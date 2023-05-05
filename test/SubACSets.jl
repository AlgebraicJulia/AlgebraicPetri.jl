using Test

using AlgebraicPetri, AlgebraicPetri.SubACSets

m1 = LabelledPetriNet(
  [:X3, :Y3, :W3, :Z3],
  :f3 => (:X3 => (:W3, :Y3, :Z3))
)

m2 = LabelledPetriNet(
  [:X4, :W4, :Y4, :Z4],
  :f4 => ((:W4, :X4, :Y4) => :Z4)
)

sub_acsets = mca(m1, m2)

@test sub_acsets == Set([
  LabelledPetriNet(
    [:X3, :Y3, :W3, :Z3],
    :f3 => (:X3 => :Z3)
  ),
  LabelledPetriNet(
    [:X3, :Y3, :W3, :Z3],
    :f3 => (:X3 => :W3)
  ),
  LabelledPetriNet(
    [:X3, :Y3, :W3, :Z3],
    :f3 => (:X3 => :Y3)
  )
])

@test Set(PetriNet.(sub_acsets)) == mca(PetriNet(m1), PetriNet(m2))
