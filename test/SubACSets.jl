module TestSubACSets

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

sub_acsets, _ = mca(m1, m2)

@test sub_acsets == [LabelledPetriNet(
    [:X3, :Y3, :W3, :Z3],
    :f3 => (:X3 => :Z3)
  )] 
  #=Set([
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
])=#

@test PetriNet.(sub_acsets) == mca(PetriNet(m1), PetriNet(m2))[1]

m3 = LabelledPetriNet(
  [:W, :X, :Y, :Z],
  :f => ((:W, :X) => (:Y, :Z))
)

sub_acsets2, _ = mca([m3, m2, m1])

@test sub_acsets2 == [LabelledPetriNet(
  [:W, :X, :Y, :Z],
  :f => (:X => :Z)
)]
#=Set([
  LabelledPetriNet(
    [:W, :X, :Y, :Z],
    :f => (:X => :Z)
  ),
  LabelledPetriNet(
    [:W, :X, :Y, :Z],
    :f => (:X => :Y)
  ),
  LabelledPetriNet(
    [:W, :X, :Y, :Z],
    :f => (:W => :Z)
  ),
  LabelledPetriNet(
    [:W, :X, :Y, :Z],
    :f => (:W => :Y)
  )
])=#

end
