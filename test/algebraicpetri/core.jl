module TestCore

using Test
using AlgebraicPetri
using Catlab.Theories
using Catlab.CategoricalAlgebra

p1 = codom(Open([1], PetriNet(1), [1]))
p2 = codom(Open([1,2], PetriNet(2), [1,2]))
pc1 = id(p1)
pc2 = id(p2)

@test dom(pc1 ⊗ pc2) == dom(pc1) ⊗ dom(pc2)
@test codom(pc1 ⊗ pc2) == codom(pc1) ⊗ codom(pc2)
@test dom(braid(p1,p2)) == p1 ⊗ p2
@test codom(braid(p1,p2)) == p2 ⊗ p1

@test dom(mcopy(p1)) == p1
@test codom(mcopy(p1)) == p1 ⊗ p1
@test dom(delete(p1)) == p1
@test codom(delete(p1)) == munit(OpenPetriNetOb)

end
