p1 = PetriCospanOb(1)
p2 = PetriCospanOb(2)
pc1 = id(p1)
pc2 = id(p2)
I = munit(PetriCospanOb)

@test dom(pc1 ⊗ pc2) == dom(pc1) ⊗ dom(pc2)
@test codom(pc1 ⊗ pc2) == codom(pc1) ⊗ codom(pc2)
@test dom(σ(p1,p2)) == p1 ⊗ p2
@test codom(σ(p1,p2)) == p2 ⊗ p1
compare_petricospan(σ(p1, p2), braid(p1,p2))

@test dom(Δ(p1)) == p1
@test codom(Δ(p1)) == p1 ⊗ p1
@test dom(◊(p1)) == p1
@test codom(◊(p1)) == I

compare_petricospan(pair(pc1,pc1), Δ(p1) ⋅ (pc1 ⊗ pc1))
compare_petricospan(copair(pc1,pc1), (pc1 ⊗ pc1) ⋅ ∇(p1))

compare_petricospan(proj1(p1,p2), pc1 ⊗ ◊(p2))
compare_petricospan(proj2(p1,p2), ◊(p1) ⊗ pc2)

compare_petricospan(coproj1(p1,p2), pc1 ⊗ □(p2))
compare_petricospan(coproj2(p1,p2), □(p1) ⊗ pc2)
