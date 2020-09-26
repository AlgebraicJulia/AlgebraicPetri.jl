f = Open([1, 2], PetriNet(4, (1,3), (2,4)), [3, 4])

g = Open([1,2], PetriNet(3, ((1,2),3)), [3])

h = f ⋅ g

h′ = Open([1,2], PetriNet(5, (1,3), (2,4), ((3,4), 5)), [5])

h_id = h ⋅ id(OpenPetriNetOb(FinSet(1)))

@test dom(f) == OpenPetriNetOb(FinSet(2))
@test codom(f) == OpenPetriNetOb(FinSet(2))
@test dom(g) == codom(f)
@test codom(g) == OpenPetriNetOb(FinSet(1))
@test dom(h) == dom(f)
@test codom(h) == codom(g)

@test h == h′
@test h == h_id
