f = PetriCospan([1, 2], PetriNet(4, (1,3), (2,4)), [3, 4])

g = PetriCospan([1,2], PetriNet(3, ((1,2),3)), [3])

h = f ⋅ g

h′ = PetriCospan([1,2], PetriNet(5, (1,3), (2,4), ((3,4), 5)), [5])

h_id = h ⋅ id(PetriCospanOb(1))

@test dom(f) == PetriCospanOb(2)
@test codom(f) == PetriCospanOb(2)
@test dom(g) == codom(f)
@test codom(g) == PetriCospanOb(1)
@test dom(h) == dom(f)
@test codom(h) == codom(g)

compare_petricospan(h, h′)
compare_petricospan(h, h_id)
