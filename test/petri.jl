f = PetriCospan([1, 2], Petri.Model(1:4, [(Dict(1=>1), Dict(3=>1)), (Dict(2=>1), Dict(4=>1))]), [3, 4])

g = PetriCospan([1,2], Petri.Model(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1))]), [3])

h = f ⋅ g

h′ = PetriCospan( [1,2], Petri.Model(1:5, [(Dict(1=>1), Dict(3=>1)), (Dict(2=>1), Dict(4=>1)), (Dict(3=>1, 4=>1), Dict(5=>1))]), [5])

h_id = h ⋅ id(PetriCospanOb(1))

@test dom(f) == PetriCospanOb(2)
@test codom(f) == PetriCospanOb(2)
@test dom(g) == codom(f)
@test codom(g) == PetriCospanOb(1)
@test dom(h) == dom(f)
@test codom(h) == codom(g)

compare_petricospan(h, h′)
compare_petricospan(h, h_id)