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

# Test open petri net notations either fully open or by specifying legs

pn = Open(PetriNet(3, ((1,2), 3)), [1], [2], [3])
pn′ = Open(PetriNet(3, ((1,2), 3)))
@test pn == pn′

lpn = Open(LabelledPetriNet([:I, :R], (:rec, :I=>:R)))
lpn′ = Open([:I], LabelledPetriNet([:I, :R], (:rec, :I=>:R)), [:R])
@test lpn == lpn′

rn = Open(ReactionNet{Number,Int}([10, 0], (.25, 1=>2)))
rn′ = Open([1], ReactionNet{Number,Int}([10, 0], (.25, 1=>2)), [2])
@test rn == rn′

lrn = Open(LabelledReactionNet{Number,Int}([:I=>10, :R=>0], ((:rec=>.25), :I=>:R)))
lrn′ = Open([:I], LabelledReactionNet{Number,Int}([:I=>10, :R=>0], ((:rec=>.25), :I=>:R)), [:R])
@test lrn == lrn′