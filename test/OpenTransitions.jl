module TestOpenTransitions

using Test
using AlgebraicPetri
using AlgebraicPetri.OpenTransitions
using Catlab, Catlab.CategoricalAlgebra, Catlab.Programs, Catlab.WiringDiagrams

# tests with unlabeled PN

# merge along 1 transition, specify legs
annihilation = OpenT(PetriNet(1, 1=>()), [1])
creation = OpenT(PetriNet(1, ()=>1), [1])

join_pattern = @relation (T,) begin
  X(T)
  Y(T)
end

joined_pn = oapply(join_pattern,
  Dict(
    :X => annihilation,
    :Y => creation
  )
)

@test ns(apex(joined_pn)) == 2
@test nt(apex(joined_pn)) == 1
@test subpart(apex(joined_pn), :it) == [1]
@test subpart(apex(joined_pn), :is) == [1]
@test subpart(apex(joined_pn), :ot) == [1]
@test subpart(apex(joined_pn), :os) == [2]
@test nparts(apex(joined_pn), :O) == 1
@test nparts(apex(joined_pn), :I) == 1

# merge along 1 transition, do not specify legs
annihilation = OpenT(PetriNet(1, 1=>()))
creation = OpenT(PetriNet(1, ()=>1))

joined_pn = oapply(join_pattern,
  Dict(
    :X => annihilation,
    :Y => creation
  )
)

@test ns(apex(joined_pn)) == 2
@test nt(apex(joined_pn)) == 1
@test subpart(apex(joined_pn), :it) == [1]
@test subpart(apex(joined_pn), :is) == [1]
@test subpart(apex(joined_pn), :ot) == [1]
@test subpart(apex(joined_pn), :os) == [2]
@test nparts(apex(joined_pn), :O) == 1
@test nparts(apex(joined_pn), :I) == 1

# merge along 2 transitions
annihilation = OpenT(PetriNet(1, 1=>(), 1=>()), [1], [2])
creation = OpenT(PetriNet(1, ()=>1, ()=>1), [1], [2])

join_pattern = @relation (T1,T2) begin
  X(T1,T2)
  Y(T1,T2)
end

joined_pn = oapply(join_pattern,
  Dict(
    :X => annihilation,
    :Y => creation
  )
)

@test length(legs(joined_pn)) == 2
@test ns(apex(joined_pn)) == 2
@test nt(apex(joined_pn)) == 2
@test subpart(apex(joined_pn), :it) == [1,2]
@test subpart(apex(joined_pn), :is) == [1,1]
@test subpart(apex(joined_pn), :ot) == [1,2]
@test subpart(apex(joined_pn), :os) == [2,2]

# tests with labeled PN

# merge along 1 transition, specify legs
annihilation = OpenT(
  LabelledPetriNet([:X],
    :T => (:X => ()) 
  ), [:T]
)

creation = OpenT(
  LabelledPetriNet([:X],
    :T => (() => :X)
  ), [:T]
)

join_pattern = @relation (T,) begin
  X(T)
  Y(T)
end

joined_pn = oapply(join_pattern,
  Dict(
    :X => annihilation,
    :Y => creation
  )
)

@test apex(joined_pn) isa LabelledPetriNet
@test ns(apex(joined_pn)) == 2
@test nt(apex(joined_pn)) == 1
@test subpart(apex(joined_pn), :it) == [1]
@test subpart(apex(joined_pn), :is) == [1]
@test subpart(apex(joined_pn), :ot) == [1]
@test subpart(apex(joined_pn), :os) == [2]
@test nparts(apex(joined_pn), :O) == 1
@test nparts(apex(joined_pn), :I) == 1

# merge along 1 transition, do not specify legs
annihilation = OpenT(
  LabelledPetriNet([:X],
    :T => (:X => ()) 
  )
)

creation = OpenT(
  LabelledPetriNet([:X],
    :T => (() => :X)
  )
)

joined_pn = oapply(join_pattern,
  Dict(
    :X => annihilation,
    :Y => creation
  )
)

@test apex(joined_pn) isa LabelledPetriNet
@test ns(apex(joined_pn)) == 2
@test nt(apex(joined_pn)) == 1
@test subpart(apex(joined_pn), :it) == [1]
@test subpart(apex(joined_pn), :is) == [1]
@test subpart(apex(joined_pn), :ot) == [1]
@test subpart(apex(joined_pn), :os) == [2]
@test nparts(apex(joined_pn), :O) == 1
@test nparts(apex(joined_pn), :I) == 1

# merge along 2 transitions
annihilation = OpenT(
  LabelledPetriNet([:X], 
    :T1 => (:X=>()), 
    :T2 => (:X=>())
  ), [:T1], [:T2]
)
creation = OpenT(
  LabelledPetriNet([:X], 
    :T1 => (()=>:X), 
    :T2 => (()=>:X)
  ), [:T1], [:T2]
)

join_pattern = @relation (T1,T2) begin
  X(T1,T2)
  Y(T1,T2)
end

joined_pn = oapply(join_pattern,
  Dict(
    :X => annihilation,
    :Y => creation
  )
)

@test apex(joined_pn) isa LabelledPetriNet
@test length(legs(joined_pn)) == 2
@test ns(apex(joined_pn)) == 2
@test nt(apex(joined_pn)) == 2
@test subpart(apex(joined_pn), :it) == [1,2]
@test subpart(apex(joined_pn), :is) == [1,1]
@test subpart(apex(joined_pn), :ot) == [1,2]
@test subpart(apex(joined_pn), :os) == [2,2]

end
