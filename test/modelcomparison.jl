SIRD = LabelledReactionNet{Float64, Float64}([:S=>0.0, :I=>0.0, :R=>0.0, :D=>0.0],
                                             (:inf=>0.0)=>((:S,:I)=>(:I,:I)),
                                             (:rec=>0.0)=>(:I=>:R),
                                             (:death=>0.0)=>(:I=>:D))

SIR  = LabelledReactionNet{Float64, Float64}([:S=>1.0, :I=>0.0, :R=>0.0],
                                             (:inf=>0.5)=>((:S,:I)=>(:I,:I)),
                                             (:rec=>0.1)=>(:I=>:R))
models = [SIR, SIRD]


AlgebraicPetri.ModelComparison.compare(A::Subobject, B::Subobject) = compare(dom(hom(A)), dom(hom(B)))
@testset "Subobjects" begin
for pn in [models,
           ReactionNet{Float64, Float64}.([SIR, SIRD]),
           LabelledPetriNet.(models),
           PetriNet.(models)]

  c_res = compare(pn[1], pn[2])
  @test apex(c_res) == pn[1]
  @test length(legs(c_res)) == 2

  @test Graph(legs(c_res)[1]) isa Graph
  @test Graph(c_res) isa Graph

  so = Subobject.(legs(c_res))
  for s in so
    @test Graph(s) isa Graph
    @test codom(hom(s)) == pn[2]
  end

  @test dom(hom(~foldl(∨, so) ∨ foldl(∨, so))) == pn[2]
  A,B = so
  @test implies(A, B) == ¬(A) ∨ B
  @test ¬(A ∧ B) == ¬(A) ∨ ¬(B)
  @test ¬(A ∧ B) != ¬(A) ∨ B
  # modus ponens holds only up to inclusion, not equality.
  @test length(compare(A ∧ implies(A,B), B)) > 0
  # this is an equivalent check because X ∧ B == X iff X ↪ B 
  @test (A ∧ implies(A,B)) == B ∧ (A ∧ implies(A,B))

  @test length(compare(B ∧ implies(B,A), A)) > 0
  @test (B ∧ implies(B,A)) == A ∧ (B ∧ implies(B,A))

  @test ¬(A ∨ (¬B)) == ¬(A) ∧ ¬(¬(B))
  @test ¬(A ∨ (¬B)) == ¬(A) ∧ B
  @test A ∧ ¬(¬(A)) == ¬(¬(A))
  @test implies((A∧B), A) == A∨B
end
end