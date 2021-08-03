SIRD = LabelledReactionNet{Float64, Float64}([:S=>0.0, :I=>0.0, :R=>0.0, :D=>0.0],
                                             (:inf=>0.0)=>((:S,:I)=>(:I,:I)),
                                             (:rec=>0.0)=>(:I=>:R),
												                     (:death=>0.0)=>(:I=>:D))

SIR  = LabelledReactionNet{Float64, Float64}([:S=>1.0, :I=>0.0, :R=>0.0],
                                             (:inf=>0.5)=>((:S,:I)=>(:I,:I)),
                                             (:rec=>0.1)=>(:I=>:R))
models = [SIR, SIRD]

for pn in [models,
           ReactionNet{Float64, Float64}.([SIR, SIRD]),
           LabelledPetriNet.(models),
           PetriNet.(models)]
  c_res = compare(pn[1], pn[2])
  @test apex(c_res) == pn[1]
  @test length(legs(c_res)) == 2

  @test Graph(legs(c_res)[1]) isa Graph
  @test Graph(c_res) isa Graph

  if Catlab.VERSION >= v"0.12.7"
    so = Subobject.(legs(c_res))
    for s in so
      @test Graph(s) isa Graph
      @test codom(hom(s)) == pn[2]
    end

    @test dom(hom(~foldl(∨, so) ∨ foldl(∨, so))) == pn[2]
  end
end
