@testset "PetriNet Visualization" begin
  sir_petri = PetriNet(3, ((1, 2), (2, 2)), (2, 3))
  sir_lpetri = LabelledPetriNet([:S, :I, :R], :inf=>((:S, :I), (:I, :I)), :rec=>(:I, :R))
  β(u,t) = 1 / sum(u)
  γ = .25
  sir_rxn = ReactionNet{Function, Int}([990, 10, 0], (β)=>((1, 2)=>(2,2)), (t->γ)=>(2=>3))
  open_sir_rxn = Open([1,2], sir_rxn, [3])
  sir_lrxn = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .001)=>((:S, :I)=>(:I,:I)), (:rec, .25)=>(:I=>:R))
  open_sir_lrxn = Open([:S,:I], sir_lrxn, [:R])

  sir_tpetri= PetriNet(TransitionMatrices(sir_petri))

  @test typeof(Graph(sir_petri)) == Graph
  @test typeof(Graph(sir_lpetri)) == Graph
  @test typeof(Graph(sir_rxn)) == Graph
  @test typeof(Graph(open_sir_rxn)) == Graph
  @test typeof(Graph(sir_lrxn)) == Graph
  @test typeof(Graph(open_sir_lrxn)) == Graph
end

@testset "Presentation Visualization" begin
  @present TheoryPetriNet(FreeSymmetricMonoidalCategory) begin
    T::Ob
    S::Ob
    I::Ob
    O::Ob

    it::Hom(I,T)
    is::Hom(I,S)
    ot::Hom(O,T)
    os::Hom(O,S)
  end

  @present TheoryMLSchema(FreeSymmetricMonoidalCategory) begin
    Files::Ob
    Images::Ob
    NeuralNet::Ob
    Accuracy::Ob
    Metadata::Ob

    extract::Hom(Files, Images)
    split::Hom(Images, Images⊗Images)
    train::Hom(NeuralNet⊗Images, NeuralNet⊗Metadata)
    evaluate::Hom(NeuralNet⊗Images, Accuracy⊗Metadata)
  end

  @test typeof(Graph(TheoryPetriNet)) == Graph
  @test typeof(Graph(TheoryMLSchema)) == Graph
end
