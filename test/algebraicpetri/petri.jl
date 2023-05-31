module TestPetri

using Test
using AlgebraicPetri
using Catlab.CategoricalAlgebra, Catlab.Theories
using Catlab.Graphs, Catlab.Graphics

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

death_petri = Open(PetriNet(1, 1=>()));
@test to_graphviz(death_petri) isa Graphics.Graphviz.Graph

# Test visualization of nested subgraphs
SIR  = LabelledReactionNet{Float64, Float64}([:S=>1.0, :I=>0.0, :R=>0.0],
                                             (:inf=>0.5)=>((:S,:I)=>(:I,:I)),
                                             (:rec=>0.1)=>(:I=>:R))
stmts1 = Vector{Graphics.Graphviz.Statement}([AlgebraicPetri.Subgraph(SIR; pre="1_")])
graph_attrs = Graphics.Graphviz.Attributes(:rankdir=>"LR")
node_attrs  = Graphics.Graphviz.Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
edge_attrs  = Graphics.Graphviz.Attributes(:splines=>"splines")
g = Graphics.Graphviz.Digraph("G", stmts1; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
stmts2 = Vector{Graphics.Graphviz.Statement}([AlgebraicPetri.tagged_subgraph(g; post="_2")])
g2 = Graphics.Graphviz.Digraph("G", stmts2; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
@test g2 isa Graphics.Graphviz.Graph
@test g2.stmts[1].stmts[1].stmts[1].name == "1_n1_2"

end
