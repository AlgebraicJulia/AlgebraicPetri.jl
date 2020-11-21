using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Edge
import Base.Iterators: flatten
using StatsBase

export Graph

function edgify(δ, transition, reverse::Bool)
  return [Edge(reverse ? ["T_$transition", "S_$k"] :
                         ["S_$k", "T_$transition"],
               Attributes(:label=>"$(δ[k])", :labelfontsize=>"6"))
           for k in collect(keys(δ)) if δ[k] != 0]
end

function Graph(p::AbstractPetriNet)
  statenodes = [Node(string("S_$(sname(p, s))"), Attributes(:shape=>"circle", :color=>"#6C9AC3")) for s in 1:ns(p)]
  transnodes = [Node(string("T_$(tname(p, k))"), Attributes(:shape=>"square", :color=>"#E28F41")) for k in 1:nt(p)]

  graph_attrs = Attributes(:rankdir=>"LR")
  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  stmts = vcat(statenodes, transnodes)

  edges = map(1:nt(p)) do k
    k_name = tname(p, k)
    vcat(edgify(countmap(map(x->sname(p,x), inputs(p, k))), k_name, false),
         edgify(countmap(map(x->sname(p, x), outputs(p, k))), k_name, true))
  end |> flatten |> collect

  stmts = vcat(stmts, edges)
  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
  return g
end

function Graph(op::Union{OpenPetriNet, OpenLabelledPetriNetUntyped, OpenReactionNet, OpenLabelledReactionNetUntyped})
    Graph(apex(op))
end