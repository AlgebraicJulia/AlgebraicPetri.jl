using Catlab.CategoricalAlgebra
using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Edge, Subgraph
import Base.Iterators: flatten
using StatsBase

export Graph

function edgify(δ::Dict{Int64, Int64}, transition, reverse::Bool; pre="")
  return [Edge(reverse ? ["\"$(pre)t$transition\"", "\"$(pre)s$k\""] :
                         ["\"$(pre)s$k\"", "\"$(pre)t$transition\""],
               Attributes(:label=>"$(δ[k])", :labelfontsize=>"6"))
           for k in collect(keys(δ)) if δ[k] != 0]
end

function edgify(δ::Dict{Tuple{Int64, Bool}, Int64}, transition, reverse::Bool; pre="")
  return [Edge(reverse ? ["\"$(pre)t$transition\"", "\"$(pre)s$(k[1])\""] :
                         ["\"$(pre)s$(k[1])\"", "\"$(pre)t$transition\""],
               Attributes(:label=>"$(δ[k])", :labelfontsize=>"6",
                          :penwidth=>(k[2] ? "3.0" : "1.0")))
           for k in collect(keys(δ)) if δ[k] != 0]
end

function Graph(p::AbstractPetriNet)
  statenodes = [Node("s$s", Attributes(:label=>"$(sname(p, s))",:shape=>"circle", :color=>"#6C9AC3")) for s in 1:ns(p)]
  transnodes = [Node("t$k", Attributes(:label=>"$(tname(p, k))", :shape=>"square", :color=>"#E28F41")) for k in 1:nt(p)]

  graph_attrs = Attributes(:rankdir=>"LR")
  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  stmts = vcat(statenodes, transnodes)

  edges = map(1:nt(p)) do k
    vcat(edgify(countmap(inputs(p, k)), k, false),
         edgify(countmap(outputs(p, k)), k, true))
  end |> flatten |> collect

  stmts = vcat(stmts, edges)
  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
  return g
end

function Subgraph(p::AbstractPetriNet; label="cluster", pre="")
   statenodes = [Node("\"$(pre)s$s\"",
                      Attributes(:label => "$(sname(p, s))", :shape=>"circle",
                                 :color=>"#6C9AC3")) for s in 1:ns(p)]

  transnodes = [Node("\"$(pre)t$k\"",
                     Attributes(:label => "$(tname(p, k))", :shape=>"square",
                                :color=>"#E28F41")) for k in 1:nt(p)]

  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  stmts = vcat(statenodes, transnodes)

  edges = map(1:nt(p)) do k
    vcat(edgify(countmap(inputs(p, k)), k, false; pre=pre),
         edgify(countmap(outputs(p, k)), k, true; pre=pre))
  end |> flatten |> collect

  stmts = vcat(stmts, edges)
  g = Graphviz.Subgraph(label, stmts; node_attrs=node_attrs, edge_attrs=edge_attrs)
end

function Graph(m::Multispan)
  apex = Subgraph(m.apex; label="clusterApex", pre="ap_")

  leg_graphs = map(enumerate(m.legs)) do (i, l)
    Subgraph(l.codom; label="clusterLeg$i", pre="leg$(i)_")
  end


  graph_attrs = Attributes(:rankdir=>"LR")
  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  t_edges = map(enumerate(m.legs)) do (i, l)
    map(1:nt(m.apex)) do t
      Edge(["\"ap_t$t\"", "\"leg$(i)_t$(l.components[:T](t))\""],
      Attributes(:constraint=>"false", :style=>"dotted"))
    end
  end

  s_edges = map(enumerate(m.legs)) do (i, l)
    map(1:ns(m.apex)) do s
      Edge(["\"ap_s$s\"", "\"leg$(i)_s$(l.components[:S](s))\""],
      Attributes(:constraint=>"false", :style=>"dotted"))
    end
  end
  stmts = vcat([apex], leg_graphs, t_edges..., s_edges...)
  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
end

function Graph(p::ACSetTransformation)
  Graph(Multispan(p.dom, [p]))
end

function Graph(so::Subobject)
  p = ob(so)
  maps = hom(so)
  states = maps.components[:S].func
  trans = maps.components[:T].func
  inps = maps.components[:I].func
  otps = maps.components[:O].func
  sns = sname(p, states)
  tns = tname(p, trans)
  statenodes = [Node("s$s",
                      Attributes(:label => "$(sname(p, s))", :shape=>"circle",
                                 :fillcolor=>"#6C9AC3", :color=>"black",
                                 :penwidth=>(s ∈ states ? "3.0" : "0.0" ))) for s in 1:ns(p)]

  transnodes = [Node("t$k",
                     Attributes(:label => "$(tname(p, k))", :shape=>"square",
                                :fillcolor=>"#E28F41", :color=>"black",
                                :penwidth=>(k ∈ trans ? "3.0" : "0.0" ))) for k in 1:nt(p)]

  graph_attrs = Attributes(:rankdir=>"LR")
  node_attrs  = Attributes(:shape=>"plain", :style=>"filled, solid", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  stmts = vcat(statenodes, transnodes)

  edges = map(1:nt(p)) do k
    vcat(edgify(countmap(map(x->(p[x, :is], x ∈ inps), incident(p, k, :it))), k, false),
         edgify(countmap(map(x->(p[x, :os], x ∈ otps), incident(p, k, :ot))), k, true))
  end |> flatten |> collect

  stmts = vcat(stmts, edges)
  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
end

function Graph(op::Union{OpenPetriNet, OpenLabelledPetriNetUntyped, OpenReactionNet, OpenLabelledReactionNetUntyped})
    Graph(apex(op))
end
