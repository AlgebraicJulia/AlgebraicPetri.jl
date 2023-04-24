using Catlab.CategoricalAlgebra
using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Subgraph
import Base.Iterators: flatten
using StatsBase

export Graph

###########################
# Base Graph Construction #
###########################
def_states(p, s; pos="") = ("s$s", Attributes(:label=>"$(sname(p, s))",
                                              :shape=>"circle",
                                              :color=>"#6C9AC3",
                                              :pos=>pos))
def_trans(p, t; pos="") = ("t$t", Attributes(:label=>"$(tname(p, t))",
                                             :shape=>"square",
                                             :color=>"#E28F41",
                                             :pos=>pos))
def_inpts(p, s, t, i) = (["s$s", "t$t"],Attributes(:labelfontsize=>"6"))
def_otpts(p, s, t, o) = (["t$t", "s$s"],Attributes(:labelfontsize=>"6"))

function Graph(p::AbstractPetriNet; make_states::Function=def_states,
               make_trans::Function=def_trans, make_inpts::Function=def_inpts,
               make_otpts::Function=def_otpts, name="G", prog="dot",
               positions=Dict(:T=>fill("", nt(p)), :S=>fill("", ns(p))), kw...)

  graph_attrs = :graph_attrs ∈ keys(kw) ? Attributes(kw[:graph_attrs]) :
                                          Attributes(:rankdir=>"LR")
  node_attrs  = :node_attrs ∈ keys(kw) ? Attributes(kw[:node_attrs]) :
                                         Attributes(:shape=>"plain",
                                                    :style=>"filled",
                                                    :color=>"white")
  edge_attrs  = :edge_attrs ∈ keys(kw) ? Attributes(kw[:edge_attrs]) :
                                         Attributes(:splines=>"splines")

  statenodes = [Node(make_states(p,s; pos=positions[:S][s])...) for s in 1:ns(p)]
  transnodes = [Node(make_trans(p,t; pos=positions[:T][t])...) for t in 1:nt(p)]
  stmts = vcat(statenodes, transnodes)

  i_edges = map(1:ni(p)) do i
    Edge(make_inpts(p, is(p, i), it(p, i), i)...)
  end |> collect
  o_edges = map(1:no(p)) do o
    Edge(make_otpts(p, os(p, o), ot(p, o), o)...)
  end |> collect
  edges = vcat(i_edges, o_edges)

  # Add count labels if no labels provided
  if all(!(:label ∈ keys(e.attrs)) for e in edges)
    edge_vals = countmap(edges)
    edges = map(filter((v)->v[2] != 0, collect(edge_vals))) do v
      attrs = v[1].attrs
      attrs[:label] = "$(v[2])"
      Edge(v[1].path, attrs)
    end |> collect
  end

  stmts = vcat(stmts, edges)
  g = Graphviz.Digraph(name, stmts; prog=prog,
                                    graph_attrs=graph_attrs,
                                    node_attrs=node_attrs,
                                    edge_attrs=edge_attrs)
  return g
end


####################
# Subgraph Drawing #
####################
add_label(n::Node, pre, post) = Node("$pre$(n.name)$post", n.attrs)
add_label(e::Edge, pre, post) = begin
  path = map(e.path) do n_id
    NodeID("$pre$(n_id.name)$post", n_id.port, n_id.anchor)
  end
  Edge(path, e.attrs)
end
add_label(g::Subgraph, pre, post) = begin
  stmts = map(g.stmts) do st
    add_label(st, pre, post)
  end
  name = "$pre$(g.name)$post"
  Subgraph(name, stmts, g.graph_attrs, g.node_attrs, g.edge_attrs)
end

function tagged_subgraph(g::Graph; pre="", post="")
  stmts = map(g.stmts) do st
    add_label(st, pre, post)
  end
  Subgraph(g.name, stmts, g.graph_attrs, g.node_attrs, g.edge_attrs)
end

function Subgraph(p::AbstractPetriNet; pre="", post="", kw...)
  tagged_subgraph(Graph(p; kw...); pre=pre, post=post)
end

######################
# Specialized Graphs #
######################
function Graph(m::Multispan{<:AbstractPetriNet})
  apex = Subgraph(m.apex; name="clusterApex", pre="ap_")

  leg_graphs = map(enumerate(m.legs)) do (i, l)
    Subgraph(l.codom; name="clusterLeg$i", pre="leg$(i)_")
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

function Graph(p::StructACSetTransformation{<:Any, <:Any, <:AbstractPetriNet, <:AbstractPetriNet})
  Graph(Multispan(p.dom, [p]))
end

function Graph(so::Subobject{<:AbstractPetriNet}; lw = 3.0, kw...)
  p = ob(so)
  maps = hom(so)
  sts = maps.components[:S].func
  trans = maps.components[:T].func
  inps = maps.components[:I].func
  otps = maps.components[:O].func
  sns = sname(p, sts)
  tns = tname(p, trans)
  mkstate(p,s; kw...) = begin
    ds = def_states(p,s; kw...)
    attrs = ds[2]
    attrs[:fillcolor] = "#6C9AC3"
    attrs[:color] = "black"
    attrs[:penwidth] = s ∈ sts ? "$lw" : "0.0"
    (ds[1], attrs)
  end
  mktrans(p,t; kw...) = begin
    dt = def_trans(p,t; kw...)
    attrs = dt[2]
    attrs[:fillcolor] = "#E28F41"
    attrs[:color] = "black"
    attrs[:penwidth] = t ∈ trans ? "$lw" : "0.0"
    (dt[1], attrs)
  end

  mkinpts(p, s, t, i) = begin
    (["s$s", "t$t"],
     Attributes(:labelfontsize=>"6",
                :penwidth=>(i ∈ inps ? "$lw" : "1.0")))
  end

  mkotpts(p, s, t, o) = begin
    (["t$t", "s$s"],
     Attributes(:labelfontsize=>"6",
                :penwidth=>(o ∈ otps ? "$lw" : "1.0")))
  end

  Graph(p; make_states=mkstate, make_trans=mktrans, make_inpts=mkinpts, make_otpts=mkotpts,
        node_attrs  = Attributes(:shape=>"plain", :style=>"filled, solid", :color=>"white"),
        kw...)
end

function Graph(op::Union{OpenPetriNet, OpenLabelledPetriNetUntyped, OpenReactionNet, OpenLabelledReactionNetUntyped})
    Graph(apex(op))
end
