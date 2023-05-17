using Catlab.CategoricalAlgebra, Catlab.Graphics.Graphviz
using Catlab.Graphs: PropertyGraph

import Catlab.Graphics.Graphviz: Subgraph
import Catlab.Graphics: to_graphviz, to_graphviz_property_graph

const GRAPH_ATTRS = Dict(:rankdir=>"LR")
const NODE_ATTRS = Dict(:shape => "plain", :style=>"filled")
const EDGE_ATTRS = Dict(:splines=>"splines")


# Single Petri Nets
###################

to_graphviz(pn::Union{
              AbstractPetriNet,
              Subobject{<:AbstractPetriNet},
              StructuredMulticospan{<:StructuredCospans.AbstractDiscreteACSet{<:AbstractPetriNet}},
            }; kw...) =
  to_graphviz(to_graphviz_property_graph(pn; kw...))

function to_graphviz_property_graph(pn::AbstractPetriNet;
    prog::AbstractString="dot", graph_attrs::AbstractDict=Dict(),
    node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(), name::AbstractString="G", kw...)
  pg = PropertyGraph{Any}(; name = name, prog = prog,
    graph = merge!(GRAPH_ATTRS, graph_attrs),
    node = merge!(NODE_ATTRS, node_attrs),
    edge = merge!(EDGE_ATTRS, edge_attrs),
  )
  S_vtx = Dict(map(parts(pn, :S)) do s
    s => add_vertex!(pg; label="$(sname(pn, s))", shape="circle", color="#6C9AC3")
  end)
  T_vtx = Dict(map(parts(pn, :T)) do t
    t => add_vertex!(pg; label="$(tname(pn, t))", shape="square", color="#E28F41")
  end)

  edges = Dict{Tuple{Int,Int}, Int}()
  map(parts(pn, :I)) do i
    edge = (S_vtx[pn[i, :is]], T_vtx[pn[i, :it]])
    edges[edge] = get(edges, edge, 0) + 1
  end
  map(parts(pn, :O)) do o
    edge = (T_vtx[pn[o, :ot]], S_vtx[pn[o, :os]])
    edges[edge] = get(edges, edge, 0) + 1
  end
  for ((src, tgt),count) in edges
    add_edge!(pg, src, tgt, label="$(count)")
  end

  pg
end

function to_graphviz_property_graph(so::Subobject{<:AbstractPetriNet};
    prog::AbstractString="dot", graph_attrs::AbstractDict=Dict(),
    node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
    name::AbstractString="G", lw::Number=3.0, kw...)
  pn = ob(so)
  comps = hom(so).components
  sts = comps[:S].func
  trans = comps[:T].func
  inps = comps[:I].func
  otps = comps[:O].func

  pg = PropertyGraph{Any}(; name = name, prog = prog,
    graph = merge!(GRAPH_ATTRS, graph_attrs),
    node = merge!(NODE_ATTRS, node_attrs),
    edge = merge!(EDGE_ATTRS, edge_attrs),
  )
  S_vtx = Dict(map(parts(pn, :S)) do s
    s => add_vertex!(pg;
      label="$(sname(pn, s))",
      shape="circle",
      fillcolor="#6C9AC3",
      penwidth = s ∈ sts ? "$lw" : "0.0"
    )
  end)
  T_vtx = Dict(map(parts(pn, :T)) do t
    t => add_vertex!(pg;
      label="$(tname(pn, t))",
      shape="square",
      fillcolor="#E28F41",
      penwidth = t ∈ trans ? "$lw" : "0.0"
    )
  end)

  edges = Dict{Tuple{Int,Int}, Tuple{Int,Number}}()
  map(parts(pn, :I)) do i
    edge = (S_vtx[pn[i, :is]], T_vtx[pn[i, :it]])
    (count, weight) = get(edges, edge, (0,1))
    edges[edge] = (count + 1, i ∈ inps ? lw : weight)
  end
  map(parts(pn, :O)) do o
    edge = (T_vtx[pn[o, :ot]], S_vtx[pn[o, :os]])
    (count, weight) = get(edges, edge, (0,1))
    edges[edge] = (count + 1, o ∈ otps ? lw : weight)
  end
  for ((src, tgt),(count,weight)) in edges
    add_edge!(pg, src, tgt, label="$(count)", penwidth="$(weight)")
  end

  pg
end

to_graphviz_property_graph(op::StructuredMulticospan{<:StructuredCospans.AbstractDiscreteACSet{<:AbstractPetriNet}}; kw...) =
  to_graphviz_property_graph(apex(op); kw...)

# Petri Nets Multi-spans
########################

function to_graphviz(m::Multispan{<:AbstractPetriNet};
    prog::AbstractString="dot", graph_attrs::AbstractDict=Dict(),
    node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
    name::AbstractString="G", kw...)
  apex = Subgraph(m.apex; name="clusterApex", pre="ap_", prog=prog, graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)

  leg_graphs = map(enumerate(m.legs)) do (i, l)
    Subgraph(l.codom; name="clusterLeg$i", pre="leg$(i)_", prog=prog, graph_attrs=graph_attrs, node_attrs=node_attrs, edge_attrs=edge_attrs)
  end

  s_edges = map(enumerate(m.legs)) do (i, l)
    map(1:ns(m.apex)) do s
      Edge(["\"ap_n$s\"", "\"leg$(i)_n$(l.components[:S](s))\""],
      Attributes(:constraint=>"false", :style=>"dotted"))
    end
  end

  t_edges = map(enumerate(m.legs)) do (i, l)
    map(1:nt(m.apex)) do t
      Edge(["\"ap_n$(t+ns(m.apex))\"", "\"leg$(i)_n$(l.components[:T](t)+ns(l.codom))\""],
      Attributes(:constraint=>"false", :style=>"dotted"))
    end
  end

  g = Digraph(name, vcat([apex], leg_graphs, t_edges..., s_edges...); prog=prog,
    graph_attrs = merge!(GRAPH_ATTRS, graph_attrs),
    node_attrs = merge!(NODE_ATTRS, node_attrs),
    edge_attrs = merge!(EDGE_ATTRS, edge_attrs)
  )
end

to_graphviz(p::StructACSetTransformation{<:Any, <:Any, <:AbstractPetriNet, <:AbstractPetriNet}; kw...) =
  to_graphviz(Multispan(p.dom, [p]); kw...)

# Subgraph Extensions
#####################

add_label(n::Node, pre, post) = Node("$pre$(n.name)$post", n.attrs)
add_label(e::Edge, pre, post) = Edge(map(e.path) do n_id
                                  NodeID("$pre$(n_id.name)$post", n_id.port, n_id.anchor)
                                end, e.attrs)
add_label(g::Subgraph, pre, post) = Subgraph("$pre$(g.name)$post",
                                             map(g.stmts) do st
                                               add_label(st, pre, post)
                                             end,
                                             g.graph_attrs, g.node_attrs, g.edge_attrs)

tagged_subgraph(g::Graph; pre="", post="") = Subgraph(g.name,
                                                      map(g.stmts) do st
                                                        add_label(st, pre, post)
                                                      end,
                                                      g.graph_attrs, g.node_attrs, g.edge_attrs)

Subgraph(p::AbstractPetriNet; pre="", post="", kw...) = tagged_subgraph(to_graphviz(p; kw...); pre=pre, post=post)
