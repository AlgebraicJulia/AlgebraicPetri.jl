using AlgebraicPetri

import Base.:+
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using Catlab.Theories
using Catlab.WiringDiagrams
using Catlab.Graphs.BasicGraphs
import Catlab.Graphs.BasicGraphs: Graph
using Catlab.Graphics.Graphviz: run_graphviz

save_fig(g, fname::AbstractString, format::AbstractString) = begin
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, g, format=format)
    end
end

matching_states(pattern, states) = collect(filter(s->(string(pattern)==first(split(string(s), "_"))), states))

(+)(a::LabelledPetriNet, b::LabelledPetriNet) = begin
  result = copy(a)
  transitions = nparts(a, :T)
  add_parts!(result, :T, nparts(b,:T), tname=subpart(b, :tname))
  add_parts!(result, :I, nparts(b, :I), it=[i + transitions for i in subpart(b, :it)],
                                        is=subpart(b, :is))
  add_parts!(result, :O, nparts(b, :O), ot=[i + transitions for i in subpart(b, :ot)],
                                        os=subpart(b, :os))
  result
end

function connection(conn_petri::LabelledPetriNet, m1::LabelledPetriNet, m2::LabelledPetriNet, ind1, ind2)
  new_conn = copy(conn_petri)
  set_subpart!(new_conn, :sname, vcat(subpart(m1, :sname), subpart(m2, :sname)))
  set_subpart!(new_conn, :tname, [Symbol("$(name)_$(ind1)→$(ind2)") for name in subpart(conn_petri, :tname)])
  Open(new_conn, subpart(m1, :sname)[1:ns(m1)], subpart(m2, :sname)[1:ns(m2)])
end

function cross_uwd(connections::Tuple{LabelledPetriNet,Graph}...)
    rel = RelationDiagram{Symbol}(0)

    # Add populations
    g = first(connections)[2]
    juncs = add_junctions!(rel, nv(g), variable=[Symbol("pop$i") for i in 1:nv(g)])

    # Add epidemiology model boxes
    boxes = add_parts!(rel, :Box, nv(g), name=[Symbol("ep$i") for i in 1:nv(g)])
    add_parts!(rel, :Port, nv(g), junction=juncs, box=boxes)
    
    for conn in 1:length(connections)
      g = connections[conn][2]
      srcs = subpart(g, :src)
      tgts = subpart(g, :tgt)
      # Add cross boxes
      boxes = add_parts!(rel, :Box, ne(g), name=[Symbol("cross_$(conn)_$(srcs[i])_$(tgts[i])") for i in 1:ne(g)])
      add_parts!(rel, :Port, ne(g), junction=srcs, box=boxes)
      add_parts!(rel, :Port, ne(g), junction=tgts, box=boxes)
    end
    rel
end

function index_petri(model::LabelledPetriNet, ind::Int)
  new_petri = copy(model)
  snames = subpart(model, :sname)
  tnames = subpart(model, :tname)
  
  set_subpart!(new_petri, :sname, [Symbol("$(name)_$ind") for name in snames])
  set_subpart!(new_petri, :tname, [Symbol("$(name)_$ind") for name in tnames])
  Open(new_petri, subpart(new_petri, :sname))
end

function diff_petri(model::LabelledPetriNet, diff_states::Array{Symbol})
  states = subpart(model, :sname)
  states2 = [Symbol("$(state)′") for state in states]
  diff = vcat([matching_states(state, states) for state in diff_states]...)

  LabelledPetriNet(vcat(states, states2),
                   [Symbol("diff_$(diff[i])")=>(diff[i]=>Symbol("$(diff[i])′")) for i in 1:length(diff)]...) 
end

function diff_petri(model::LabelledPetriNet)
  states = unique([Symbol(first(split("$name", '_'))) for name in subpart(model, :sname)])
  diff_petri(model, collect(states))
end

function dem_petri(model::LabelledPetriNet, sus_state::Symbol, 
                                            exp_state::Symbol,
                                            inf_states::Array{Symbol})
  states1 = subpart(model, :sname)
  states2 = [Symbol("$(state)′") for state in states1]

  sus1 = matching_states(sus_state, states1)

  # This assumes that the susceptible and corresponding exposed states 
  # are identical after the first '_'
  get_exp(sus) = Symbol(join(vcat("$(exp_state)′",split(string(sus), "_")[2:end]), '_'))
  inf1 = vcat([matching_states(inf, states1) for inf in inf_states]...)

  transitions = vcat([[Symbol("crx_$(sus)′_$(inf)")=>((Symbol("$(sus)′"), inf)=>(inf, get_exp(sus))) for sus in sus1] for inf in inf1]...)

  LabelledPetriNet(vcat(states1, states2),
                   transitions...) 
end

function stratify(model::LabelledPetriNet, connections::Tuple{LabelledPetriNet, Graph}...)
    conn_uwd = cross_uwd(connections...)

    # Calls diffusion_petri for each edge as (src, tgt)
    g = first(connections)[2]
    ep_map = Dict{Symbol, OpenLabelledPetriNet}([Symbol("ep$i")=>index_petri(model, i) for i in 1:nv(g)])

    for conn in 1:length(connections)
      p = connections[conn][1]
      g = connections[conn][2]
      srcs = subpart(g, :src)
      tgts = subpart(g, :tgt)
      for i in 1:ne(g)
        ep_map[Symbol("cross_$(conn)_$(srcs[i])_$(tgts[i])")] =
              connection(p, apex(ep_map[Symbol("ep$(srcs[i])")]),
                            apex(ep_map[Symbol("ep$(tgts[i])")]),
                            srcs[i], tgts[i])
      end
    end

    apex(oapply(conn_uwd, ep_map))
end



test_petri = LabelledPetriNet([:S,:I,:R], :inf=>((:S,:I)=>(:I,:I)),:rec=>(:I=>:R))

clique(n) = begin
  c = Graph(n)
  for i in 1:n
    for j in 1:n
      if i != j
        add_edges!(c, [i],[j])
      end
    end
  end
  c
end

cycle(n) = begin
  c = Graph(n)
  for i in 1:n
    add_edges!(c, [i],[(i)%n+1])
  end
  c
end
