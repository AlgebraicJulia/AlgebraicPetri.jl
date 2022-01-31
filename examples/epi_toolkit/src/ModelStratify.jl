module ModelStratify

# The available functions, all have necessary docstrings
export dem_strat, diff_strat, diff_petri, dem_petri, stratify,
       serialize, deserialize, save_petri, save_json, save_model,
       ScaleGraph
using AlgebraicPetri
using Semagrams

using JSON

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs.BasicGraphs
using Catlab.Graphs
using Catlab.Graphs.BasicGraphs: Graph, TheoryGraph
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics.Graphviz: run_graphviz
using Catlab.Graphics.GraphvizGraphs
using AlgebraicPetri: scale_locs
import AlgebraicPetri: locations

import Base.:+
import Base.convert
Base.convert(::Type{Symbol}, str::String) = Symbol(str)

function Base.setindex!(pn::LabelledReactionNet, val, label)
  if label ∈ tnames(pn)
    pn[first(incident(pn, label, :tname)), :rate] = val
  elseif label ∈ snames(pn)
    pn[first(incident(pn, label, :sname)), :concentration] = val
  else
    error("$label is not a valid state or transition name")
  end
end
function Base.getindex(pn::LabelledReactionNet, label)
  if label ∈ tnames(pn)
    pn[first(incident(pn, label, :tname)), :rate]
  elseif label ∈ snames(pn)
    pn[first(incident(pn, label, :sname)), :concentration]
  else
    error("$label is not a valid state or transition name")
  end
end
Base.setindex!(pn::LabelledReactionNet, val, label::Symbol, city) =
  setindex!(pn, val, Symbol(label, "@", city))

Base.getindex(pn::LabelledReactionNet, label::Symbol, city) =
  getindex(pn, Symbol(label, "@", city))


# Define helper functions for defining the two types of
# reactions in an epidemiology Model. Either a state
# spontaneously changes, or one state causes another to change

##################
# Weighted Graph #
##################

@present TheoryScaleGraph <: TheoryGraph begin
  Scale::Data
  Label::Data
  crossexposure::Attr(E, Scale)
  pop_frac::Attr(V, Scale)
  exposure::Attr(V, Scale)
  node_label::Attr(V, Label)
end

const AbstractScaleGraph = AbstractACSetType(TheoryScaleGraph)
const ScaleGraph = ACSetType(TheoryScaleGraph, index=[:src,:tgt])

iscale(s::ScaleGraph, i::Int) = s[i, :exposure]
cscale(s::ScaleGraph, i::Int) = s[i, :pop_frac]
escale(s::ScaleGraph, i::Int) = s[i, :crossexposure]
label(s::ScaleGraph, i::Int) = s[i, :node_label]
infection_scale(s::ScaleGraph) =
  1/sum(vcat([cscale(s, src(s, e))*cscale(s, tgt(s,e))*escale(s,e)
               for e in 1:ne(s)],
             [iscale(s,v) * cscale(s,v)^2 for v in 1:nv(s)]...))

labels(s::ScaleGraph) = s[:node_label]

matching_states(pattern, states) = collect(filter(s->(string(pattern)==first(split(string(s), "@"))), states))

(+)(a::LabelledPetriNet, b::LabelledPetriNet) = begin
  result = copy(a)

  b_tran = add_parts!(result, :T, nparts(b,:T), tname=subpart(b, :tname))
  add_parts!(result, :I, nparts(b, :I), it=b_part[subpart(b,:it)],
                                        is=subpart(b, :is))
  add_parts!(result, :O, nparts(b, :O), ot=b_part[subpart(b, :ot)],
                                        os=subpart(b, :os))
  result
end


# This function preserves the concentation of the states from the first argument
(+)(a::LabelledReactionNet, b::LabelledReactionNet) = begin
  result = copy(a)

  b_tran = add_parts!(result, :T, nparts(b,:T), tname=subpart(b, :tname), rate=subpart(b, :rate))
  add_parts!(result, :I, nparts(b, :I), it=b_part[subpart(b,:it)],
                                        is=subpart(b, :is))
  add_parts!(result, :O, nparts(b, :O), ot=b_part[subpart(b, :ot)],
                                        os=subpart(b, :os))
  result
end


# Given a base connection net along with two networks to connect, this will create a special
# connection petrinet which is consistent with the attributes on the networks to connect
function connection(conn_petri::LabelledPetriNet, m1::LabelledPetriNet, m2::LabelledPetriNet, ind1, ind2)
  new_conn = copy(conn_petri)
  set_subpart!(new_conn, :sname, vcat(subpart(m1, :sname), subpart(m2, :sname)))
  set_subpart!(new_conn, :tname, [Symbol("$(name)_$(ind1)_$(ind2)") for name in subpart(conn_petri, :tname)])
  Open(new_conn, subpart(m1, :sname)[1:ns(m1)], subpart(m2, :sname)[1:ns(m2)])
end

function connection(conn_petri::LabelledReactionNet, m1::LabelledReactionNet,
                    m2::LabelledReactionNet, ind1, ind2, rate_fact, cscale1, cscale2)
  new_conn = copy(conn_petri)
  set_subpart!(new_conn, :sname, vcat(subpart(m1, :sname), subpart(m2, :sname)))
  set_subpart!(new_conn, :concentration, vcat(subpart(m1, :concentration), subpart(m2, :concentration)))
  set_subpart!(new_conn, :tname, [Symbol("$(name)_$(ind1)_$(ind2)") for name in subpart(conn_petri, :tname)])

  set_subpart!(new_conn, :rate, subpart(conn_petri,:rate) .* rate_fact)
  Open(new_conn, subpart(m1, :sname)[1:ns(m1)], subpart(m2, :sname)[1:ns(m2)])
end

# Creates a UWD which provides the composition pattern of the stratification
function cross_uwd(connections::Tuple{Union{LabelledPetriNet, LabelledReactionNet},
                                      Union{Graph, ScaleGraph}}...)
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

# Creates a copy of a petrinet with indices added to the names to create
# uniquely named states
function index_petri(model::LabelledPetriNet, ind::Int)
  new_petri = copy(model)
  snames = subpart(model, :sname)
  tnames = subpart(model, :tname)

  set_subpart!(new_petri, :sname, [Symbol("$(name)@$ind") for name in snames])
  set_subpart!(new_petri, :tname, [Symbol("$(name)@$ind") for name in tnames])
  Open(new_petri, subpart(new_petri, :sname))
end

# Creates a copy of a petrinet with indices added to the names to create
# uniquely named states
function index_petri(model::LabelledReactionNet, ind, rscale::Number, cscale::Number; scaled_transitions::Union{Nothing, Array{Symbol}}=nothing)
  new_petri = copy(model)
  snames = subpart(model, :sname)
  tnames = subpart(model, :tname)

  rscales = [(isnothing(scaled_transitions) ||
              new_petri[i, :tname] ∈ scaled_transitions) ? rscale : 1
              for i in 1:nparts(new_petri, :T)]

  set_subpart!(new_petri, :sname, [Symbol("$(name)@$ind") for name in snames])
  set_subpart!(new_petri, :tname, [Symbol("$(name)@$ind") for name in tnames])

  set_subpart!(new_petri, :rate, rscales .* subpart(new_petri, :rate))
  set_subpart!(new_petri, :concentration, cscale*subpart(new_petri, :concentration))
  Open(new_petri, subpart(new_petri, :sname))
end

# Creates a diffusion petrinet between two petrinets of the shape `model`
function diff_petri(model::LabelledPetriNet, diff_states::Array{Symbol})
  states = subpart(model, :sname)
  states2 = [Symbol("$(state)′") for state in states]
  diff = vcat([matching_states(state, states) for state in diff_states]...)

  LabelledPetriNet(vcat(states, states2),
                   [Symbol("diff_$(diff[i])")=>(diff[i]=>Symbol("$(diff[i])′")) for i in 1:length(diff)]...)
end

# Creates a diffusion petrinet between two petrinets of the shape `model`
function diff_petri(model::LabelledReactionNet{R,C}, diff_states::Array{<:Pair{Symbol, <:Number}}) where {R,C}
  states = subpart(model, :sname)
  concs  = concentrations(model)
  states2 = [Symbol("$(state)′") for state in states]
  diff = vcat([[(m_state, state[2]) for m_state in matching_states(state[1], states)] for state in diff_states]...)

  new_states = vcat([state=>concs[state] for state in states],
                    [Symbol("$(state)′")=>concs[state] for state in states])
  transitions = [(Symbol("diff_$(t[1])")=>t[2],(t[1]=>Symbol("$(t[1])′"))) for t in diff]

  LabelledReactionNet{R, C}(new_states,
                                      transitions...)
end

function diff_petri(model::LabelledPetriNet)
  states = unique([Symbol(first(split("$name", '@'))) for name in subpart(model, :sname)])
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
  get_exp(sus) = Symbol(join(vcat("$(exp_state)′",split(string(sus), "@")[2:end]), '@'))
  inf1 = vcat([matching_states(inf, states1) for inf in inf_states]...)

  transitions = vcat([[Symbol("crx_$(sus)′_$(inf)")=>((Symbol("$(sus)′"), inf)=>(inf, get_exp(sus))) for sus in sus1] for inf in inf1]...)

  LabelledPetriNet(vcat(states1, states2),
                   transitions...)
end

function dem_petri(model::LabelledReactionNet{R,C}, transitions::Vector) where {R, C}
  states1 = subpart(model, :sname)
  concs = concentrations(model)
  transitions = [(Symbol("crx_$(sus)′_$(inf)")=>rate)=>((Symbol(sus, "′"), inf)=>(Symbol(exp, "′"), inf))
                 for (((sus, inf), exp), rate) in transitions]

  states = vcat([state=>concs[state] for state in states1],
                [Symbol("$(state)′")=>concs[state] for state in states1])
  LabelledReactionNet{R, C}(states, transitions...)

end

function dem_petri(model::LabelledReactionNet{R,C}, sus_state::Symbol,
                                               exp_state::Symbol,
                                               inf_states::Array{<:Pair{Symbol, <:Number}}) where {R,C}
  states1 = subpart(model, :sname)
  concs = concentrations(model)

  sus1 = matching_states(sus_state, states1)

  # This assumes that the susceptible and corresponding exposed states
  # are identical after the first '_'
  get_exp(sus) = Symbol(join(vcat("$(exp_state)",split(string(sus), "@")[2:end]), '@'))
  inf1 = vcat([[(m_state, state[2]) for m_state in matching_states(state[1], states1)] for state in inf_states]...)

  transitions = vcat([[((sus, inf[1])=>get_exp(sus))=>inf[2]
                       for sus in sus1] for inf in inf1]...)
  dem_petri(model, transitions)
end

infection(model, t) = begin
  inp, otp = (inputs(model, t), outputs(model, t))
  (length(inp) == 2 && length(otp) == 2) || return nothing
  inf_ind = findfirst(op -> op ∈ inp, otp)
  inf = otp[inf_ind]
  isnothing(inf) && return nothing
  exp_ind = 3 - inf_ind
  exp = otp[exp_ind]
  sus = inp[3 - findfirst(ip -> ip == inf, inp)]
  return ((sname(model, sus), sname(model, inf))=>sname(model, exp))=>rate(model, t)
end

infections(model) = [infection(model, t) for t in 1:nt(model)]

function dem_petri(model::LabelledReactionNet{R,C}) where {R,C}
  transitions = infections(model)
  filter!(x -> !isnothing(x), transitions)
  dem_petri(model, transitions)
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

function stratify(model::LabelledReactionNet, connections::Tuple{LabelledReactionNet, ScaleGraph}...; scaled_transitions::Union{Nothing, Array{Symbol}}=nothing)
    conn_uwd = cross_uwd(connections...)

    # Calls diffusion_petri for each edge as (src, tgt)
    g = first(connections)[2]
    inf_scale = infection_scale(g)
    scaled_transitions = tnames(model)[.!(isnothing.(infections(model)))]
    ep_map = Dict{Symbol, OpenLabelledReactionNet}(
                 [Symbol("ep$i") =>
                  index_petri(model, label(g,i), inf_scale * iscale(g,i), cscale(g,i);
                              scaled_transitions=scaled_transitions) for i in 1:nv(g)])

    for conn in 1:length(connections)
      p = connections[conn][1]
      g = connections[conn][2]
      srcs = subpart(g, :src)
      tgts = subpart(g, :tgt)
      for i in 1:ne(g)
        ep_map[Symbol("cross_$(conn)_$(srcs[i])_$(tgts[i])")] =
              connection(p, apex(ep_map[Symbol("ep$(srcs[i])")]),
                            apex(ep_map[Symbol("ep$(tgts[i])")]),
                            srcs[i], tgts[i], inf_scale*escale(g,i),
                            cscale(g, srcs[i]), cscale(g, tgts[i]))
      end
    end

    apex(oapply(conn_uwd, ep_map))
end

""" diff_strat(epi_model, connection_graph)
  This function takes in a LabelledPetriNet and a graph which describes
  geographical connections. It returns a LabelledPetriNet which models
  diffusion between geographic populations described by the given graph.
"""
function diff_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph)
  diff_conn = diff_petri(epi_model)
  stratify(epi_model, (diff_conn, connection_graph))
end

function diff_strat(epi_model::LabelledReactionNet, connection_graph::ScaleGraph,
                    labels::Array{<:Pair{Symbol, <:Number}}; kw...)
  diff_conn = diff_petri(epi_model, labels)
  stratify(epi_model, (diff_conn, connection_graph); kw...)
end

function diff_strat(epi_model::LabelledReactionNet{R,C}, connection_graph::ScaleGraph,
                    labels::Array{Symbol}; kw...) where {R,C}
  new_labels = [l=>zero(C) for l in labels]
  diff_conn = diff_petri(epi_model, new_labels)
  stratify(epi_model, (diff_conn, connection_graph); kw...)
end

function diff_strat(epi_model::LabelledReactionNet{R,C}, connection_graph::ScaleGraph; kw...) where {R,C}
  labels = [l=>zero(C) for l in snames(epi_model)]
  diff_conn = diff_petri(epi_model, labels)
  stratify(epi_model, (diff_conn, connection_graph); kw...)
end


""" dem_strat(epi_model, connection_graph, sus_state, exp_state, inf_states)
  This function takes in a LabelledPetriNet and a graph which describes
  infection connections between populations. It also takes in the symbol used
  for susceptible states, the symbol used for the exposed state, and an array
  of symbols for states that can expose susceptible states. It returns a
  LabelledPetriNet which models diffusion between populations described by the
  given graph.
"""
function dem_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph,
                   sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol})
  dem_conn = dem_petri(epi_model, sus_state, exp_state, inf_states)
  stratify(epi_model, (dem_conn, connection_graph))
end

function dem_strat(epi_model::LabelledReactionNet{R,C}, connection_graph::ScaleGraph,
                   sus_state::Symbol, exp_state::Symbol, inf_states::Array{<:Pair{Symbol, <:Number}}; kw...) where {R,C}
  dem_conn = dem_petri(epi_model, sus_state, exp_state, inf_states)
  stratify(epi_model, (dem_conn, connection_graph); kw...)
end

function dem_strat(epi_model::LabelledReactionNet{R,C}, connection_graph::ScaleGraph,
                   sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol}; kw...) where {R,C}
  new_inf_states = [l=>zero(C) for l in inf_states]
  dem_conn = dem_petri(epi_model, sus_state, exp_state, new_inf_states)
  stratify(epi_model, (dem_conn, connection_graph); kw...)
end

function dem_strat(epi_model::LabelledReactionNet{R,C}, connection_graph::ScaleGraph; kw...) where {R,C}
  dem_conn = dem_petri(epi_model)
  stratify(epi_model, (dem_conn, connection_graph); kw...)
end

""" Serialize an ACSet object to a JSON string
"""
serialize(x::ACSet; io=stdout) = JSON.print(io, x.tables, 2)

serialize_string(x::ACSet) = JSON.json(x.tables, 2)

""" Deserialize a dictionary from a parsed JSON string to an object of the given ACSet type
"""
function deserialize(input::Dict, type)
    out = type()
    for (k,v) ∈ input
        add_parts!(out, Symbol(k), length(v))
    end
    for l ∈ values(input)
        for (i, j) ∈ enumerate(l)
            for (k,v) ∈ j
                set_subpart!(out, i, Symbol(k), v)
            end
        end
    end
    out
end

""" Deserialize a JSON string to an object of the given ACSet type
"""
deserialize(input::String, type) = deserialize(JSON.parse(input), type)

""" Save Petri net as an svg image
"""
save_petri(g, fname::AbstractString, format::AbstractString) =
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, AlgebraicPetri.Graph(g), format=format)
    end

""" Save Graph as an svg image
"""
save_graph(g, fname::AbstractString, format::AbstractString) =
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, g, format=format)
    end

""" Show Graph
"""
show_graph(sg::ScaleGraph) = begin
	g = BasicGraphs.Graph()
	copy_parts!(g, sg)
  vprops(i) = Dict(:label=>"$(sg[i, :node_label])\n pop: $(sg[i, :pop_frac])\n exp: $(sg[i, :exposure])")
	eprops(i) = Dict(:label=>"$(sg[i, :crossexposure])")
	p = PropertyGraph{Any}(g, vprops, eprops; prog="dot",
	                            graph=Dict(:rankdir => "LR"),
	                            node = GraphvizGraphs.default_node_attrs(true),
	                            edge = Dict(:arrowsize => "0.5"))
	to_graphviz(p)
end

show_graph(g) = to_graphviz(g, node_labels=true)

""" Save serialization to json file
"""
save_json(C, fname::AbstractString) =
    open(string(fname, ".json"),"w") do f
        serialize(C; io=f)
    end

""" Save both Petri graph as svg and json file
"""
function save_model(model, fname::AbstractString)
    save_json(model, fname)
    save_petri(model, fname, "svg");
end

@semagramschema ScaleGraphSchema(TheoryScaleGraph) begin
  @box V Circle label=:node_label
  @wire E(src,tgt)
  @data Scale Numeric
  @data Label Stringlike
end;

function locations(sg::ScaleGraph, lsd::LocatedSemagramData)
  blocs = lsd.boxlocs
  sclocs = scale_locs(collect(values(blocs)); scale = 1.0)
  for (i,k) in enumerate(keys(blocs))
    blocs[k] = sclocs[i]
  end

  boxes = sort(collect(lsd.sg.boxes); by=first)
  v_boxes = filter(b -> b[2].ty == :V, boxes)

  v_label2id = Dict(b[2].weights[:node_label]=>b[1] for b in v_boxes)

  Dict(:V=>[blocs[v_label2id["$v"]] for v in labels(sg)])
end

function strat_location(strat_model, pn, sg, pnlocs, sglocs; scale = 1.5)

	all_pnlocs = vcat(pnlocs[:T], pnlocs[:S])
  w, h = (maximum(first.(all_pnlocs)) .- minimum(first.(all_pnlocs)) + 10.0,
          maximum(last.(all_pnlocs)) .- minimum(last.(all_pnlocs)) + 10.0)

  # Calculate scale for scaled-graph coordinates
  sc_dists = [abs.(sglocs[:V][i] .- sglocs[:V][j]) ./ (w,h) for i in 1:length(sglocs[:V]), j in 1:length(sglocs[:V]) if i != j]
  ratio = 1/minimum(maximum.(sc_dists)) * scale

  offset_dists = map(p -> ratio .* p, sglocs[:V])

  tname2loc = Dict("$t"=>pnlocs[:T][i] for (i,t) in enumerate(tnames(pn)))
  sname2loc = Dict("$s"=>pnlocs[:S][i] for (i,s) in enumerate(snames(pn)))
  vname2loc = Dict("$v"=>offset_dists[i] for (i,v) in enumerate(labels(sg)))

  tlocs = map(tnames(strat_model)) do t
    if length(split("$t", "@")) == 2
			p, c = split("$t", "@")
			tname2loc[p] .+ vname2loc[c]
    else
      nothing
    end
  end
  slocs = map(snames(strat_model)) do s
    if length(split("$s", "@")) == 2
      p, c = split("$s", "@")
      sname2loc[p] .+ vname2loc[c]
    else
      nothing
    end
  end
  Dict(:T=>tlocs, :S=>slocs)
end

end
