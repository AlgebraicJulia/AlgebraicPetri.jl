module BilayerNetworks

using AlgebraicPetri
using Catlab
import Catlab.CategoricalAlgebra: migrate!

export ThBilayerNetwork, AbstractBilayerNetwork, BilayerNetwork,
    ThLabelledBilayerNetwork, AbstractLabelledBilayerNetwork, LabelledBilayerNetwork

@present ThBilayerNetwork(FreeSchema) begin
    (Qin, Qout, Win, Wn, Wa, Box)::Ob
    arg::Hom(Win, Qin)
    call::Hom(Win, Box)
    influx::Hom(Wa, Box)
    infusion::Hom(Wa, Qout)
    efflux::Hom(Wn, Box)
    effusion::Hom(Wn, Qout)
end

@present ThLabelledBilayerNetwork <: ThBilayerNetwork begin
    Name::AttrType

    parameter::Attr(Box, Name)
    variable::Attr(Qin, Name)
    tanvar::Attr(Qout, Name)
end

@abstract_acset_type AbstractBilayerNetwork
@acset_type BilayerNetwork(ThBilayerNetwork) <: AbstractBilayerNetwork

@abstract_acset_type AbstractLabelledBilayerNetwork <: AbstractBilayerNetwork
@acset_type LabelledBilayerNetworkUntyped(ThLabelledBilayerNetwork) <: AbstractLabelledBilayerNetwork
const LabelledBilayerNetwork = LabelledBilayerNetworkUntyped{Symbol}

function balance!(bn::AbstractBilayerNetwork)
  for t in 1:nparts(bn, :Box)
    for s in 1:nparts(bn, :Qin)
      n_in = count(i->(bn[i, :arg] == s), incident(bn, t, :call))
      n_neg = count(i->(bn[i, :effusion] == s), incident(bn, t, :efflux))
      n_pos = count(i->(bn[i, :infusion] == s), incident(bn, t, :influx))
      n_diff = n_pos - n_neg

      req_pairs = n_neg - n_in

      req_pairs <= n_pos ||
        error("Mass action does not allow species $(bn[s, :variable]) to be "*
              "removed without contributing to the rate.")

      if req_pairs < 0
        add_parts!(bn, :Wn, -req_pairs, effusion=s, efflux=t)
        add_parts!(bn, :Wa, -req_pairs, infusion=s, influx=t)
      elseif req_pairs > 0
        neg = filter(i -> s == bn[i, :effusion], incident(bn, t, :efflux)) |> sort
        pos = filter(i -> s == bn[i, :infusion], incident(bn, t, :influx)) |> sort
        rem_parts!(bn, :Wn, neg[1:req_pairs])
        rem_parts!(bn, :Wa, pos[1:req_pairs])
      end
    end
  end
  bn
end

function migrate!(bn::AbstractBilayerNetwork, pn::AbstractPetriNet)
    T,S,I,O, it, is, ot, os = generators(SchPetriNet)
    (Qin, Qout, Win, Wn, Wa, Box, arg, call, influx, infusion, efflux, effusion
     ) = generators(ThBilayerNetwork)
    migrate!(bn, pn,
             Dict(Qin=>S, Qout=>S, Box=>T, Win=>I, Wn=>I, Wa=>O),
             Dict(arg=>is,
                  call=>it,
                  efflux=>it,
                  effusion=>is,
                  influx=>ot,
                  infusion=>os))
end

function migrate!(bn::AbstractLabelledBilayerNetwork, pn::AbstractPetriNet)
    S = Presentation(acset_schema(pn)) # could be SchPetri or SchLabeledPetri
    (Qin, Qout, Win, Wn, Wa, Box, arg, call, influx, infusion, efflux, effusion,
     Name, parameter, variable, tanvar) = generators(ThLabelledBilayerNetwork)

    gens = Dict{Any,Any}(Qin=>S[:S], Qout=>S[:S], Box=>S[:T], Win=>S[:I], Wn=>S[:I], Wa=>S[:O])
    parts = Dict{Any,Any}(arg=>S[:is],
                 call=>S[:it],
                 efflux=>S[:it],
                 effusion=>S[:is],
                 influx=>S[:ot],
                 infusion=>S[:os])
    if :Name in attrtypes(acset_schema(pn))
        gens[Name] = S[:Name]
    end
    if has_subpart(pn, :sname)
        parts[variable] = S[:sname]
        parts[tanvar] = S[:sname]
    end
    if has_subpart(pn, :tname)
        parts[parameter] = S[:tname]
    end
    migrate!(bn, pn, gens, parts)
end


function migrate!(pn::AbstractPetriNet, bn::AbstractBilayerNetwork)
    P, N = Presentation.(acset_schema.([pn, bn])) 
    bnc = copy(bn)
    balance!(bnc)
    migrate!(pn,bnc,
         Dict(P[:S]=>N[:Qin], P[:T]=>N[:Box], P[:I]=>N[:Win], P[:O]=>N[:Wa]),
         Dict(P[:is]=>N[:arg],
              P[:it]=>N[:call],
              P[:ot]=>N[:influx],
              P[:os]=>N[:infusion]))
end

function migrate!(pn::AbstractPetriNet, bn::AbstractLabelledBilayerNetwork)
    P, N = Presentation.(acset_schema.([pn, bn])) 
    bnc = copy(bn)
    balance!(bnc)
    generators = Dict{Any,Any}(P[:S]=>N[:Qin], P[:T]=>N[:Box], P[:I]=>N[:Win], P[:O]=>N[:Wa])
    parts = Dict{Any,Any}(P[:is]=>N[:arg],
                 P[:it]=>N[:call],
                 P[:ot]=>N[:influx],
                 P[:os]=>N[:infusion])
    if :Name in attrtypes(acset_schema(pn))
        generators[P[:Name]] = N[:Name]
    end
    if has_subpart(pn, :sname)
        parts[P[:sname]] = N[:variable]
    end
    if has_subpart(pn, :tname)
        parts[P[:tname]] = N[:parameter]
    end
    migrate!(pn, bnc, generators, parts)
end


function propertygraph(bn::AbstractBilayerNetwork;
                       prog::AbstractString="dot", graph_attrs::AbstractDict=Dict(),
                       node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
                       node_labels::Bool=true, edge_labels::Bool=false)

    Qinshift = nparts(bn, :Box)
    Qoutshift = nparts(bn, :Box) + nparts(bn, :Qin)

    colors = ["#E28F41", "#6C9AC3", "#80BE63", "#A8DCD9", "#E0E5CD"]

    # vertex partition helpers V = coproduct(Box, Qin, Qout)
    isbox(v) = v <=Qinshift
    isqin(v) = Qinshift < v <= Qoutshift
    isqout(v) = Qoutshift < v

    # edge partition helpers E = coproduct(Win, W, Wn)
    calledge(e) = e <= nparts(bn, :Win)
    influxedge(e) = nparts(bn, :Win) < e <= nparts(bn, :Win) + nparts(bn, :Wa)
    effluxedge(e) = nparts(bn, :Win) + nparts(bn, :Wa) < e

    node_labeler(v) = if isbox(v)
        Dict(:label=>"p[$v]", :style=>"rounded, filled", :width=>"0.75", :height=>"0.5", :fixedsize=>"true", :shape=>"rectangle", :color=>colors[1])
    elseif isqin(v)
        return Dict(:label=>"u[$(v-Qinshift)]", :style=>"filled", :color=>colors[2])
    elseif isqout(v)
        return Dict(:label=>"u[$(v-Qoutshift)]'", :style=>"filled", :color=>colors[3])
    else
        return Dict(:label=>"$v", :style=>"filled", :color=>colors[4])
    end

    edge_labeler(e) = if calledge(e)
        return Dict(:label => edge_labels ? string(e) : "", :color=>"black", :style=>"dashed")
    elseif effluxedge(e)
        return Dict(:label => edge_labels ? string(e) : "", :color=>"red")
    else
        return Dict(:label => edge_labels ? string(e) : "", :color=>"black")
    end

    g′ = Catlab.Graphs.Graph()
    add_parts!(g′, :E, nparts(bn, :Win) + nparts(bn, :Wa) + nparts(bn, :Wn))
    add_parts!(g′, :V, nparts(bn, :Box) + nparts(bn, :Qin) + nparts(bn, :Qout))

    srcs = vcat(bn[:arg].+Qinshift, bn[:influx], bn[:efflux])
    tgts = vcat(bn[:call], bn[:infusion].+Qoutshift, bn[:effusion].+Qoutshift)

    g′[:, :src] = srcs
    g′[:, :tgt] = tgts
    Catlab.Graphs.PropertyGraph{Any}(g′, node_labeler, edge_labeler;
                                     prog = prog,
                                     graph = merge!(Dict(:rankdir => "TB"), graph_attrs),
                                     node = merge!(Graphics.GraphvizGraphs.default_node_attrs(node_labels), node_attrs),
                                     edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
                                     )
end

function GraphvizGraphs.to_graphviz(g::AbstractBilayerNetwork; kwargs...)
    to_graphviz(propertygraph(g; kwargs...))
end
function propertygraph(bn::AbstractLabelledBilayerNetwork;
                       prog::AbstractString="dot", graph_attrs::AbstractDict=Dict(),
                       node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
                       node_labels::Bool=true, edge_labels::Bool=false)

    Qinshift = nparts(bn, :Box)
    Qoutshift = nparts(bn, :Box) + nparts(bn, :Qin)

    colors = ["#E28F41", "#6C9AC3", "#80BE63", "#A8DCD9", "#E0E5CD"]

    # vertex partition helpers V = coproduct(Box, Qin, Qout)
    isbox(v) = v <=Qinshift
    isqin(v) = Qinshift < v <= Qoutshift
    isqout(v) = Qoutshift < v

    # edge partition helpers E = coproduct(Win, W, Wn)
    calledge(e) = e <= nparts(bn, :Win)
    influxedge(e) = nparts(bn, :Win) < e <= nparts(bn, :Win) + nparts(bn, :Wa)
    effluxedge(e) = nparts(bn, :Win) + nparts(bn, :Wa) < e

    node_labeler(v) = if isbox(v)
        Dict(:label=>"$(bn[v, :parameter])", :style=>"rounded, filled", :width=>"0.75", :height=>"0.5", :fixedsize=>"true", :shape=>"rectangle", :color=>colors[1])
    elseif isqin(v)
        return Dict(:label=>"$(bn[v-Qinshift, :variable])", :style=>"filled", :color=>colors[2])
    elseif isqout(v)
        return Dict(:label=>"$(bn[v-Qoutshift, :tanvar])'", :style=>"filled", :color=>colors[3])
    else
        return Dict(:label=>"$v", :style=>"filled", :color=>colors[4])
    end

    edge_labeler(e) = if calledge(e)
        return Dict(:label => edge_labels ? string(e) : "", :color=>"black", :style=>"dashed")
    elseif effluxedge(e)
        return Dict(:label => edge_labels ? string(e) : "", :color=>"red")
    else
        return Dict(:label => edge_labels ? string(e) : "", :color=>"black")
    end

    g′ = Catlab.Graphs.Graph()
    add_parts!(g′, :E, nparts(bn, :Win) + nparts(bn, :Wa) + nparts(bn, :Wn))
    add_parts!(g′, :V, nparts(bn, :Box) + nparts(bn, :Qin) + nparts(bn, :Qout))

    srcs = vcat(bn[:arg].+Qinshift, bn[:influx], bn[:efflux])
    tgts = vcat(bn[:call], bn[:infusion].+Qoutshift, bn[:effusion].+Qoutshift)

    g′[:, :src] = srcs
    g′[:, :tgt] = tgts
    Catlab.Graphs.PropertyGraph{Any}(g′, node_labeler, edge_labeler;
                                     prog = prog,
                                     graph = merge!(Dict(:rankdir => "TB"), graph_attrs),
                                     node = merge!(Graphics.GraphvizGraphs.default_node_attrs(node_labels), node_attrs),
                                     edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
                                     )
end

function GraphvizGraphs.to_graphviz(g::AbstractLabelledBilayerNetwork; kwargs...)
    to_graphviz(propertygraph(g; kwargs...))
end

function evaluate!(du, ϕ, bn::AbstractLabelledBilayerNetwork, state; params...)
    du.= 0.0
    ϕ .= 1.0
    for i in parts(bn, :Win)
        ϕ[bn[i, :call]] *= state[bn[i, :arg]]
    end
    for i in parts(bn, :Box)
        ϕ[i] *= params[bn[i, :parameter]]
    end

    for i in parts(bn, :Wn)
        du[bn[i, :effusion]] -= ϕ[bn[i, :efflux]]
    end
    for i in parts(bn, :Wa)
        du[bn[i, :infusion]] += ϕ[bn[i, :influx]]
    end
    return du
end

evaluate(bn::AbstractLabelledBilayerNetwork, state; params...) = evaluate!(zeros(length(state)), ones(nparts(bn, :Box)), bn, state; params...)

function paramexps(bn::AbstractLabelledBilayerNetwork, ϕ::Symbol, params::Symbol)
    map(parts(bn, :Box)) do i
        p = bn[i, :parameter]
        :($ϕ[$i] *= $params[$(Meta.quot(p))])
    end
end

function paramexps(bn::AbstractLabelledBilayerNetwork, ϕ::Symbol, params)
    map(parts(bn, :Box)) do i
        p = bn[i, :parameter]
        β = params[p]
        :($ϕ[$i] *= $β)
    end
end

funcwrap(du, ϕ, state, params::Symbol, body::Expr) = :(f!($du, $ϕ, $state, $params, t) = $body)
funcwrap(du, ϕ, state, params, body::Expr) = :(f!($du, $ϕ, $state, t) = $body)


function compile(bn::Union{AbstractLabelledBilayerNetwork, AbstractBilayerNetwork}, du::Symbol, ϕ::Symbol, state::Symbol, params)
    body = quote
        $du.= 0.0
        $ϕ .= 1.0
    end

    ϕs = map(parts(bn, :Win)) do i
        j = bn[i, :arg]
        k = bn[i, :call]
        :($ϕ[$k] *= $state[$j])
    end
    append!(body.args, ϕs)
    ps = paramexps(bn, ϕ, params)
    append!(body.args, ps)

    effs = map(parts(bn, :Wn)) do i
        j = bn[i, :efflux]
        k = bn[i, :effusion]
        :($du[$k] -= $ϕ[$j])
    end
    append!(body.args, effs)

    infs = map(parts(bn, :Wa)) do i
        j = bn[i,:influx]
        k = bn[i,:infusion]
        :($du[$k] += $ϕ[$j])
    end
    append!(body.args, infs)
    push!(body.args, :(return $du))
    return funcwrap(du, ϕ, state, params, body)
end

compile(bn, du, ϕ, state; params...) = compile(bn, du, ϕ, state, params)

end
