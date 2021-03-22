module BilayerNetworks

using AlgebraicPetri
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
import Catlab.CategoricalAlgebra: migrate!

export ThBilayerNetwork, AbstractBilayerNetwork, BilayerNetwork

@present ThBilayerNetwork(FreeSchema) begin
    (Qin, Qout, Win, Wn, Wa, Box)::Ob
    arg::Hom(Win, Qin)
    call::Hom(Win, Box)
    influx::Hom(Wa, Box)
    infusion::Hom(Wa, Qout)
    efflux::Hom(Wn, Box)
    effusion::Hom(Wn, Qout)
end


const AbstractBilayerNetwork = AbstractACSetType(ThBilayerNetwork)
const BilayerNetwork = ACSetType(ThBilayerNetwork)

function migrate!(bn::AbstractBilayerNetwork, pn::AbstractPetriNet)
    migrate!(bn, pn,
             Dict(:Qin=>:S, :Qout=>:S, :Box=>:T, :Win=>:I, :Wn=>:I, :Wa=>:O),
             Dict(:arg=>:is,
                  :call=>:it,
                  :efflux=>:it,
                  :effusion=>:is,
                  :influx=>:ot,
                  :infusion=>:os))
end

function migrate!(pn::AbstractPetriNet, bn::AbstractBilayerNetwork)
    migrate!(pn,bn,
         Dict(:S=>:Qin, :T=>:Box, :I=>:Win, :O=>:Wa),
         Dict(:is=>:arg,
              :it=>:call,
              :ot=>:influx,
              :os=>:infusion))
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

end
