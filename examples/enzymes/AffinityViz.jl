using Colors
using AlgebraicPetri
using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using LabelledArrays

@present ThAffinityNet(FreeSchema) begin
  (Rxn, Species)::Ob
  (Rate, Label)::Data
  binder::Hom(Rxn, Species)
	bindee::Hom(Rxn, Species)

  affinity::Attr(Rxn, Rate)

  family::Attr(Species, Label)
  form::Attr(Species, Label)
end

const AbstractAffinityNet = AbstractACSetType(ThAffinityNet)
const AffinityNet = ACSetType(ThAffinityNet){Real, Symbol}

function migrate_species!(an::AffinityNet, ln::LabelledReactionNet, index::Int64)
  labels = Symbol.(split(string(ln[index, :sname]), "_"))
  if length(labels) == 1
    push!(labels, :norm)
  end
  @assert (length(labels) == 2)  begin
    error("Label \"$(ln[index,:sname])\" is not in a valid format (family_form)")
  end
  add_part!(an, :Species, family=labels[1], form=labels[2])
end

function AffinityNet(ln::LabelledReactionNet{R, C}) where {R,C}
  # Collect associative/dissociative pairs
  assoc  = Dict{Int64, Array{Tuple{Int64, Int64, Int64},1}}()
  dissoc = Dict{Int64, Array{Tuple{Int64, Int64, Int64},1}}()
  for o in 1:nparts(ln, :T)
    outputs = ln[incident(ln, o, :ot), :os]
    inputs = ln[incident(ln, o, :it), :is]

    if length(outputs) == 2 && length(inputs) == 1
      if !(inputs[1] in keys(dissoc))
        dissoc[inputs[1]] = Array{Tuple{Int64, Int64, Int64},1}()
      end
      push!(dissoc[inputs[1]], outputs[1] < outputs[2] ? (outputs[1], outputs[2], o) :
                                                         (outputs[2], outputs[1], o))
    elseif length(inputs) == 2 && length(outputs) == 1
      if !(outputs[1] in keys(assoc))
        assoc[outputs[1]] = Array{Tuple{Int64, Int64, Int64},1}()
      end
      push!(assoc[outputs[1]], inputs[1] < inputs[2] ? (inputs[1], inputs[2], o) :
            (inputs[2], inputs[1], o))
    end
  end

  # Generate affinity net from pairs
  an = AffinityNet()
  s2i = Dict{Int64, Int64}()
  for prod in keys(assoc)
    if prod in keys(dissoc)
      for j in 1:length(assoc[prod])
        for i in 1:length(dissoc[prod])
          if assoc[prod][j][[1,2]] == dissoc[prod][i][[1,2]]
            ls, rs = ln[incident(ln, assoc[prod][j][3], :it), :is]
            ar, dr = (ln[assoc[prod][j][3], :rate], ln[dissoc[prod][i][3], :rate])
            if !(ls in keys(s2i))
              s2i[ls] = migrate_species!(an, ln, ls)
            end
            if !(rs in keys(s2i))
              s2i[rs] = migrate_species!(an, ln, rs)
            end
            add_part!(an, :Rxn, binder=s2i[ls], bindee=s2i[rs], affinity=ar/dr)
          end
        end
      end
    end
  end
  an
end

function propertygraph(an::AffinityNet;
                       prog::AbstractString="neato", graph_attrs::AbstractDict=Dict(),
                       node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(:len=>"2.0"),
                       node_labels::Bool=true, edge_labels::Bool=false)

  shapes = ["box", "circle", "triangle", "diamond", "pentagon", "hexagon", "septagon", "octagon"]

  families = unique(an[:family])
  forms =	unique(an[:form])

  colors = string.(hex.(distinguishable_colors(length(forms), [RGB(1,1,1),RGB(0,0,0)], dropseed=true)))

  fam2shp = Dict(families[i]=>shapes[i] for i in 1:length(families))
  frm2clr = Dict(forms[i]=>colors[i] for i in 1:length(forms))

  rates = an[:affinity]
  minrate, maxrate = (minimum(rates), maximum(rates))
  rate_to_width(rate) = begin
    0.1 + (log(rate) - log(minrate)) * (2/(log(maxrate)-log(minrate)))
  end

	node_labeler(v) = begin
    Dict(:label=>"$(an[v,:family])$(an[v,:form]==:norm ? "" : Symbol("_",an[v,:form]))",
         :style=>"filled",
         :width=>"0.75", :height=>"0.75", :fixedsize=>"false",
         :shape=>fam2shp[an[v,:family]],
         :fillcolor=>"#$(frm2clr[an[v,:form]])")
  end

  edge_labeler(e) = begin
    return Dict(:color=>"black", :penwidth=>"$(round(rate_to_width(an[e,:affinity]), digits=1))")
  end

  g = Catlab.Graphs.Graph()
  migrate!(g, an, Dict(:V=>:Species, :E=>:Rxn), Dict(:tgt=>:bindee, :src=>:binder))

  Catlab.Graphs.PropertyGraph{Any}(g, node_labeler, edge_labeler;
                                  prog = prog,
                                  graph = merge!(Dict(:rankdir => "TB"), graph_attrs),
                                  node = merge!(Graphics.GraphvizGraphs.default_node_attrs(node_labels), node_attrs),
                                  edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
                                  )
end

function Graphics.to_graphviz(an::AffinityNet; kwargs...)
  propertygraph(an; kwargs...) |> to_graphviz
end
