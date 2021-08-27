using ...Semagrams
import ...Semagrams.UI: load, Semagram
using ...Semagrams.Data: IDGen, LocatedSemagramData
using Catlab
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs
using JSON
export LabelledReactionSema, locations

@semagramschema PetriNetSema(TheoryPetriNet) begin
  @box S Circle
  @box T Square
  @wire I(is,it)
  @wire O(ot,os)
end;

@semagramschema LabelledPetriNetSema(TheoryLabelledPetriNet) begin
  @box S Circle label = :sname
  @box T Square label = :tname
  @wire I(is,it)
  @wire O(ot,os)
  @data Name Stringlike
end;

@semagramschema ReactionNetSema(TheoryReactionNet) begin
  @box S Circle label = :concentration
  @box T Square label = :rate
  @wire I(is,it)
  @wire O(ot,os)
  @data Rate Numeric
  @data Concentration Numeric
end;

@semagramschema LabelledReactionNetSema(TheoryLabelledReactionNet) begin
  @box S Circle label = :sname
  @box T Square label = :tname
  @wire I(is,it)
  @wire O(ot,os)
  @data Name Stringlike
  @data Rate Numeric
  @data Concentration Numeric
end;
load(sg, pn::AbstractPetriNet; kw...) = load(sg, LocatedSemagramData(pn; kw...))

Semagram(a::T; kw...) where {T} = Semagram{T}(LocatedSemagramData(a; kw...))

LocatedSemagramData(pn::PetriNet; kw...) =
  LocatedSemagramData(pn, PetriNetSema; kw...)
LocatedSemagramData(pn::LabelledPetriNet; kw...) =
  LocatedSemagramData(pn, LabelledPetriNetSema; kw...)
LocatedSemagramData(pn::ReactionNet; kw...) =
  LocatedSemagramData(pn, ReactionNetSema; kw...)
LocatedSemagramData(pn::LabelledReactionNet; kw...) =
  LocatedSemagramData(pn, LabelledReactionNetSema; kw...)

function LocatedSemagramData(pn::AbstractPetriNet, schema::SemagramSchema;
                             y_center = 300, scale=2.0, translate=(80.0, 0.0),
                             kw...)
  pg = read(run_graphviz(Graph(pn; kw...),format="json"), String) |>
                                                JSON.parse |>
                                                parse_graphviz
  b_ind = Dict(zip(get_vprop.([pg], 1:nv(pg), :name),
                   Tuple.(get_vprop.([pg], 1:nv(pg), :position))))

  y_locs = last.(values(b_ind))
  y_shift = y_center
  y_mean = (maximum(y_locs) + minimum(y_locs)) / 2

  # Scale returned coordinates
  for (k,v) in b_ind
    #First flip y-axis of y-centered graph
    b_ind[k] = (1,-1) .* (b_ind[k] .- (0, y_mean))
    b_ind[k] = scale .* b_ind[k] .+ (0.0, y_shift) .+ translate
  end

  locs = Dict(:T=>[b_ind["t$t"] for t in 1:nt(pn)],
              :S=>[b_ind["s$s"] for s in 1:ns(pn)])
  LocatedSemagramData(pn, schema, locs)
end

function LocatedSemagramData(pn, schema::SemagramSchema,
	locs::Dict{Symbol, Vector{Tuple{Float64, Float64}}})
  l = -1
  uid() = (l += 1)

  k_offset = Dict{Symbol, Int64}()
  boxes = Dict(vcat(map(collect(schema.box_types)) do (k, v)
    k_offset[k] = l
    [(uid(), Box(k, Dict{Symbol, String}(w=>"$(pn[s, w])" for w in last.(v.weights)),
                    Dict{Symbol, Semagrams.Data.Entity}(),
                    Dict{Int, Port}())) for s in 1:nparts(pn, k)]
  end...))

  loc_dict = Dict(vcat(map(collect(schema.box_types)) do (k, v)
    [(k_offset[k] + i, loc) for (i, loc) in enumerate(locs[k])]
  end...))

  wires = Dict(vcat(map(collect(schema.wire_types)) do (k, v)
    [(uid(), Wire(k, Dict{Symbol, String}(w=>"$(pn[s, w])" for w in last.(v.weights)),
                     Dict{Symbol, Semagrams.Data.Entity}(),
                     Semagrams.Data.BoxEntity(pn[i, v.src_map] + k_offset[v.src[2]]),
                     Semagrams.Data.BoxEntity(pn[i, v.tgt_map] + k_offset[v.tgt[2]])))
     for i in 1:nparts(pn, k)]
  end...))

  LocatedSemagramData(SemagramData(boxes, wires, Semagrams.Data.IDGen(uid()),
                                   schema), loc_dict)
end

function scale_locs(locs::Vector{Tuple{Float64, Float64}}; scale=72 * 2)
  # Scale the locs
  locs = map(v->v ./ scale, locs)

  # Rotate over y-axis
  y_mean = (maximum(last.(locs)) + minimum(last.(locs))) / 2.0
  locs = map(locs) do l
    (1,-1) .* (l .- (0, y_mean))
  end

  # Shift back to 0,0
  shift = (minimum(first.(locs)), minimum(last.(locs)))
  map(v-> v .- shift, locs)
end

function locations(pn::AbstractPetriNet, lsd::LocatedSemagramData; kw...)
  blocs = lsd.boxlocs
  sclocs = scale_locs(collect(values(lsd.boxlocs)); kw...)
  for (i,k) in enumerate(keys(blocs))
    blocs[k] = sclocs[i]
  end

  boxes = sort(collect(lsd.sg.boxes); by=first)
  locs = Dict{Symbol, Vector{Tuple{Float64, Float64}}}(:T=>[], :S=>[])

  for (id, b) in boxes
    push!(locs[b.ty], locs[id])
  end
  locs
end

function locations(pn::Union{LabelledReactionNet, LabelledPetriNet},
                   lsd::LocatedSemagramData; kw...)
  blocs = lsd.boxlocs
  sclocs = scale_locs(collect(values(lsd.boxlocs)); kw...)
  for (i,k) in enumerate(keys(blocs))
    blocs[k] = sclocs[i]
  end

  boxes = sort(collect(lsd.sg.boxes); by=first)
  t_boxes = filter(b -> b[2].ty == :T, boxes)
  s_boxes = filter(b -> b[2].ty == :S, boxes)

  t_label2id = Dict(b[2].weights[:tname]=>b[1] for b in t_boxes)
  s_label2id = Dict(b[2].weights[:sname]=>b[1] for b in s_boxes)

  Dict(:T=>[blocs[t_label2id["$t"]] for t in tnames(pn)],
       :S=>[blocs[s_label2id["$s"]] for s in snames(pn)])
end

Graph(pn::AbstractPetriNet, sg::Semagram; kw...) =
  Graph(pn, save(sg); kw...)
function Graph(pn::AbstractPetriNet, lsd::LocatedSemagramData; scale=72*2, kw...)
  Graph(pn; positions = locations(pn, lsd; scale=scale), kw...)
end
