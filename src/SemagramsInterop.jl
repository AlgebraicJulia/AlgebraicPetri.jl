using ...Semagrams
import ...Semagrams.UI: load, Semagram
using ...Semagrams.Data: IDGen, LocatedSemagramData
using Catlab
using JSON
export LabelledReactionSema

@semagramschema PetriNetSema(TheoryPetriNet) begin
  @box S Circle
  @box T Square
  @wire I(is,it)
  @wire O(ot,os)
end;

@semagramschema LabelledPetriNetSema(TheoryLabelledPetriNet) begin
  @box S Circle :sname
  @box T Square :tname
  @wire I(is,it)
  @wire O(ot,os)
  @data Name Stringlike
end;

@semagramschema ReactionNetSema(TheoryReactionNet) begin
  @box S Circle :concentration
  @box T Square :rate
  @wire I(is,it)
  @wire O(ot,os)
  @data Rate Numeric
  @data Concentration Numeric
end;

@semagramschema LabelledReactionNetSema(TheoryLabelledReactionNet) begin
  @box S Circle :sname
  @box T Square :tname
  @wire I(is,it)
  @wire O(ot,os)
  @data Name Stringlike
  @data Rate Numeric
  @data Concentration Numeric
end;
load(sg, pn::AbstractPetriNet) = load(sg, LocatedSemagramData(pn))

Semagram(a::T) where {T} = Semagram{T}(LocatedSemagramData(a))

LocatedSemagramData(pn::PetriNet) =
  LocatedSemagramData(pn, PetriNetSema)
LocatedSemagramData(pn::LabelledPetriNet) =
  LocatedSemagramData(pn, LabelledPetriNetSema)
LocatedSemagramData(pn::ReactionNet) =
  LocatedSemagramData(pn, ReactionNetSema)
LocatedSemagramData(pn::LabelledReactionNet) =
  LocatedSemagramData(pn, LabelledReactionNetSema)

function LocatedSemagramData(pn::AbstractPetriNet, schema::SemagramSchema)
  g_boxes = JSON.parse(read(Catlab.Graphics.Graphviz.run_graphviz(Graph(pn),
                                                                  format="json"),
                            String))["objects"];
  b_ind = Dict(b["name"] =>
               Tuple{Float64, Float64}(parse.(Float64, split(b["pos"], ",")))
               for b in g_boxes)
  y_locs = last.(values(b_ind))
  y_shift = 300 - (maximum(y_locs) + minimum(y_locs))

  # Scale returned coordinates
  for (k,v) in b_ind
    b_ind[k] = 2 .* v .+ (80.0, y_shift)
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
                         Dict{Int, Port}())) for s in 1:nparts(pn, k)]
  end...))

  loc_dict = Dict(vcat(map(collect(schema.box_types)) do (k, v)
    [(k_offset[k] + i, loc) for (i, loc) in enumerate(locs[k])]
  end...))

  wires = Dict(vcat(map(collect(schema.wire_types)) do (k, v)
    [(uid(), Wire(k, Dict{Symbol, String}(w=>"$(pn[s, w])" for w in last.(v.weights)),
          Semagrams.Data.AttachBox(pn[i, v.src_map] + k_offset[v.src[2]]),
          Semagrams.Data.AttachBox(pn[i, v.tgt_map] + k_offset[v.tgt[2]])))
     for i in 1:nparts(pn, k)]
  end...))

  LocatedSemagramData(SemagramData(boxes, wires, Semagrams.Data.IDGen(uid()),
                                   schema), loc_dict)
end
