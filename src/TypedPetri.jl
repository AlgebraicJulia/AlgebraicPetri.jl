module TypedPetri
export prim_petri, strip_names, prim_cospan, oapply_typed

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using AlgebraicPetri
using AlgebraicPetri: LabelledPetriNetUntyped, OpenLabelledPetriNetUntyped

# Dependencies for my_oapply
using Catlab.Present, Catlab.Theories
using Catlab.Schemas: dom_nums, codom_nums, attr
using Catlab.CategoricalAlgebra
import Catlab.CategoricalAlgebra.CSets: homomorphisms, homomorphism, is_homomorphic
using Catlab.WiringDiagrams.UndirectedWiringDiagrams

""" Compose morphisms according to UWD.

The morphisms corresponding to the boxes, and optionally also the objects
corresponding to the junctions, are given by dictionaries indexed by
box/junction attributes. The default attributes are those compatible with the
`@relation` macro.
"""
function my_oapply(composite::UndirectedWiringDiagram, hom_map::AbstractDict,
                ob_map::Union{AbstractDict,Nothing}=nothing;
                hom_attr::Symbol=:name, ob_attr::Symbol=:variable)
  # XXX: Julia should be inferring these vector eltypes but isn't on v1.7.
  homs = valtype(hom_map)[ hom_map[name] for name in composite[hom_attr] ]
  obs = isnothing(ob_map) ? nothing :
    valtype(ob_map)[ ob_map[name] for name in composite[ob_attr] ]
  my_oapply(composite, homs, obs)
end

function my_oapply(composite::UndirectedWiringDiagram,
                cospans::AbstractVector{<:StructuredMulticospan{L}},
                junction_feet::Union{AbstractVector,Nothing}=nothing) where L
  @assert nboxes(composite) == length(cospans)
  if isnothing(junction_feet)
    junction_feet = Vector{first(dom(L))}(undef, njunctions(composite))
  else
    @assert njunctions(composite) == length(junction_feet)
  end

  # Create bipartite free diagram whose vertices of types 1 and 2 are the UWD's
  # junctions and boxes, respectively.
  diagram = BipartiteFreeDiagram{codom(L)...}()
  add_vertices₁!(diagram, njunctions(composite))
  add_vertices₂!(diagram, nboxes(composite), ob₂=map(apex, cospans))
  for (b, cospan) in zip(boxes(composite), cospans)
    for (p, leg, foot) in zip(ports(composite, b), legs(cospan), feet(cospan))
      j = junction(composite, p)
      add_edge!(diagram, j, b, hom=leg)
      if isassigned(junction_feet, j)
        foot′ = junction_feet[j]
        foot == foot′ || error("Feet of cospans are not equal: $foot != $foot′")
      else
        junction_feet[j] = foot
      end
    end
  end
  for (j, foot) in enumerate(junction_feet)
    diagram[j, :ob₁] = L(foot)
  end

  # Find, or if necessary create, an outgoing edge for each junction. The
  # existence of such edges is an assumption for colimits of bipartite diagrams.
  # The edges are also needed to construct inclusions for the outer junctions.
  junction_edges = map(junctions(composite)) do j
    out_edges = incident(diagram, j, :src)
    if isempty(out_edges)
      x = ob₁(diagram, j)
      v = add_vertex₂!(diagram, ob₂=x)
      add_edge!(diagram, j, v, hom=id(x))
    else
      first(out_edges)
    end
  end

  # The composite multicospan is given by the colimit of this diagram.
  colim = colimit(diagram)
  outer_js = junction(composite, outer=true)
  outer_legs = map(junction_edges[outer_js]) do e
    hom(diagram, e) ⋅ legs(colim)[tgt(diagram, e)]
  end
  outer_feet = junction_feet[outer_js]
  (colim, StructuredMulticospan{L}(Multicospan(ob(colim), outer_legs), outer_feet))
end

"""
Takes in a petri net and a transition in that petri net,
constructs a petri net with just that transition, returns
the acsettransformation from that into the type system.
"""
function prim_petri(type_system, transition)
  prim = typeof(type_system)()
  t = add_transition!(prim)
  s_map = Int[]
  for i in inputs(type_system, transition)
    s = is(type_system, i)
    s′ = add_species!(prim)
    push!(s_map, s)
    add_input!(prim, 1, s′)
  end
  for i in outputs(type_system, transition)
    s = os(type_system, i)
    s′ = add_species!(prim)
    push!(s_map, s)
    add_output!(prim, 1, s′)
  end
  ACSetTransformation(
    prim,
    type_system,
    S=s_map,
    T=[transition],
    O=outputs(type_system, transition),
    I=inputs(type_system, transition),
  )
end

"""
Produces the structured cospan of the transition over
different species.
"""
function prim_cospan(type_system, transition)
  prim = prim_petri(type_system, transition)
  p = dom(prim)
  feet = [FinSet(1) for s in 1:ns(p)]
  structured_feet = [PetriNet(1) for s in 1:ns(p)]
  legs = [
    ACSetTransformation(
      structured_feet[s],
      p,
      S = [s],
      T = Int[],
      I = Int[],
      O = Int[],
    )
    for s in 1:ns(p)
  ]
  (OpenPetriNet(Multicospan(p, legs), feet), prim)
end

"""
Takes in a labelled petri net and an undirected wiring diagram,
where each of the boxes is labeled by a symbol that matches the label
of a transition in the petri net. Then produces a petri net given by
colimiting the transitions together, and returns the ACSetTransformation
from that Petri net to the type system.
"""
function oapply_typed(type_system::LabelledPetriNet, uwd)
  type_system′ = PetriNet(type_system)
  prim_cospan_data = Dict(
    tname(type_system, t) => prim_cospan(type_system′, t)
    for t in 1:nt(type_system)
  )
  prim_cospans = Dict(t => prim_cospan_data[t][1] for t in keys(prim_cospan_data))
  (colim, petri) = my_oapply(uwd, prim_cospans)
  universal(
    colim,
    Multicospan(
      type_system′,
      [prim_cospan_data[subpart(uwd, b, :name)][2] for b in 1:nboxes(uwd)]
    )
  )
end

end
