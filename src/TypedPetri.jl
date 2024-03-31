module TypedPetri
export prim_petri, strip_names, prim_cospan, oapply_typed, typed_product,
  add_params, add_reflexives, pairwise_id_petri, pairwise_id_typed_petri

using Catlab
using AlgebraicPetri
using AlgebraicPetri: LabelledPetriNetUntyped, OpenLabelledPetriNetUntyped, LabelledReactionNetUntyped


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
    i′ = add_species!(prim)
    push!(s_map, i)
    add_input!(prim, 1, i′)
  end
  for o in outputs(type_system, transition)
    o′ = add_species!(prim)
    push!(s_map, o)
    add_output!(prim, 1, o′)
  end
  LooseACSetTransformation(
    (S=s_map, T=[transition], O=incident(type_system, transition, :ot), I=incident(type_system, transition, :it)),
    (),
    prim,
    type_system
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
    LooseACSetTransformation(
      (S = [s], T = Int[], I = Int[], O = Int[]),
      (),
      structured_feet[s],
      p
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

The `tnames` argument renames the transitions in the resulting typed Petri net.
"""
function oapply_typed(type_system::LabelledPetriNet, uwd, tnames::Vector{Symbol})
  if junction(uwd, outer=true) != junctions(uwd)
    # XXX: This could be considered a user error, but for the sake of backwards
    # compatibility, we will fix it for them.
    uwd = copy(uwd)
    rem_parts!(uwd, :OuterPort, parts(uwd, :OuterPort))
    add_parts!(uwd, :OuterPort, njunctions(uwd), outer_junction=junctions(uwd))
  end
  type_system′ = PetriNet(type_system)
  prim_cospan_data = Dict(
    tname(type_system, t) => prim_cospan(type_system′, t)
    for t in 1:nt(type_system)
  )
  prim_cospans = Dict(t => prim_cospan_data[t][1] for t in keys(prim_cospan_data))
  petri, colim = oapply(uwd, prim_cospans; return_colimit=true)
  unlabelled_map = universal(
    colim,
    Multicospan(
      type_system′,
      [prim_cospan_data[uwd[b, :name]][2] for b in boxes(uwd)]
    )
  )
  state_labels = Vector{Symbol}(undef, njunctions(uwd))
  for (j, leg) in zip(junctions(uwd), legs(petri))
    state_labels[only(collect(leg[:S]))] = uwd[j, :variable]
  end
  labelled_petri = LabelledPetriNet(dom(unlabelled_map), state_labels, tnames)
  LooseACSetTransformation(
    components(unlabelled_map),
    (Name=x->nothing,),
    labelled_petri,
    type_system
  )
end

"""
Modify a typed petri net to add "reflexive transitions". These are transitions which go from one species to itself, so they don't change the mass action semantics, but they are important for stratification.

The idea behind this is similar to the fact that the product of the graph *----* with itself is

*    *
    /
   /
  /
 /
*    *

One needs to add self-edges to each of the vertices in *----* to get

*----*
|   /|
|  / |
| /  |
|/   |
*----*

Example:
```julia
add_reflexives(sird_typed, [[:strata], [:strata], [:strata], []], infectious_ontology)
```
"""
function add_reflexives(
  typed_petri::ACSetTransformation,
  reflexive_transitions,
  type_system::AbstractPetriNet
)
  petri = deepcopy(dom(typed_petri))
  type_comps = Dict(k=>collect(v) for (k,v) in pairs(deepcopy(components(typed_petri))))
  for (s_i,cts) in enumerate(reflexive_transitions)
    for ct in cts
      type_ind = findfirst(==(ct), type_system[:tname])
      is, os = [incident(type_system, type_ind, f) for f in [:it, :ot]]
      new_t = add_part!(petri, :T; tname=ct)
      if has_subpart(petri, :rate)
        set_subpart!(petri, new_t, :rate, 1.0)
      end
      add_parts!(petri, :I, length(is); is=s_i, it=new_t)
      add_parts!(petri, :O, length(os); os=s_i, ot=new_t)
      push!(type_comps[:T], type_ind)
      append!(type_comps[:I], is); append!(type_comps[:O], os);
    end
  end
  LooseACSetTransformation(
    type_comps,
    (Name=x->nothing, (has_subpart(petri, :rate) ? [:Rate=>x->nothing] : [])..., (has_subpart(petri, :concentration) ? [:Concentration=>x->nothing] : [])...),
    petri,
    codom(typed_petri)
  )
end

function add_params(
  typed_petri::ACSetTransformation,
  initial_concentrations::AbstractDict{Symbol,Float64},
  rates::AbstractDict{Symbol,Float64}
)
  domain = begin
    pn = LabelledReactionNet{Float64,Float64}()
    copy_parts!(pn, typed_petri.dom)
    for (tlabel, rate) in rates
      set_subpart!(pn, only(incident(typed_petri.dom, tlabel, :tname)), :rate, rate)
    end
    for (slabel, ic) in initial_concentrations
      set_subpart!(pn, only(incident(typed_petri.dom, slabel, :sname)), :concentration, ic)
    end
    pn
  end
  typing = begin
    pn = LabelledReactionNetUntyped{Nothing, Nothing, Nothing}()
    copy_parts!(pn, PetriNet(typed_petri.codom))
    for t in parts(pn, :T)
      set_subpart!(pn, t, :tname, nothing)
      set_subpart!(pn, t, :rate, nothing)
    end
    for s in parts(pn, :S)
      set_subpart!(pn, s, :sname, nothing)
      set_subpart!(pn, s, :concentration, nothing)
    end
    pn
  end
  LooseACSetTransformation(
    NamedTuple(map(objects(acset_schema(domain))) do ob
      ob => components(typed_petri)[ob]
    end),
    NamedTuple(map(attrtypes(acset_schema(domain))) do attrtype
      attrtype => x->nothing
    end),
    domain,
    typing,
  )
end

"""
This takes the "typed product" of two typed Petri nets p1 and p2, which has

- a species for every pair of a species in p1 and a species in p2 with the same type
- a transition for every pair of a transitions in p1 and a species in p2 with the same type

This is the "workhorse" of stratification; this is what actually does the stratification, though you may have to "prepare" the petri nets first with `add_reflexives`

Returns a typed model, i.e. a map in Petri.
"""
function typed_product(p1, p2)
  pb = pullback(p1, p2; product_attrs=true)
  return first(legs(pb)) ⋅ p1
end
function typed_product(ps::AbstractVector)
  pb = pullback(ps; product_attrs=true)
  return first(legs(pb)) ⋅ first(ps)
end

""" Make typed Petri net with 'identity' transformation between species pairs.

Assumes a single species type and a single transition type. For each pair of places, generate a transition
consuming 1 input for each place and producing 1 output in each place.
"""
function pairwise_id_typed_petri(type_net, stype, ttype, args...;
                                 codom_net=nothing) # TODO: this keyword parameter is a workaround. It should be removed and fixed later.
  net = pairwise_id_petri(args...)
  type_components = NamedTuple(map(attrtypes(acset_schema(net))) do type
    type => x->nothing
  end)
  s = only(incident(type_net, stype, :sname))
  t = only(incident(type_net, ttype, :tname))
  LooseACSetTransformation(
    (S = repeat([s], ns(net)), T = repeat([t], nt(net)), I = repeat(incident(type_net, t, :it), nt(net)), O = repeat(incident(type_net, t, :ot), nt(net))),
    type_components,
    net,
    isnothing(codom_net) ? type_net : codom_net
  )
end

""" Make Petri net with 'identity' transformation between all species pairs.
"""
function pairwise_id_petri(names::AbstractVector{Symbol})
  net = LabelledPetriNet()
  add_pairwise_id_transitions!(net, names)
  net
end
function pairwise_id_petri(names::AbstractVector{Symbol},
                           concentration_vec::AbstractVector{C},
                           contact_mat::AbstractMatrix{R}) where {R,C}
  n = length(names)
  @assert length(concentration_vec) == n
  @assert size(contact_mat) == (n, n)

  net = LabelledReactionNet{R,C}()
  add_pairwise_id_transitions!(net, names)
  net[:, :concentration] = concentration_vec
  net[:, :rate] = reshape(contact_mat, :)
  net
end

function add_pairwise_id_transitions!(net::AbstractPetriNet,
                                      names::AbstractVector{Symbol})
  n = length(names)
  species = add_species!(net, n, sname=names)
  for i in 1:n, j in 1:n
    tname = Symbol(string(names[i], "_", names[j]))
    t = add_transition!(net, tname=tname)
    add_input!(net, t, species[i])
    add_input!(net, t, species[j])
    add_output!(net, t, species[i])
    add_output!(net, t, species[j])
  end
  species
end

end
