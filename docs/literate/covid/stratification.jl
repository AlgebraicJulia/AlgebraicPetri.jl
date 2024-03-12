# # [Stratification of Epidemiological Models](@id stratification_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/covid/stratification.ipynb)

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using DisplayAs, Markdown

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

# In this tutorial we show how to use this package to define basic epidemiological models and 
# model stratification, when a submodel is replicated across a variable and allowed to interact
# in some specified way (e.g.; age-structured infection dynamics). For more details on the underlying
# theory, please see the research article ["An algebraic framework for structured epidemic modelling"](https://doi.org/10.1098/rsta.2021.0309).

# ## Define basic epidemiology model

# One of the fundamental parts of getting model stratification to work is defining a "type system"
# for the models one is interested in, which ensures that spurious or erroneous transitions are
# not generated when we construct the stratified model (e.g.; no vector to vector transmission in a malaria model).
# In the research article this tutorial follows, the Petri net that defines the type system is ``P_{type}``.
# We define ``P_{infectious}`` below.

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
  :disease => (:Pop => :Pop),
  :strata => (:Pop => :Pop)
)

to_graphviz(infectious_ontology)

# More details are given in the article, but briefly, ``P_{infectious}`` has a single place, which represents
# a potentially stratified population, and 3 "types" of transitions: "infect" is for infections, "disease" for changes
# in disease state (e.g.; from E to I), and "strata" for changes not in disease state (e.g.; from region ``i`` to region ``j``).
# 
# A typed Petri net is a morphism ``\phi: P \to P_{type}``. Here we define a Petri net representing the SIRD model
# typed over ``P_{infectious}``. We start by defining an undirected wiring diagram (UWD) which represents the processes
# and states in the SIRD model; note that the name of each box corresponds to one of the transitions in ``P_{infectious}``,
# and that we specify the types of the junctions as `Pop` in the `where` clause.
# 
# Next we use `oapply_typed` which produces a typed Petri net by composing transitions based on ``P_{infectious}``.
# The typed Petri net is modified with reflexive transitions typed as `:strata` to indicate which states can be stratified
# Here we add reflexive transitions to the susceptible, infected, and recovered populations but we leave out the dead
# population as no individuals entering that compartment may leave.

sird_uwd = @relation (S,I,R,D) where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S, I, I, I)
  disease(I, R)
  disease(I, D)
end

sird_model = oapply_typed(infectious_ontology, sird_uwd, [:infection, :recovery, :death])
sird_model = add_reflexives(sird_model, [[:strata], [:strata], [:strata], []], infectious_ontology)

to_graphviz(dom(sird_model))

# ## Define intervention models

# ### Masking model

# Let's say we want to stratify the populations in the SIRD model by masking or non-masking; or, generally
# by any behavioral change which is reversible and does not preclude the possibility of infection in either
# state. To generate our stratification scheme, we begin with a UWD describing the strata levels we are
# interested in (Masked and UnMasked), and the transitions possible in each level. Note that we give
# transitions of type `strata` which allow masking or unmasking. Unmasked individuals may transmit
# disease to either masked or unmasked individuals, but masked persons cannot transmit disease.
# 
# We again use `oapply_typed` to generate the stratification scheme ``\phi^{'}:P^{'}\to P_{infectious}``,
# adding reflexive `disease` transitions to the typed Petri net as changes in disease status may occur
# in either stratum.

masking_uwd = @relation (M,UM) where (M::Pop, UM::Pop) begin
  strata(M, UM)
  strata(UM, M)
  infect(M, UM, M, UM)
  infect(UM, UM, UM, UM)
end
mask_model = oapply_typed(infectious_ontology, masking_uwd, [:unmask, :mask, :infect_um, :infect_uu])
mask_model = add_reflexives(mask_model, [[:disease], [:disease]], infectious_ontology)

to_graphviz(dom(mask_model))

# To generate the stratified SIRD model, we use `typed_product`. As described in the article, a stratification of
# a model can be seen as a pullback of ``\phi`` and ``\phi^{'}``, or a product in the slice category ``Petri/P_{type}``.
# Stratify our SIRD model on this masking model to get a model of SIRD with masking:

typed_product(sird_model, mask_model) |> dom |> to_graphviz

# ### Vaccine model

vax_uwd = @relation (UV,V) where (UV::Pop, V::Pop) begin
  strata(UV, V)
  infect(V, V, V, V)
  infect(V, UV, V, UV)
  infect(UV, V, UV, V)
  infect(UV, UV, UV, UV)
end
vax_model = oapply_typed(infectious_ontology, vax_uwd, [:vax, :infect_vv, :infect_uv, :infect_vu, :infect_uu])
vax_model = add_reflexives(vax_model, [[:disease], [:disease]], infectious_ontology)

to_graphviz(dom(vax_model))

# Stratify our SIRD model on this vaccine model to get a model of SIRD with a vaccination rate:

typed_product(sird_model, vax_model) |> dom |> to_graphviz

# ### Mask-Vax Model

mask_vax_uwd = @relation (UV_UM,UV_M,V_UM,V_M) where (UV_UM::Pop, UV_M::Pop, V_UM::Pop, V_M::Pop) begin
  strata(UV_UM, UV_M)
  strata(UV_M, UV_UM)
  strata(V_UM, V_M)
  strata(V_M, V_UM)
  strata(UV_UM, V_UM)
  strata(UV_M, V_M)
  infect(V_UM, V_UM, V_UM, V_UM)
  infect(V_UM, UV_UM, V_UM, UV_UM)
  infect(UV_UM, V_UM, UV_UM, V_UM)
  infect(UV_UM, UV_UM, UV_UM, UV_UM)
  infect(V_M, V_UM, V_M, V_UM)
  infect(V_M, UV_UM, V_M, UV_UM)
  infect(UV_M, V_UM, UV_M, V_UM)
  infect(UV_M, UV_UM, UV_M, UV_UM)
end
mask_vax_model = oapply_typed(
  infectious_ontology,
  mask_vax_uwd,
  [:mask_uv, :unmask_uv, :mask_v, :unmask_v, :vax_um, :vax_m, :infect_um_vv, :infect_um_uv, :infect_um_vu, :infect_um_uu, :infect_m_vv, :infect_m_uv, :infect_m_vu, :infect_m_uu]
)
mask_vax_model = add_reflexives(mask_vax_model, [[:disease], [:disease], [:disease], [:disease]], infectious_ontology)

to_graphviz(dom(mask_vax_model))

# Stratify our SIRD model on this mask + vaccine model to get a model of SIRD with a vaccination rate and masking policies:

typed_product(sird_model, mask_vax_model) |> dom |> to_graphviz

# ## Define geographic models

# ### Travel model between $N$ regions

# For this model we can use a julia function to programmatically build up our undirected wiring diagram for defining this model.
# Here we want there to be $N$ regions in which people can travel between each region and people within the same region are able
# to infect other people in the same region.

function travel_model(n)
  uwd = RelationDiagram(repeat([:Pop], n))
  junctions = Dict(begin
    variable = Symbol("Region$(i)")
    junction = add_junction!(uwd, :Pop, variable=variable)
    set_junction!(uwd, port, junction, outer=true)
    variable => junction
  end for (i, port) in enumerate(ports(uwd, outer=true)))

  pairs = filter(x -> first(x) != last(x), collect(Iterators.product(keys(junctions), keys(junctions))))
  for pair in pairs
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:strata)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end

  act = oapply_typed(infectious_ontology, uwd, [Symbol("$(a)_$(b)") for (a, b) in pairs])
  add_reflexives(act, repeat([[:infect, :disease]], n), infectious_ontology)
end

to_graphviz(dom(travel_model(2)))

# Stratify our SIRD model on this travel model with two regions:

typed_product(sird_model, travel_model(2)) |> dom |> to_graphviz

# ### Simple Trip model between $N$ regions

# For this model we can use a julia function to programmatically build up our model where people have the property of living somewhere
# and we are modelling them travelling between locations while maintaining the status of where they live.  Here we can actually just
# define the model of having a "Living" status and stratify it with the previously defined travel model to get a model of someone taking a simple trip.

function living_model(n)
  typed_living = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, [Symbol("Living$(i)") for i in 1:n])
  add_reflexives(typed_living, repeat([[:disease, :strata]], n), infectious_ontology)
end

to_graphviz(dom(living_model(2)))

# The resulting simple trip model:

simple_trip_model = typed_product(travel_model(2), living_model(2))
to_graphviz(dom(simple_trip_model))

# Stratify our SIRD model on this simple trip model between two regions:

typed_product(sird_model, simple_trip_model) |> dom |> to_graphviz

# ## Stratification of COVID models

# we set up a simple helper function to connect the undirected wiring diagrams to our
# infectious disease type system and add the necessary reflexive transitions for stratification.
function oapply_mira_model(uwd) 
  model = oapply_typed(infectious_ontology, uwd, [Symbol("t$(n)") for n in 1:nboxes(uwd)])
  add_reflexives(model, [repeat([[:strata]], njunctions(uwd)-3)..., [], [:strata],[]], infectious_ontology)
end

# ### SIDARTHE Model
#
# BIOMD0000000955_miranet

m1_model = (@relation (S,I,D,A,R,T,H,E) where (S::Pop, I::Pop, D::Pop, A::Pop, R::Pop, T::Pop, H::Pop, E::Pop) begin
  infect(S, D, I, D)
  infect(S, A, I, A)
  infect(S, R, I, R)
  infect(S, I, I, I)
  disease(I, D)
  disease(I, A)
  disease(I, H)
  disease(D, R)
  disease(D, H)
  disease(A, R)
  disease(A, H)
  disease(A, T)
  disease(R, T)
  disease(R, H)
  disease(T, E)
  disease(T, H)
end) |> oapply_mira_model

to_graphviz(dom(m1_model))

# ### SEIAHRD Model
#
# BIOMD0000000960_miranet

m2_model = (@relation (S,E,I,A,H,R,D) where (S::Pop, E::Pop, I::Pop, A::Pop, H::Pop, R::Pop, D::Pop) begin
  infect(S, I, E, I)
  infect(S, A, E, A)
  infect(S, H, E, H)
  disease(E, I)
  disease(E, A)
  disease(I, H)
  disease(I, R)
  disease(I, D)
  disease(A, R)
  disease(A, D)
  disease(H, D)
  disease(H, R)
end) |> oapply_mira_model

to_graphviz(dom(m2_model))

# ### SEIuIrQRD Model
#
# BIOMD0000000983_miranet

m3_model = (@relation (S,E,Iu,Ir,Q,R,D) where (S::Pop, E::Pop, Iu::Pop, Ir::Pop, Q::Pop, R::Pop, D::Pop) begin
  infect(S, Ir, E, Ir)
  infect(S, Iu, E, Iu)
  infect(S, Ir, Q, Ir)
  infect(S, Iu, Q, Iu)
  disease(Q, Ir)
  disease(E, Ir)
  disease(E, Iu)
  disease(Ir, R)
  disease(Iu, R)
  disease(Q, S)
  disease(Ir, D)
end) |> oapply_mira_model

to_graphviz(dom(m3_model))

# ### Enumerating stratified models
#
# Next we can take all of our epidemiology, intervention, and geography models and easily enumerate
# and calculate the models for every possible stratification combination and investigate the resulting models
# with some set number of regions for the geographic stratification models.

num_rgns = 5
disease_models = [("SIRD", sird_model), ("SIDARTHE", m1_model), ("SEIAHRD", m2_model), ("SEIuIrQRD", m3_model)]
policy_models = [nothing, ("Vaccination", vax_model), ("Masking", mask_model), ("Masking + Vaccination", mask_vax_model)]
travel_models = [nothing, ("Travel", travel_model(num_rgns)), ("Simple Trip", typed_product(travel_model(num_rgns), living_model(num_rgns)))]

table = ["| Model | Intervention | Geography ($(num_rgns) regions) | # of States | # of Transitons |","|:--|$(repeat(":-:|", 4))"]

for pieces in Iterators.product(disease_models, policy_models, travel_models)
  petri = typed_product(last.(collect(filter(x -> !isnothing(x), pieces)))) |> dom
  push!(table, "|$(join([isnothing(piece) ? "N/A" : first(piece) for piece in pieces],"|"))|$(ns(petri))|$(nt(petri))|")
end

Markdown.parse(join(table, "\n")) |> DisplayAs.HTML

# ## Performance
#
# As we increase the number of regions in our geographic stratification models, the
# number of states and transitions increase polynomially which causes the execution time
# for calculating the final stratified model to also increase polynomially.
#
# ![Runtime vs. Number of Georgraphic Regions](../../assets/runtime_vs_num_rgns.svg)
