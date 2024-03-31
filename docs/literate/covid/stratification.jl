# # [Stratification of Epidemiological Models](@id epidemiology_stratification)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/covid/stratification.ipynb)

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab
using DisplayAs, Markdown

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

# This tutorial describes typed Petri nets and model stratification. 
# Model stratification is the process by which a submodel is replicated along some axis and 
# the submodels are allowed to interact in some specified way (e.g.; age-structured models).
# Specifically, methods from the section "Type systems for open Petri nets" of
# [[Libkind 2022](https://doi.org/10.1098/rsta.2021.0309)] are presented, which should
# be consulted for more information.

# ## Typed Petri nets

# An important part of model stratification is making sure that when combining e.g.; a mode of ageing
# and a model of disease transmission, that spurious or erroneous states or transitions are not generated
# (for example, simultaneously ageing and becoming infected). Category theory provides a way to
# specify formal stratification schemes that will respect domain knowledge via a type system for Petri nets.
# A typed Petri net is a *morphism* (generalized function) from a Petri net to another, where the codomain
# Petri net represents the type system under consideration, given as ``\phi: P \to P_{type}``.
# The morphism does not need to be injective, that is, multiple places can be mapped to the same place (likewise for transitions, etc).

# We denote the Petri net that defines a type system as ``P_{type}``. Let's consider a very simple type
# system that will nontheless be very useful for a class of models of directly transmitted diseases, called
# ``P_{infectious}``. There is a single type of place (population) and three types of transitions. Infectious
# transitions require two inputs and produce two outputs. Disease transitions are those that move a single
# individual amongst disease classes (e.g.; from ``E`` to ``I``). Strata transitions are those that move
# an individual amongst classes of the stratifying variable (e.g.; between age classes). Because all transitions
# in the generated (stratified Petri nets) must have a map into ``P_{infectious}``, we already see that spurious
# transitions like our example earlier are excluded.

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
  :disease => (:Pop => :Pop),
  :strata => (:Pop => :Pop)
)

to_graphviz(infectious_ontology)

# Now that we have ``P_{infectious}``, we define a SIRD model with that type system
# We start by defining an undirected wiring diagram (UWD) which represents the processes
# and states in the SIRD model; note that the name of each box corresponds to one of the transitions in ``P_{infectious}``,
# and that we specify the types of the junctions as `Pop` in the `where` clause.
# 
# Next we use `oapply_typed` which produces a typed Petri net by composing transitions based on ``P_{infectious}``. The
# first two arguments are the type Petri net and UWD, and the vector of symbols gives the names of the transitions in the
# new typed Petri net.
# We then add reflexive transitions typed as `:strata` to indicate which states can be stratified
# Here we add reflexive transitions to the susceptible, infected, and recovered populations but we leave out the dead
# population as no individuals entering that compartment may leave.
# 
# Calling `to_graphviz` on the resulting object will display $\phi : P \to P_{infectious}$. We could apply `dom` (domain)
# or `codom` (codomain) on the object to extract $P$ or $P_{infectious}$, respectively.

sird_uwd = @relation (S,I,R,D) where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S, I, I, I)
  disease(I, R)
  disease(I, D)
end

sird_model = oapply_typed(infectious_ontology, sird_uwd, [:infection, :recovery, :death])
sird_model = add_reflexives(sird_model, [[:strata], [:strata], [:strata], []], infectious_ontology)

to_graphviz(sird_model)

# ## Define intervention models

# ### Masking model

# Let's say we want to stratify the populations in the SIRD model by masking or non-masking; or, generally
# by any reversible behavioral change.
# To generate our stratification scheme, we begin with a UWD describing the strata levels we are
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

# To generate the stratified SIRD model, we use `typed_product`. As described in [[Libkind 2022](https://doi.org/10.1098/rsta.2021.0309)],
# stratification of a model can be seen as a pullback of ``\phi`` and ``\phi^{'}``, or as a product in the slice category ``Petri/P_{type}``.

typed_product(sird_model, mask_model) |> dom |> to_graphviz

# ### Vaccine model

# By changing ``P^{'}`` slightly, we can generate a stratification scheme to model vaccination. In this case,
# the transition between strata is irreversible (non-leaky vaccine), infection may occur between any
# pair of individuals (imperfect protection), and change in disease state may occur in any strata.

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

# Once again, the typed product, or pullback gives the SIRD model stratified by vaccination.

typed_product(sird_model, vax_model) |> dom |> to_graphviz

# ### Mask-Vax Model

# By using a more complicated UWD, we can have a model that combines masking and vaccination.

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

# ### Flux model between $N$ regions

# In many cases, it is not practical to construct by hand the composition diagram used to build a complex typed Petri net. 
# Here we demonstrate how to use the imperative interface provided by Catlab to construct a UWD describing a model of a population spread out amongst $N$
# regions. Individuals can travel between regions, but disease transmission is only possible between individuals in the same region.
# This is known as a flux model because individuals have the same behavior, so that the model can correctly simulate
# net "flux" of individuals between regions but not specific travel patterns.

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

# As explored in [[Citron 2021](https://doi.org/10.1073/pnas.2007488118)], the flux model is not sufficient when individuals
# travel behavior differs by their home region. A simple extension is the "simple trip model". We build this model in two parts. First we develop a "living model"
# which classifies populations according to what place their home is. The function `pairwise_id_typed_petri` is used to
# generate the "infection" transitions between each pair of populations. Reflexive transitions are added for strata and
# disease state transitions. 

function living_model(n)
  typed_living = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, [Symbol("Living$(i)") for i in 1:n])
  add_reflexives(typed_living, repeat([[:disease, :strata]], n), infectious_ontology)
end

to_graphviz(dom(living_model(2)))

# Taking the typed product with the previous travel model will now stratify that model by living status.

simple_trip_model = typed_product(travel_model(2), living_model(2))
to_graphviz(dom(simple_trip_model))

# Stratify our SIRD model on this simple trip model between two regions:

typed_product(sird_model, simple_trip_model) |> dom |> to_graphviz

# ## Model stratification and performance

# Here we investigate performance of stratified models. We first set up a simple helper function to connect the undirected wiring diagrams to our
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
