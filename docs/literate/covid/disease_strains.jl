# # [Multiple-Strain Model of COVID and Vaccines](@id multistrain_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/covid/disease_strains.ipynb)

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using DisplayAs, Markdown

# This example presents models incorporating multiple strains of disease and vaccine type.
# Importantly, it shows why stratification by disease strain is different from other stratifications, e.g. geography or age, and requires using a different type system.
# If you are unfamiliar with stratification, we reccomend first reading the [stratification tutorial](@ref stratification_example).

# ## Define basic epidemiology model

# We start by defining our basic type system, $P_{infectious}$ for infectious disease models.

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
  :disease => (:Pop => :Pop),
  :strata => (:Pop => :Pop)
)

to_graphviz(infectious_ontology)

# We define a simple SIRD model with reflexive transitions typed as "strata" to indicate which stratified
# states will interact with transitions amongst strata. All but the dead population will have reflexive 
# transitions amongst strata because no individuals may leave any death state.

sird_uwd = @relation (S,I,R,D) where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S, I, I, I)
  disease(I, R)
  disease(I, D)
end

sird_model = oapply_typed(infectious_ontology, sird_uwd, [:infection, :recovery, :death])
sird_model = add_reflexives(sird_model, [[:strata], [:strata], [:strata], []], infectious_ontology)

to_graphviz(dom(sird_model))

# ## Define a model of multiple vaccine types

# We also define a model of vaccination with multiple vaccine types, typed over $P_{infectious}$.
# Each stratum has reflexive "disease" transitions. Transitions which represent the administration
# of a vaccine to an individual are typed as "strata". Infection (typed as "infect") between individuals of different
# vaccination status is possible (i.e.; we do not assume a perfect vaccine).

function vax_model(n)
  uwd = RelationDiagram(repeat([:Pop], n+1))

  variable = :Unvaxxed
  junction = add_junction!(uwd, :Pop, variable=variable)
  port = ports(uwd, outer=true)[1]
  set_junction!(uwd, port, junction, outer=true)
  junctions = Dict(variable => junction)
  for i in 1:n  
    variable = Symbol("VaxType$(i)")
    junction = add_junction!(uwd, :Pop, variable=variable)
    port = ports(uwd, outer=true)[i+1]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end

  strains = filter((x) -> x != Symbol("Unvaxxed"), keys(junctions))
  for s in strains 
    pair = (:Unvaxxed, s)
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:strata)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end
  tnames = [Symbol("vax_$(b)") for b in strains]

  pairs = collect(Iterators.product(keys(junctions), keys(junctions)))
  for pair in pairs
    ins_outs = (pair[1], pair[2], pair[1], pair[2])
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in ins_outs], name=:infect)
    for (rgn, port) in zip(ins_outs, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("inf_$(pair[1])_$(pair[2])"))
  end
  act = oapply_typed(infectious_ontology, uwd, tnames)
  add_reflexives(act, repeat([[:disease]], n+1), infectious_ontology)
end

to_graphviz(dom(vax_model(2)))

# ## Stratify the SIRD model with vaccinations for two vaccine types

# We can now stratify the two typed models to get a model of SIRD with vaccination by multiple possible vaccine types. 

typed_product(sird_model, vax_model(2)) |> dom |> to_graphviz

# ## Define a model of multiple disease strains

# Here we define a model of multiple strains of disease infection.
# In this model, uninfected individuals can become infected by indivuals carrying one of the strains.
# This model can be typed  the `infectious_ontology`.
# We add reflexives of `:disease` and `:strata` for the strain states but only `:strata` for the uninfected state.

function strain_model′(n)
  uwd = RelationDiagram(repeat([:Pop], n+1))
  variable = :Uninfected
  junction = add_junction!(uwd, :Pop, variable=variable)
  port = ports(uwd, outer=true)[1]
  set_junction!(uwd, port, junction, outer=true)
  junctions = Dict(variable => junction)
  for i in 1:n  
    variable = Symbol("Strain$(i)")
    junction = add_junction!(uwd, :Pop, variable=variable)
    port = ports(uwd, outer=true)[i+1]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end
  strains = filter((x) -> x != Symbol("Uninfected"), keys(junctions))
  for s in strains 
    pair = (:Uninfected, s, s, s)
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:infect)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end
  act = oapply_typed(infectious_ontology, uwd, [Symbol("inf_$(b)") for b in strains])
  add_reflexives(act, vcat([[:strata]],repeat([[:disease,:strata]], n)), infectious_ontology)
end

to_graphviz(dom(strain_model′(2)))

# ## Stratify the SIRD model for two strains

# Unfortunately, stratification of these models does not produce the desired result, leading
# to many states and transitions. The problem is that it does not make sense for the
# uninfected population to be stratified over the full set of states of the SIRD model
# (i.e.; uninfected persons can never have any disease state other than "S"). 
# We can address this by changing the type system.

typed_product(sird_model, strain_model′(2)) |> dom |> to_graphviz

# ## Define a new type system and corresponding disease and multi-strain models

# The new type system has separate states for uninfected and infected to account for the asymmetry in their role in infection
# and for the corresponding asymmetry in the intended stratification. 
# Accordingly, the `:Inf` state has `:disease` and `:strataI` transitions. 
# The `:Uninf` state only has an additional "strata" transition, `:strataU`, but note this transition is distinct from that for the `:Inf` state.

const strain_ontology = LabelledPetriNet(
  [:Uninf, :Inf],
  :infect => ((:Uninf, :Inf) => (:Inf, :Inf)),
  :disease => (:Inf => :Inf),
  :strataI => (:Inf => :Inf),
  :strataU => (:Uninf => :Uninf)
)

to_graphviz(strain_ontology)

# We now reform the SIRD model using the new type system. Note that the `where` clause is used to define
# the types for boxes in the UWD.
sird_for_strains_uwd = @relation (S,I,R,D) where (S::Uninf, I::Inf, R::Inf, D::Inf) begin
  infect(S, I, I, I)
  disease(I, R)
  disease(I, D)
end
sird_for_strains_model = oapply_typed(strain_ontology, sird_for_strains_uwd, [:infection, :recovery, :death])
sird_for_strains_model = add_reflexives(sird_for_strains_model, [[:strataU], [:strataI], [:strataI], []], strain_ontology)

to_graphviz(dom(sird_for_strains_model))

# And similarly reform the multi-strain model.

function strain_model(n)
  uwd = RelationDiagram(vcat([:Uninf],repeat([:Inf], n)))
  variable = :Uninfected
  junction = add_junction!(uwd, :Uninf, variable=variable)
  port = ports(uwd, outer=true)[1]
  set_junction!(uwd, port, junction, outer=true)
  junctions = Dict(variable => junction)
  for i in 1:n  
    variable = Symbol("Strain$(i)")
    junction = add_junction!(uwd, :Inf, variable=variable)
    port = ports(uwd, outer=true)[i+1]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end
  strains = filter((x) -> x != Symbol("Uninfected"), keys(junctions))
  for s in strains 
    pair = (:Uninfected, s, s, s)
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:infect)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end
  act = oapply_typed(strain_ontology, uwd, [Symbol("inf_$(b)") for b in strains])
  add_reflexives(act, vcat([[:strataU]],repeat([[:disease,:strataI]], n)), strain_ontology)
end

to_graphviz(dom(strain_model(2)))

# When we now stratify we get the desired model. Note however that there are extraneous "strata"
# transitions; this is because there were no non-trivial stratum transitions in either $\phi$ or $\phi^{'}$
# (note that only "disease" and "infect" transition types were needed to fully define the multiple strain dynamics).

sird_strain = typed_product(sird_for_strains_model, strain_model(2)) 

to_graphviz(dom(sird_strain))

# ## Post-composition: Retyping the type system

# In some instances, we may want to relate models typed to different type systems. 
# For example, we usually type our `simple_trip` model of geographic regions to the `infectious_ontology` such that we can stratify a disease model by geographic regions, 
# but the multi-strain disease model above is typed by the new `strain_ontology`.

# Crucially, we can accomplish this IF there is an appropriate morphism (map) between the type systems 
# ($P_{strain}$ and $P_{infectious}$) because post-composition by a morphism of type systems is functorial.
# In this case, there is a morphism from `strain_ontology` to `infectious_ontology`, so we can form the morphism

# ### Morphism from `strain_ontology` to `infectious_ontology`

# We use `oapply_typed` on a UWD representation of the `strain_ontology`, but note that we could also directly
# define the map $P_{strain}\to P_{infectious}$ with `ACSetTransformation`.

strain_ont_uwd = @relation (Uninf,Inf) where (Uninf::Pop, Inf::Pop) begin
  infect(Uninf, Inf, Inf, Inf)
  disease(Inf, Inf)
  strata(Inf, Inf)
  strata(Uninf, Uninf)
end
strain_ont_act = oapply_typed(infectious_ontology,strain_ont_uwd,[:infect,:disease,:strataI,:strataU])

to_graphviz(strain_ont_act)

# ### Define simple-trip geographic model of $N$ regions

# To demonstrate stratification utilizing post-composition to re-type the models, we use the simple-trip geographic model.
# This model is comprises a travel model and a living model.

# In this model, there are $N$ regions which people can travel between. People within the same region are able
# to infect other people in the same region. It is typed by `infectious_ontology` ($P_{infectious}$).

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

# This model could itself be stratified with the SIRD model, but we want to model
# persons travelling between locations while maintaining the status of where they live.  

# In this model, people have the property of "Living" somewhere.
# It is also typed by `infectious_ontology` ($P_{infectious}$).

function living_model(n)
  typed_living = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, [Symbol("Living$(i)") for i in 1:n])
  add_reflexives(typed_living, repeat([[:disease, :strata]], n), infectious_ontology)
end

to_graphviz(dom(living_model(2)))

# We can stratify the living model with the travel model to get a model of someone taking a simple trip.

simple_trip_model = typed_product(travel_model(2), living_model(2))

to_graphviz(dom(simple_trip_model))

# ### Stratify multi-strain SIRD and simple-trip models

# Now, to stratify our multi-strain SIRD model with the simple-trip, we first retype the multi-strain model 
# to the `infectious_ontology` by composing with the morphism we defined. Mathematically, because
# the multi-strain SIRD model is a typed Petri net, it is a morphism $\phi : P \to P_{strain}$.
# The object `strain_ont_act` we made earlier is a morphism between $P_{strain} \to P_{infectious}$,
# so when we post-compose, we get a new morphism between $P \to P_{infectious}$.

sird_strain_retyped = compose(sird_strain,strain_ont_act)

# We can now take the typed product to get the multi-strain SIRD model stratified over the simple trip model
# of movement between different regions.

sird_strain_trip = typed_product(sird_strain_retyped,simple_trip_model)

to_graphviz(dom(sird_strain_trip))

# ## Define a multi-strain SIRD model with vaccination by multiple vaccine types

# Because the multi-vaccine model was typed according to $P_{infectious}$, our retyped
# multi-strain SIRD can be stratified according to the multiple vaccine model in a
# similar way.

sird_strain_vax = typed_product(sird_strain_retyped,vax_model(2))

to_graphviz(dom(sird_strain_vax))

# ## Re-stratify the multi-strain multi-vax SIRD model with the simple trip model

# If we would like to re-stratify our SIRD-strain-vax model with the simple trip model, we again face a difficulty.
# Both the "vaccination" transitions of the first model and the "travel" transitions of the second 
# are currently typed to the `:strata` transition of the `infectious_ontology` type system.
# Naively stratifying would thus produce transitions in which persons traveled and were vaccinated simultaneously.

# To appropriately stratify, we need an additional "strata" transition to distinguish 
# between the two types of transitions. 
# We can again use post-compostion with the morphism between type systems to reuse our existing models.

# ### Define an augmented version of the `infectious_ontology` type system with an additional "strata" transition

const aug_inf_ontology = LabelledPetriNet(
  [:nPop],
  :ninfect => ((:nPop, :nPop) => (:nPop, :nPop)),
  :ndisease => (:nPop => :nPop),
  :nstrata => (:nPop => :nPop),
  :nstrata2 => (:nPop => :nPop)
)

to_graphviz(aug_inf_ontology)

# ### Define morphisms from the original type system to the new augmented type system

# We form one morphism that maps the `:strata` transition to `:nstrata`.
# This morphism will serve to re-type the SIRD-strain-vax model.

function retype_inf_ont(strata_map)
  uwd = RelationDiagram([:nPop])
  variable = :Pop
  junction = add_junction!(uwd, :nPop, variable=variable)
  port = ports(uwd, outer=true)[1]
  set_junction!(uwd, port, junction, outer=true)
  junctions = Dict(variable => junction)

  boxes = [:ninfect, :ndisease, strata_map]
  for bname in boxes
    if bname == :ninfect
      pair = (:Pop, :Pop, :Pop, :Pop)
    else 
      pair = (:Pop, :Pop)
    end
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=bname)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end

  act = oapply_typed(aug_inf_ontology, uwd, [:infect, :disease, :strata])
end

inf_ont_act = retype_inf_ont(:nstrata)

to_graphviz(inf_ont_act)

# We form another morphism that maps the `:strata` transition to `:nstrata2`.
# This morphism will serve to re-type the simple trip model.

rgn_ont_act = retype_inf_ont(:nstrata2)

to_graphviz(rgn_ont_act)

# ### Add reflexive transitions

# To finish preparing for stratification, we need to add the new reflexive transitions to the component models.
# To the SIRD-strain-vax model, we add an `:nstrata2` tranisiton to each state that does not represent 
# a portion of the population that is deceased (because deceased individuals cannot travel). 

sird_strain_vax_retyped = flatten_labels(compose(sird_strain_vax,inf_ont_act))
reflx = [[:nstrata2]]
for ii in 2:ns(dom(sird_strain_vax_retyped))
  if split(String(dom(sird_strain_vax_retyped)[ii,:sname]),"_")[1] == "D" 
    push!(reflx,[])
  else
    push!(reflx,[:nstrata2])  
  end
end
aug_sird_strain_vax = add_reflexives(sird_strain_vax_retyped, reflx, aug_inf_ontology);

# To the simple trip model, we add an `:nstrata` tranisiton for each state.

simple_trip_retyped = flatten_labels(compose(simple_trip_model,rgn_ont_act))
aug_trip = add_reflexives(simple_trip_retyped, repeat([[:nstrata]],ns(dom(simple_trip_retyped))), aug_inf_ontology);

# ### Stratify the SIRD-strain-vax and simple trip models

sird_strain_vax_trip = typed_product([aug_sird_strain_vax,aug_trip])

to_graphviz(dom(sird_strain_vax_trip))
