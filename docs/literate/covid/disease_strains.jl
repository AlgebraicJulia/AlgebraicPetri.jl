# # [Multiple-Strain Model of COVID and Vaccines](@id multistrain_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/covid/disease_strains.ipynb)

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using DisplayAs, Markdown

# This example presents models incorporating multiple strains of disease and vaccine type.
# Importantly, it shows why stratification by strain is different from other stratifications, e.g. geography, age, or other strata, and requires using a different type system.

# ## Define basic epidemiology model

# We start by defining our basic type system for infectious disease models

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
  :disease => (:Pop => :Pop),
  :strata => (:Pop => :Pop)
)

to_graphviz(infectious_ontology)

# We define a simple SIRD model with reflexive transitions typed as `:strata` to indicate which states can be stratified
# Here we add reflexive transitions to the susceptible, infected, and recovered populations but we leave out the dead
# population because they cannot do things such as get vaccinated or travel between regions.

sird_uwd = @relation () where {(S::Pop, I::Pop, R::Pop, D::Pop)} begin
  infect(S, I, I, I)
  disease(I, R)
  disease(I, D)
end

sird_model = oapply_typed(infectious_ontology, sird_uwd, [:infection, :recovery, :death])
sird_model = add_reflexives(sird_model, [[:strata], [:strata], [:strata], []], infectious_ontology)

to_graphviz(dom(sird_model))

# ## Define a model of multiple vaccine types

# We define a model of vaccination with multiple vaccine types. 
# In this model, vaccination transitions are typed as `:strata`.
# Note that the `:infect` transitions must be included to enable cross-infection between different vax types.

function vax_model(n)
  uwd = RelationalPrograms.TypedUnnamedRelationDiagram{Symbol,Symbol,Symbol}()
  junction = :Unvaxxed
  junctions = Dict(junction => add_junction!(uwd, :Pop, variable=junction))
  for i in 1:n  
    junction = Symbol("VaxType$(i)")
    junctions[junction] = add_junction!(uwd, :Pop, variable=junction)
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

# We can now stratify the two typed models. 

typed_product(sird_model, vax_model(2)) |> dom |> to_graphviz

# ## Define a multi-strain model

# Here we define a model of multiple strains of disease infection.
# In this model, uninfected individuals can become infected by indivuals carrying one of the strains.
# This model can be typed  the `infectious_ontology`.
# We add reflexives of `:disease` and `:strata` for the strain states but only `:strata` for the uninfected state.

function strain_model′(n)
  uwd = RelationalPrograms.TypedUnnamedRelationDiagram{Symbol,Symbol,Symbol}()
  junction = :Uninfected
  junctions = Dict(junction => add_junction!(uwd, :Pop, variable=junction))
  for i in 1:n  
    junction = Symbol("Strain$(i)")
    junctions[junction] = add_junction!(uwd, :Pop, variable=junction)
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

typed_product(sird_model, strain_model′(2)) |> dom |> to_graphviz

# Unfortunately, stratification of these models does not produce the desired result. 
# There are quite a few extraneous states and transitions. 
# The primary issue is the asymmetry in the role of the uninfected population. 
# We can address this be changing the type system.

# ## Define a new type system and corresponding epidemiology and multi-strain models

# The new type system has separate states for uninfected and infected to account for the asymmetry in their role in infection
# and in the corresponding asymmetry in the intended stratification. 
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

# We now reform the SIRD and multi-strain models using the new type system 
sird_for_strains_uwd = @relation () where {(S::Uninf, I::Inf, R::Inf, D::Inf)} begin
  infect(S, I, I, I)
  disease(I, R)
  disease(I, D)
end
sird_for_strains_model = oapply_typed(strain_ontology, sird_for_strains_uwd, [:infection, :recovery, :death])
sird_for_strains_model = add_reflexives(sird_for_strains_model, [[:strataU], [:strataI], [:strataI], []], strain_ontology)

to_graphviz(dom(sird_for_strains_model))

function strain_model(n)
  uwd = RelationalPrograms.TypedUnnamedRelationDiagram{Symbol,Symbol,Symbol}()
  junction = :Uninfected
  junctions = Dict(junction => add_junction!(uwd, :Pop, variable=junction))
  for i in 1:n  
    junction = Symbol("Strain$(i)")
    junctions[junction] = add_junction!(uwd, :Inf, variable=junction)
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

# When we now stratify we get the desired model.

sird_strain = typed_product(sird_for_strains_model, strain_model(2)) |> dom |> to_graphviz

to_graphviz(dom(strain_model(2)))

# ## Post composition: Typing the type system

# In some instances we may want to relate models typed to different type systems. 
# For example, we usually type our `simple_trip` model of geographic regions to the `infectious_ontology` s.t. we can stratify a disease model by geographic regions, 
# but the multi-strain disease model above is typed by the new `strain_ontology`.

# Crucially, we can accomplish this IF there is an appropriate morphism (map) between the type systems because post-composition by a morphism of type systems is functorial.
# In this case there is a morphism from `strain_ontology` to `infectious_ontology`, so we can form the morphism

# ### Morphism from `strain_ontology` to `infectious_ontology`

strain_ont_uwd = @relation () where {(Uninf::Pop, Inf::Pop)} begin
  infect(Uninf, Inf, Inf, Inf)
  disease(Inf, Inf)
  strata(Inf, Inf)
  strata(Uninf, Uninf)
end
strain_ont_act = oapply_typed(infectious_ontology,strain_ont_uwd,[:infect,:disease,:strataI,:strataU])

# ### Define geographic models

# ### Travel model between $N$ regions

# For this model we can use a julia function to programmatically build up our undirected wiring diagram for defining this model.
# Here we want there to be $N$ regions in which people can travel between each region and people within the same region are able
# to infect other people in the same region.

function travel_model(n)
  uwd = RelationalPrograms.TypedUnnamedRelationDiagram{Symbol,Symbol,Symbol}()
  junctions = Dict(begin
    junction = Symbol("Region$(i)")
    junction => add_junction!(uwd, :Pop, variable=junction)
  end for i in 1:n)

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


test = typed_product(sird_for_strains_model, strain_model(2))
test2 = typed_product(compose(test,strain_ont_act),simple_trip_model)

