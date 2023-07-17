# # [Model of COVID Incorporating Wastewater](@id wastewater_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/covid/wastewater.ipynb)

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.Programs, Catlab.Graphics
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using DisplayAs, Markdown

# This example presents models incorporating wastewater as well as various forms of 
# model modification and stratifications. 
# From a modeling perspective, it would likely be easier and more natural to
# include wastewater via an observation model, but for demonstration purposes
# we include it structurally.

# ## Define basic epidemiology model

# We start by defining our basic type system for infectious disease models.

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
  :disease => (:Pop => :Pop),
  :strata => (:Pop => :Pop)
)

to_graphviz(infectious_ontology)

# We define a simple SEIRD model with reflexive transitions typed as `:strata` to indicate which states can be stratified.
# Here we add reflexive transitions to the susceptible, exposed, infected, and recovered populations, 
# but we leave out the dead population because they cannot do things such as get vaccinated or travel between regions.

seird_uwd = @relation (S,I,E,R,D) where (S::Pop, I::Pop, E::Pop, R::Pop, D::Pop) begin
  infect(S, I, E, I)
  disease(E, I)
  disease(I, R)
  disease(I, D)
end

seird_model = oapply_typed(infectious_ontology, seird_uwd, [:infection, :progress, :recovery, :death])
seird_model = add_reflexives(seird_model, [[:strata], [:strata], [:strata], [:strata], []], infectious_ontology)

to_graphviz(dom(seird_model))

# ## Add masking to the SEIRD model via stratification

# We define a model of masking where masked individuals can be infected but do not infect

masking_uwd = @relation (M,Um) where (M::Pop, Um::Pop) begin
  strata(M, Um)
  strata(Um, M)
  infect(M, Um, M, Um)
  infect(Um, Um, Um, Um)
end
mask_model = oapply_typed(infectious_ontology, masking_uwd, [:unmask, :mask, :infect_um, :infect_uu])
mask_model = add_reflexives(mask_model, [[:disease], [:disease]], infectious_ontology)

to_graphviz(dom(mask_model))

# Stratify our SEIRD model on this masking model to get a model of SEIRD with masking:

typed_product(seird_model, mask_model) |> dom |> to_graphviz

# ## Add the possibility of reinfection to the SEIRD model

# We can add the possibility of of becoming reinfected to the model by adding another transition
# from the `R` state to the `S` state. 
# For convenience, we achieve this by respecifying the UWD.

seirds_uwd = @relation (S,I,E,R,D) where (S::Pop, I::Pop, E::Pop, R::Pop, D::Pop) begin
  infect(S, I, E, I)
  disease(E, I)
  disease(I, R)
  disease(I, D)
  disease(R, S)
end

seirds_model = oapply_typed(infectious_ontology, seirds_uwd, [:infection, :progress, :recovery, :death, :wane])
seirds_model = add_reflexives(seirds_model, [[:strata], [:strata], [:strata], [:strata], []], infectious_ontology)

to_graphviz(dom(seirds_model))

# We again stratify to add masking 

seirds_mask = typed_product(seirds_model, mask_model)

seirds_mask |> dom |> to_graphviz

# ## Incorporate wastewater into the disease model

# To incorporate wastewater we need to form a new type system to capture its formation.

const ww_ontology = LabelledPetriNet(
  [:Pop,:Food,:Waste],
  :infect => ((:Pop, :Pop) => (:Pop, :Pop)),
  :disease => (:Pop => :Pop),
  :strata => (:Pop => :Pop),
  :grow => (() => :Food),
  :dispose => (:Waste => ()),
  :digest => ((:Pop, :Food) => (:Pop, :Waste))
)

to_graphviz(ww_ontology)

# We can form a model of wastewater production

ww_uwd = @relation (P,F,W) where (P::Pop, F::Food, W::Waste) begin
  digest(P, F, P, W)
end

ww_model = oapply_typed(ww_ontology, ww_uwd, [:digest])
ww_model = add_reflexives(ww_model, [[:infect,:disease,:strata], [:grow], [:dispose]], ww_ontology)

to_graphviz(dom(ww_model))

# It is more convenient to add the new states and transitions to our SEIRDS and masking models directly.

function seirds_ww_model()
  uwd = RelationDiagram(vcat(repeat([:Pop], 5),repeat([:Food, :Waste], 4)))
  disease_states = ["S","I","E","R","D"]
  junctions = Dict(begin
    variable = Symbol(s)
    junction = add_junction!(uwd, :Pop, variable=variable)
    port = ports(uwd, outer=true)[i]
    set_junction!(uwd, port, junction, outer=true)
    variable => junction
  end for (i, s) in enumerate(disease_states))
  digest_states = ["S","I","E","R"] # omit dead state from digesting
  for (i, s) in enumerate(digest_states)
    variable = Symbol("F_"*s)
    junction = add_junction!(uwd, :Food, variable=variable)
    port = ports(uwd, outer=true)[1+(i-1)*2+5]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction

    variable = Symbol("W_"*s)
    junction = add_junction!(uwd, :Waste, variable=variable)
    port = ports(uwd, outer=true)[i*2+5]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end  

  pairs = [(:S,:I,:E,:I),(:E,:I),(:I,:R),(:I,:D),(:R,:S)]
  tnames = [:infection, :progress, :recovery, :death, :wane]
  for pair in pairs
    if length(pair)==4
      tmp_name = :infect
    else
      tmp_name = :disease
    end
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=tmp_name)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end

  for s in digest_states
    pair = (Symbol(s), Symbol("F_"*s), Symbol(s), Symbol("W_"*s))
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:digest)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("digest_"*s))
  end
  
  act = oapply_typed(ww_ontology, uwd, tnames)
  add_reflexives(act, vcat(repeat([[:strata]], 4), repeat([[]],9)), ww_ontology)
end

to_graphviz(dom(seirds_ww_model()))


function mask_ww_model()
  uwd = RelationDiagram(vcat(repeat([:Pop], 2),repeat([:Food, :Waste], 2)))
  disease_states = ["M","Um"]
  junctions = Dict(begin
    variable = Symbol(s)
    junction = add_junction!(uwd, :Pop, variable=variable)
    port = ports(uwd, outer=true)[i]
    set_junction!(uwd, port, junction, outer=true)
    variable => junction
  end for (i, s) in enumerate(disease_states))
  for (i, s) in enumerate(disease_states)
    variable = Symbol("F_"*s)
    junction = add_junction!(uwd, :Food, variable=variable)
    port = ports(uwd, outer=true)[1+2*i]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction

    variable = Symbol("W_"*s)
    junction = add_junction!(uwd, :Waste, variable=variable)
    port = ports(uwd, outer=true)[2+2*i]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end  

  pairs = [(:M,:Um),(:Um,:M),(:M,:Um,:M,:Um),(:Um,:Um,:Um,:Um)]
  tnames = [:unmask, :mask, :infect_um, :infect_uu]
  for pair in pairs
    if length(pair)==4
      tmp_name = :infect
    else
      tmp_name = :strata
    end
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=tmp_name)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end

  for s in disease_states
    pair = (Symbol(s), Symbol("F_"*s), Symbol(s), Symbol("W_"*s))
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:digest)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("digest_"*s))
  end
  
  act = oapply_typed(ww_ontology, uwd, tnames)
  add_reflexives(act, vcat(repeat([[:disease]], 2), repeat([[]],4)), ww_ontology)
end

to_graphviz(dom(mask_ww_model()))

# Now we can stratify the models with wastewater production

seirds_mask_ww = typed_product([seirds_ww_model(), mask_ww_model()])

seirds_mask_ww |> dom |> to_graphviz

# ## Compare the structure of the models thus far

# mca = maximum_common_subobject([dom(seird_model),dom(seirds_model),dom(seirds_mask),dom(seirds_mask_ww)])

# ## Explore further modifications and stratifications

# We now explore other model modifications and stratifications that may be of interest. 
# Because a full model incorporating all components will require forming a new type system and re-typing, 
# we first form the component models and demonstrate them by stratifying with our SEIRDS-WW model.

# ### Define an age model

# Here we set up a model with `n` age cohorts

function age_ww_model(n)
  uwd = RelationDiagram(vcat(repeat([:Pop], n),repeat([:Food, :Waste], n)))
  junctions = Dict(begin
    variable = Symbol("Age$(i)")
    junction = add_junction!(uwd, :Pop, variable=variable)
    port = ports(uwd, outer=true)[i]
    set_junction!(uwd, port, junction, outer=true)
    variable => junction
  end for i in 1:n)
  for i in 1:n
    variable = Symbol("F_age$(i)")
    junction = add_junction!(uwd, :Food, variable=variable)
    port = ports(uwd, outer=true)[1+(i-1)*2+n]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction

    variable = Symbol("W_age$(i)")
    junction = add_junction!(uwd, :Waste, variable=variable)
    port = ports(uwd, outer=true)[i*2+n]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end  

  combos = collect(Iterators.product(1:n, 1:n))
  tnames = Vector{Symbol}()
  for combo in combos
    pair = (Symbol("Age$(combo[1])"),Symbol("Age$(combo[2])"),Symbol("Age$(combo[1])"),Symbol("Age$(combo[2])"))
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:infect)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("infect_$(combo[2])$(combo[1])"))
  end

  for i in 1:n
    pair = (Symbol("Age$(i)"), Symbol("F_age$(i)"), Symbol("Age$(i)"), Symbol("W_age$(i)"))
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:digest)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("digest$(i)"))
  end
  
  act = oapply_typed(ww_ontology, uwd, tnames)
  add_reflexives(act, vcat(repeat([[:disease]], n), repeat([[]],2*n)), ww_ontology)
end

to_graphviz(dom(age_ww_model(2)))

# We can stratify this with our SEIRDS-WW model

seirds_age_ww = typed_product([seirds_ww_model(), age_ww_model(2)])

seirds_age_ww |> dom |> to_graphviz

# ### Define a model of multiple vaccine types

# We also define a model of vaccination with multiple vaccine types. 
# In this model, vaccination transitions are typed as `:strata`.
# Note that the `:infect` transitions must be included to enable cross-infection between different vax types.

function vax_ww_model(n)
  uwd = RelationDiagram(vcat(repeat([:Pop], n+1),repeat([:Food, :Waste], n+1)))

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

  vax_names = deepcopy(keys(junctions))

  for (i, v) in enumerate(vax_names)
    variable = Symbol("F_"*String(v))
    junction = add_junction!(uwd, :Food, variable=variable)
    port = ports(uwd, outer=true)[1+(i-1)*2+n+1]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction

    variable = Symbol("W_"*String(v))
    junction = add_junction!(uwd, :Waste, variable=variable)
    port = ports(uwd, outer=true)[i*2+n+1]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end  

  for (i, v) in enumerate(vax_names)
    pair = (Symbol(v), Symbol("F_"*String(v)), Symbol(v), Symbol("W_"*String(v)))
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:digest)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("digest_"*String(v)))
  end

  act = oapply_typed(ww_ontology, uwd, tnames)
  add_reflexives(act, repeat([[:disease]], n+1), ww_ontology)
end

to_graphviz(dom(vax_ww_model(1)))

# We can now stratify the SEIRDS with vaccination by multiple possible vaccine types. 

seirds_vax_ww = typed_product(seirds_ww_model(), vax_ww_model(1))

seirds_vax_ww |> dom |> to_graphviz

# ### Define a model incorporating diagnosis 

function seirdsd_ww_model()
  uwd = RelationDiagram(vcat(repeat([:Pop], 6),repeat([:Food, :Waste], 5)))
  disease_states = ["S","I","E","Diag","R","D"]
  junctions = Dict(begin
    variable = Symbol(s)
    junction = add_junction!(uwd, :Pop, variable=variable)
    port = ports(uwd, outer=true)[i]
    set_junction!(uwd, port, junction, outer=true)
    variable => junction
  end for (i, s) in enumerate(disease_states))
  digest_states = ["S","I","Diag","E","R"] # omit dead state from digesting
  for (i, s) in enumerate(digest_states)
    variable = Symbol("F_"*s)
    junction = add_junction!(uwd, :Food, variable=variable)
    port = ports(uwd, outer=true)[1+(i-1)*2+6]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction

    variable = Symbol("W_"*s)
    junction = add_junction!(uwd, :Waste, variable=variable)
    port = ports(uwd, outer=true)[i*2+6]
    set_junction!(uwd, port, junction, outer=true)
    junctions[variable] = junction
  end  

  pairs = [(:S,:I,:E,:I),(:E,:I),(:I,:Diag),(:I,:R),(:Diag,:R),(:I,:D),(:Diag,:D),(:R,:S)]
  tnames = [:infection, :progress, :diagnosis, :recovery, :diag_recov, :death, :diag_death, :wane]
  for pair in pairs
    if length(pair)==4
      tmp_name = :infect
    else
      tmp_name = :disease
    end
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=tmp_name)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end

  for s in digest_states
    pair = (Symbol(s), Symbol("F_"*s), Symbol(s), Symbol("W_"*s))
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=:digest)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
    push!(tnames,Symbol("digest_"*s))
  end
  
  act = oapply_typed(ww_ontology, uwd, tnames)
  add_reflexives(act, vcat(repeat([[:strata]], 5), repeat([[]],11)), ww_ontology)
end

to_graphviz(dom(seirdsd_ww_model()))

# ## Re-stratify the multi-strain multi-vax SIRD model with the simple trip model

# If we would like to re-stratify our SIRD-strain-vax model with the simple trip model, we again face a difficulty.
# Both the "vaccination" transitions of the first model and the "travel" transitions of the second 
# are currently typed to the `:strata` transition of the `infectious_ontology` type system.
# To appropriately stratify, we need an additional "strata" transition to distinguish 
# between the two types of transitions. 
# We can again use post-compostion to overcome the difficulty.

# ### Define an augmented version of the `infectious_ontology` type system with an additional "strata" transition

const aug_ww_ontology = LabelledPetriNet(
  [:nPop,:nFood,:nWaste],
  :ninfect => ((:nPop, :nPop) => (:nPop, :nPop)),
  :ndisease => (:nPop => :nPop),
  :nstrata_mask => (:nPop => :nPop),
  :nstrata_vax => (:nPop => :nPop),
  :nstrata_travel => (:nPop => :nPop),
  :ngrow => (() => :nFood),
  :ndispose => (:nWaste => ()),
  :ndigest => ((:nPop, :nFood) => (:nPop, :nWaste))
)

to_graphviz(aug_ww_ontology)

# ### Define morphisms from the original type system to the new augmented type system

# We first form a function to generate the several re-typing morphisms to the different 
# "strata" states for the models we have made thus far. 

function retype_ww_ont(strata_map)
  uwd = RelationDiagram([:nPop,:nFood,:nWaste])
  vs = [:Pop,:Food,:Waste]
  junctions = Dict( begin
    variable = vs[i]
    junction = add_junction!(uwd, Symbol("n"*String(vs[i])), variable=variable)
    set_junction!(uwd, port, junction, outer=true)
    variable => junction
  end for (i, port) in enumerate(ports(uwd, outer=true)))
  
  boxes = [:ninfect, :ndisease, :ndigest, strata_map]
  for bname in boxes
    if bname == :ninfect
      pair = (:Pop, :Pop, :Pop, :Pop)
    elseif bname == :ndigest
      pair = (:Pop, :Food, :Pop, :Waste)
    else 
      pair = (:Pop, :Pop)
    end
    box = add_box!(uwd, [junction_type(uwd, junctions[p]) for p in pair], name=bname)
    for (rgn, port) in zip(pair, ports(uwd, box))
      set_junction!(uwd, port, junctions[rgn])
    end
  end

  act = oapply_typed(aug_ww_ontology, uwd, [:infect, :disease, :digest, :strata])
end

# With this, we form one morphism that maps the `:strata` transition to `:nstrata_vax`.
# This morphism will serve to re-type the vax model.

ww_vax_ont_act = retype_ww_ont(:nstrata_vax)

to_graphviz(ww_vax_ont_act)

# We form another morphism that maps the `:strata` transition to `:nstrata_mask`.
# This morphism will serve to re-type the mask model.

ww_mask_ont_act = retype_ww_ont(:nstrata_mask)

to_graphviz(ww_mask_ont_act)

# Lastly, We can form a morphism that maps the `:strata` transition to `:nstrata_travel`.
# This morphism will serve to re-type a travel model.

ww_travel_ont_act = retype_ww_ont(:nstrata_travel)

to_graphviz(ww_travel_ont_act)

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




# ### Define simple-trip geographic model of $N$ regions

# The simple-trip geographic model comprises a travel model and a living model.

# **Travel model between $N$ regions**: 
# In this model, there are $N$ regions which people can travel between. People within the same region are able
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

# **Living model of $N$ regions**: 
# In this model, people have the property of "Living" somewhere. 

function living_model(n)
  typed_living = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, [Symbol("Living$(i)") for i in 1:n])
  add_reflexives(typed_living, repeat([[:disease, :strata]], n), infectious_ontology)
end

to_graphviz(dom(living_model(2)))

# **Simple trip model of $N$ regions**: 
# We can stratify the living model with the travel model to get a model of someone taking a simple trip.

simple_trip_model = typed_product(travel_model(2), living_model(2))

to_graphviz(dom(simple_trip_model))

# ### Stratify SIRD-multi-strain and simple-trip models

# Now, to stratify our multi-strain SIRD model with the simple-trip, we first retype the multi-strain model 
# to the `infectious_ontology` by composing with the morphism we defined.

sird_strain_retyped = compose(sird_strain,strain_ont_act)

# We can now stratify.

sird_strain_trip = typed_product(sird_strain_retyped,simple_trip_model)

to_graphviz(dom(sird_strain_trip))

# ## Define a multi-strain SIRD model with vaccination by multiple vaccine types

# We can similarly stratify the multi-strain SIRD model with the multi-vax model.

sird_strain_vax = typed_product(sird_strain_retyped,vax_model(2))

to_graphviz(dom(sird_strain_vax))
