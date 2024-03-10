# # [Lotka-Volterra Model](@id predation_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/predation/lotka-volterra.ipynb)

using AlgebraicPetri

using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Programs.RelationalPrograms

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

# #### Introduction

# In this tutorial we use the [Lotka-Volterra model](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations), a classic
# mathematical model from ecology, to demonstrate basic conceps of compositional model building using AlgebraicPetri.

# #### Step 1: Define the building block Petri nets needed to construct the model.
# 
# These are the basic Petri nets which will be substituted as concrete models in the compositional syntax.
# There are three processes which we need to represent explicitly: birth, predation, and death.
# Because these models will be composed, we use `Open` to generate an open Petri net. When given only a single argument,
# `Open` generates a structured multicospan, which is an object containing the original Petri net ``S`` (accessed by `apex`),
# a list of finite sets ``A_{1},\dots,A_{n}`` (accessed by `feet`), and morphisms ``A_{1}\to S,\dots,A_{n}\to S`` (accessed by `legs`).
# The feet are where the Petri net may interact with other systems, justifying the moniker open Petri nets.
# 
# The basic Petri nets are generated and the apex of each open Petri net is displayed.
birth_petri = Open(PetriNet(1, 1=>(1,1)));
to_graphviz(birth_petri)
#-
predation_petri = Open(PetriNet(2, (1,2)=>(2,2)));
to_graphviz(predation_petri)
#-
death_petri = Open(PetriNet(1, 1=>()));
to_graphviz(death_petri)

# #### Step 2: Specify interaction between submodels using a relational syntax

# To stitch independent open Petri nets together into a combined model, we need to specify an interaction pattern (composition syntax).
# This is specified as an undirected wiring diagram (UWD) using the `@relation` macro from Catlab.jl, please see the documentation
# there for more information.
# 
# The UWD produced here consists of 3 "boxes" and 2 "junctions". Boxes connect to junctions via wires based on the function-like syntax
# in the macro. Wires that lead off the display are "outer ports" indicating how output may be collected from the combined system.
# 
# We use `oapply` to substitute the open Petri nets constructed earlier into each box of the UWD, where the feet of each open Petri net
# are identified with the wires from each box to shared junctions. The result of `oapply` is, once again, a multicospan, where the legs
# now give the map between outer ports of the composed model to places. 

lotka_volterra = @relation (wolves, rabbits) begin
  birth(rabbits)
  predation(rabbits, wolves)
  death(wolves)
end
display_uwd(lotka_volterra)
#-
lv_dict = Dict(:birth=>birth_petri, :predation=>predation_petri, :death=>death_petri);
lotka_petri = apex(oapply(lotka_volterra, lv_dict))
to_graphviz(lotka_petri)

# We can now define an initial state `u0`, reaction rate constants for each transition `p`, and use the method `vectorfield`
# to create a function we can pass to a differential equation solver to simulate and plot a trajectory of our
# composed Lotka-Volterra model.

u0 = [100, 10];
p = [.3, .015, .7];
prob = ODEProblem(vectorfield(lotka_petri),u0,(0.0,100.0),p);
sol = solve(prob,Tsit5(),abstol=1e-8);
plot(sol, labels=["Rabbits" "Wolves"])

# #### Step 3: Extend your model to handle more complex phenomena

# By making a larger UWD to describe, we can describe a small food chain between little fish, big fish, and sharks.

dual_lv = @relation (fish, Fish, Shark) begin
  birth(fish)
  predation(fish, Fish)
  death(Fish)
  predation(Fish, Shark)
  death(Shark)
end
display_uwd(dual_lv)
#-
dual_lv_petri = apex(oapply(dual_lv, lv_dict))
to_graphviz(dual_lv_petri)

# Generate a new solver, provide parameters, and analyze results

u0 = [100, 10, 2];
p = [.3, .015, .7, .017, .35];
prob = ODEProblem(vectorfield(dual_lv_petri),u0,(0.0,100.0),p);
sol = solve(prob,Tsit5(),abstol=1e-6);
plot(sol, label=["Little fish" "Big fish" "Sharks"])
