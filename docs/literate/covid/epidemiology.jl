# # [Basic Epidemiology Models](@id epidemiology_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/covid/epidemiology.ipynb)

using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Programs.RelationalPrograms

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

# #### SIR Model

# We use the `@relation` macro from Catlab to create an undirected wiring diagram (UWD)
# which describes the composition syntax for the SIR model. Briefly, boxes (labeled ovals)
# represent processes which may depend on (consume or produce) resources represented by
# junctions (labeled dots).
sir = @relation (s,i,r) begin
    infection(s,i)
    recovery(i,r)
end
display_uwd(sir)

# To generate a Petri net (PN) from our compositional syntax, we need to apply a composition
# syntax, which assigns concrete mathematical models to each process. To compose
# with PNs, each box in the UWD will correspond to an open PN whose feet attach to
# the junctions that box is connected to. The composite PN is then constructed by
# gluing the component PNs along at the shared junctions (places).
# 
# The method `oapply_epi` is used to do this composition. It is a wrapper for the `oapply`
# method, which can be used with general models, and substitutes PNs for boxes
# in the UWD according to the box name. For a list of the allowed box labels/PNs,
# please see the help file for the function.
# 
# The returned object is a `MultiCospan`, where the feet are each outer port of the UWD
# and legs are morphisms which identify each outer port to a place in the composed PN.
# The apex of the multicospan is thus the composed model.
#-
p_sir = apex(oapply_epi(sir))
to_graphviz(p_sir)

# Define initial states and transition rates, then
# create, solve, and visualize ODE problem.

u0 = LVector(S=10, I=1, R=0);
p = LVector(inf=0.4, rec=0.4);

# The `vectorfield` method interprets the PN as describing mass-action kinetics
# with a rate constant associated to each transition, which can be used to
# simulate ODEs associated to the PN.

prob = ODEProblem(vectorfield(p_sir),u0,(0.0,7.5),p);
sol = solve(prob,Tsit5())

plot(sol, labels=hcat(string.(p_sir[:,:sname])...))

# #### SEIR Model:

# define model
seir = @relation (s,e,i,r) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
end
display_uwd(seir)
#-
p_seir = apex(oapply_epi(seir))
to_graphviz(p_seir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = LVector(S=10, E=1, I=0, R=0);
p = LVector(exp=.9, ill=.2, rec=.5);

prob = ODEProblem(vectorfield(p_seir),u0,(0.0,15.0),p);
sol = solve(prob,Tsit5())

plot(sol)

# #### SEIRD Model:

# define model
seird = @relation (s,e,i,r,d) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
    death(i,d)
end
display_uwd(seird)
#-
p_seird = apex(oapply_epi(seird))
to_graphviz(p_seird)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = LVector(S=10, E=1, I=0, R=0, D=0);
p = LVector(exp=0.9, ill=0.2, rec=0.5, death=0.1);

prob = ODEProblem(vectorfield(p_seird),u0,(0.0,15.0),p);
sol = solve(prob,Tsit5())

plot(sol)
