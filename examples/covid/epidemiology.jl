# # [Basic Epidemiology Models](@id epidemiology_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/epidemiology.ipynb)

using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using Petri: Model, Graph
using OrdinaryDiffEq
using Plots

using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# #### SIR Model:

# define model

sir = transmission ⋅ recovery

# get resulting petri net as a C-Set

cset_sir = apex(F_epi(sir));
display_wd(sir)
#-

# Use Petri.jl to visualize the C-Set

Graph(Model(cset_sir))

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0];
p = [0.4, 0.4];

# The C-Set representation has direct support for generating a DiffEq vector field

prob = ODEProblem(vectorfield(cset_sir),u0,(0.0,7.5),p);
sol = solve(prob,Tsit5())

plot(sol)

# #### SEIR Model:

# define model

seir = @program InfectiousDiseases (s::S,i::I) begin
    i2 = illness(exposure(s,i))
    return recovery([i,i2])
end
seir = to_hom_expr(FreeBiproductCategory, seir)

# here we convert the C-Set decoration to a Petri.jl model
# to use its StochasticDifferentialEquations support

p_seir = apex(F_epi(seir));

display_wd(seir)
#-
Graph(Model(p_seir))

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0, 0];
p = [.9, .2, .5];

prob = ODEProblem(vectorfield(p_seir),u0,(0.0,15.0),p);
sol = solve(prob,Tsit5())

plot(sol)

# #### SEIRD Model:

# define model

seird = @program InfectiousDiseases (s::S,i::I) begin
    i_all = [i, illness(exposure(s,i))]
    return recovery(i_all), death(i_all)
end
seird = to_hom_expr(FreeBiproductCategory, seird)

# get resulting petri net and visualize model

p_seird = apex(F_epi(seird));

display_wd(seird)
#-
Graph(Model(p_seird))

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0, 0, 0];
p = [0.9, 0.2, 0.5, 0.1];

prob = ODEProblem(vectorfield(p_seird),u0,(0.0,15.0),p);
sol = solve(prob,Tsit5())

plot(sol)
