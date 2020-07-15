# # [Basic Epidemiology Models](@id epidemiology_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/epidemiology.ipynb)

using AlgebraicPetri.Epidemiology

using Petri
using StochasticDiffEq
using Plots

using Catlab.Theories
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.Graphics

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# #### SIR Model:

# define model

sir = transmission ⋅ recovery

# get resulting petri net and visualize model

p_sir = decoration(F_epi(sir));
display_wd(sir)
#-
Graph(p_sir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0];
p = [0.4, 0.4];

prob,cb = SDEProblem(p_sir,u0,(0.0,7.5),p);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)

# #### SEIR Model:

# define model

seir = exposure ⋅ illness ⋅ recovery

# get resulting petri net and visualize model

p_seir = decoration(F_epi(seir));

display_wd(seir)
#-
Graph(p_seir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0, 0];
p = [0.9, 0.2, 0.5];

prob,cb = SDEProblem(p_seir,u0,(0.0,15.0),p);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)

# #### SEIRD Model:

# define model

seird = exposure ⋅ illness ⋅ Δ(I) ⋅ (death ⊗ recovery)

# get resulting petri net and visualize model

p_seird = decoration(F_epi(seird));

display_wd(seird)
#-
Graph(p_seird)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0, 0, 0];
p = [0.9, 0.2, 0.5, 0.1];

prob,cb = SDEProblem(p_seird,u0,(0.0,15.0),p);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)