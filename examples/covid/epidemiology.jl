# # [Basic Epidemiology Models](@id epidemiology_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/epidemiology.ipynb)

using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using Petri
using StochasticDiffEq
using Plots

using Catlab
using Catlab.Programs
using Catlab.Theories
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.Graphics

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# #### SIR Model:

# define model

sir = transmission ⋅ recovery

# get resulting petri net and visualize model

pob = PetriCospanOb(1)
spontaneous_petri(rate, u0) = PetriCospan([1], PetriWithRates(1:2, [(Dict(1=>1), Dict(2=>1))], [rate], u0), [2])
transmission_petri(rate, u0) = PetriCospan([1,2], PetriWithRates(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))], [rate], u0), [2])
exposure_petri(rate, u0) = PetriCospan([1, 2], PetriWithRates(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))], [rate], u0), [3])

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(S=>pob, E=>pob, I=>pob, R=>pob, D=>pob,
        transmission=>transmission_petri(0.4, [10,1]), exposure=>exposure_petri(.9, [10,0,1]),
        illness=>spontaneous_petri(.2, [0,0]), recovery=>spontaneous_petri(0.4,[0,0]), death=>spontaneous_petri(0.1,[0,0])))

p_sir = decoration(F(sir));
display_wd(sir)
#-
Graph(p_sir.m)

u0 = [10.0, 1, 0];
p = [0.4, 0.4];

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

prob,cb = SDEProblem(p_sir.m,p_sir.u0,(0.0,7.5),p_sir.rates);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)

# #### SEIR Model:

# define model

seir = @program InfectiousDiseases (s::S,i::I) begin
  e = exposure(s,i)
  i2 = illness(e)
  i_all = [i,i2]
  return recovery(i_all)
end
seir = to_hom_expr(FreeBiproductCategory, seir)

# get resulting petri net and visualize model

p_seir = decoration(F(seir));

display_wd(seir)
#-
Graph(p_seir.m)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

prob,cb = SDEProblem(p_seir.m,p_seir.u0,(0.0,15.0),p_seir.rates);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)

# #### SEIRD Model:

# define model

seird = @program InfectiousDiseases (s::S,i::I) begin
  e = exposure(s,i)
  i2 = illness(e)
  i_all = [i,i2]
  return recovery(i_all), death(i_all)
end
seird = to_hom_expr(FreeBiproductCategory, seird)

# get resulting petri net and visualize model

p_seird = decoration(F(seird));

display_wd(seird)
#-
Graph(p_seird)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

prob,cb = SDEProblem(p_seird.m,p_seird.u0,(0.0,15.0),p_seird.rates);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)