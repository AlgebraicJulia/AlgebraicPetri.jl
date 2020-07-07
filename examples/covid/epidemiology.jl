# # [Basic Epidemiology Models](@id epidemiology_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/epidemiology.ipynb)

using AlgebraicPetri

using Petri
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Theories
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.Graphics
using Catlab.Graphics.Graphviz: Graph

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# ### Step 1: Define building block Petri Net models

ob = PetriCospanOb(1);
Graph(decoration(id(ob)))
#-
spontaneous_petri = PetriCospan([1], Petri.Model(1:2, [(Dict(1=>1), Dict(2=>1))]), [2]);
Graph(decoration(spontaneous_petri))
#-
transmission_petri = PetriCospan([1], Petri.Model(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))]), [2]);
Graph(decoration(transmission_petri))
#-
exposure_petri = PetriCospan([1, 2], Petri.Model(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]), [3, 2]);
Graph(decoration(exposure_petri))

# ### Step 2: Define a strongly type presentation of the Free Biproduct Category for the desired domain

@present Epidemiology(FreeBiproductCategory) begin
    S::Ob
    E::Ob
    I::Ob
    R::Ob
    D::Ob
    transmission::Hom(S⊗I, I)
    exposure::Hom(S⊗I, E⊗I)
    illness::Hom(E,I)
    recovery::Hom(I,R)
    death::Hom(I,D)
end;

# Create the generators
S,E,I,R,D,transmission,exposure,illness,recovery,death = generators(Epidemiology);

# Define a functor from the generators to the building block Petri Nets
F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
        S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri));

# ### Step 3: Create, visualize, and solve possible models

# #### SIR Model:

# define model

sir = transmission ⋅ recovery

# get resulting petri net and visualize model

p_sir = decoration(F(sir));
display_wd(sir)
#-
Graph(p_sir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0];
p = [0.4, 0.4];

prob = ODEProblem(toODE(p_sir),u0,(0.0,7.5),p);
sol = OrdinaryDiffEq.solve(prob,Tsit5());

plot(sol)

# #### SEIR Model:

# define model

sei = exposure ⋅ (illness ⊗ id(I)) ⋅ ∇(I)

seir = sei ⋅ recovery

# get resulting petri net and visualize model

p_seir = decoration(F(seir));

display_wd(seir)
#-
Graph(p_seir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0, 0];
p = [0.9, 0.2, 0.5];

prob = ODEProblem(toODE(p_seir),u0,(0.0,15.0),p);
sol = OrdinaryDiffEq.solve(prob,Tsit5());

plot(sol)

# #### SEIRD Model:

# define model

seird = sei ⋅ Δ(I) ⋅ (death ⊗ recovery)

# get resulting petri net and visualize model

p_seird = decoration(F(seird));

display_wd(seird)
#-
Graph(p_seird)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1, 0, 0, 0];
p = [0.9, 0.2, 0.5, 0.1];

prob = ODEProblem(toODE(p_seird),u0,(0.0,15.0),p);
sol = OrdinaryDiffEq.solve(prob,Tsit5());

plot(sol)
