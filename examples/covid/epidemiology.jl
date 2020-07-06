# -*- coding: utf-8 -*-
# + {}
using Petri
using OrdinaryDiffEq
using Plots
using AlgebraicPetri
using Catlab
using Catlab.Theories
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.Graphics.Graphviz: Graph

# A helper function
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);
# -

# ### Step 1: Define building block Petri Net models

ob = PetriCospanOb(1)
Graph(decoration(id(ob)))

# +
spontaneous_petri = PetriCospan(
        Cospan(FinOrdFunction([1], 2),
               FinOrdFunction([2], 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1), Dict(2=>1))]))

Graph(decoration(spontaneous_petri))

# +
transmission_petri = PetriCospan(
        Cospan(FinOrdFunction([1], 2),
               FinOrdFunction([2], 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1, 2=>1), Dict(2=>2))]))

Graph(decoration(transmission_petri))

# +
exposure_petri = PetriCospan(
        Cospan(FinOrdFunction([1, 2], 3),
               FinOrdFunction([3, 2], 3)
        ), id(PetriFunctor), Petri.Model([1, 2, 3], [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]))

Graph(decoration(exposure_petri))
# -

# ### Step 2: Define a strongly type presentation of the Free Biproduct Category for the desired domain

# +
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
end

# Create the generators
S,E,I,R,D,transmission,exposure,illness,recovery,death = generators(Epidemiology)

# Define a functor from the generators to the building block Petri Nets
F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
        S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri))
map(display, generators(Epidemiology));
# -

# ### Step 3: Create, visualize, and solve possible models

# #### SIR Model:

# define model
sir = transmission ⋅ recovery

# get resulting petri net
p_sir = decoration(F(sir))
# display wiring diagram and petri net visualization
display(display_wd(sir))
display(Graph(p_sir))

# define initial states and transition rates
u0 = [10.0, 1, 0]
p = [0.4, 0.4]
# create and solve ODE problem
prob = ODEProblem(toODE(p_sir),u0,(0.0,7.5),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)

# #### SEIR Model:

# define model
sei = exposure ⋅ (illness ⊗ id(I)) ⋅ ∇(I)

seir = sei ⋅ recovery

# get resulting petri net
p_seir = decoration(F(seir))
# display wiring diagram and petri net visualization
display(display_wd(seir))
display(Graph(p_seir))

# define initial states and transition rates
u0 = [10.0, 1, 0, 0]
p = [0.9, 0.2, 0.5]
# create and solve ODE problem
prob = ODEProblem(toODE(p_seir),u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)

# #### SEIRD Model:

# define model
seird = sei ⋅ Δ(I) ⋅ (death ⊗ recovery)

# get resulting petri net
p_seird = decoration(F(seird))
# display wiring diagram and petri net visualization
display(display_wd(seird))
display(Graph(p_seird))

# define initial states and transition rates
u0 = [10.0, 1, 0, 0, 0]
p = [0.9, 0.2, 0.5, 0.1]
# create and solve ODE problem
prob = ODEProblem(toODE(p_seird),u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)
