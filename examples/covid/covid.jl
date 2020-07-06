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

import Catlab.Theories: id

# A few helper functions
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true)
id(args...) = foldl((x,y)->id(x) ⊗ id(y), args);
# -

# #### Step 1: Define building block Petri Net models

ob = PetriCospanOb(1);

spontaneous_petri = PetriCospan([1], Petri.Model(1:2, [(Dict(1=>1), Dict(2=>1))]), [2]);

transmission_petri = PetriCospan([1], Petri.Model(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))]), [2]);

exposure_petri = PetriCospan([1, 2], Petri.Model(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]), [3, 2]);
travel_petri = PetriCospan(
        [1,2,3],
        Petri.Model(1:6, [(Dict(1=>1), Dict(4=>1)), (Dict(2=>1), Dict(5=>1)), (Dict(3=>1), Dict(6=>1))]),
        [4,5,6]);

# #### Step 2: Define a strongly type presentation of the Free Biproduct Category for the desired domain

# +
@present EpidemiologyWithTravel(FreeBiproductCategory) begin
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
    travel::Hom(S⊗E⊗I,S⊗E⊗I)
end

# Create the generators
S,E,I,R,D,transmission,exposure,illness,recovery,death,travel = generators(EpidemiologyWithTravel)

# Define a functor from the generators to the building block Petri Nets
F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
        S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri,travel=>travel_petri));
# -

# #### Step 3: Create, visualize, and solve model

# ### COVID-19 Multi-City MODEL:
#
# SEIRD City Model with travel between cities as $S ⊗ E ⊗ I → S ⊗ E ⊗ I$
#
# $seird_{city} = (((Δ(S) ⊗ id(E)) ⋅ (id(S) ⊗ σ(S,E))) ⊗ id(I)) ⋅ (id(S, E) ⊗ exposure) ⋅ (id(S) ⊗ (∇(E) ⋅ Δ(E)) ⊗ id(I)) ⋅ (id(S, E) ⊗ ((illness ⊗ id(I)) ⋅ (∇(I) ⋅ Δ(I)) ⋅ (id(I) ⊗ (Δ(I) ⋅(recovery ⊗ death))))) ⋅ (travel ⊗ ◊(R) ⊗ ◊(D))$
#       
# This is a very complicated interaction, so we can use the program interface for easier model definition

# +
seird_city = @program EpidemiologyWithTravel (s::S, e::E, i::I) begin
    e2, i2 = exposure(s, i)
    i3 = illness(e2)
    d = death(i3)
    r = recovery(i3)
    return travel(s, [e, e2], [i2, i3])
end
seird_city = to_hom_expr(FreeBiproductCategory, seird_city)

display_wd(seird_city)
# -

Graph(decoration(F(seird_city)))

# +
# function to compose n city models together
ncities(city,n::Int) = compose([city for i in 1:n]...)

# create a 3 city SEIRD models
seird_3 = ncities(seird_city, 3)
pc_seird_3 = F(seird_3)
p_seird_3 = decoration(pc_seird_3)
display(display_wd(seird_3))
display(Graph(p_seird_3))
# -

# Define time frame, 2 months
tspan = (0.0,60.0)
# Define initial states
u0 = zeros(Float64, base(pc_seird_3).n)
u0[1]  = 10000
u0[6]  = 10000
u0[11] = 10000
u0[2]  = 1
# Define transition rates
seirdparams(n::Int, k::Number) = begin
    βseird = [10/sum(u0), 1/2, 1/5, 1/16]
    βtravel = [1/2, 1/200, 1/2]/100k
    β = vcat(βseird, βtravel)
    return foldl(vcat, [β for i in 1:n])
end
params = seirdparams(3, 2)
# Generate and solve resulting ODE
prob = ODEProblem(toODE(p_seird_3),u0,tspan,params)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)
