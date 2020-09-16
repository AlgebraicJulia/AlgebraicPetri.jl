# # [Lotka-Volterra Model](@id predation_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/predation/lotka-volterra.ipynb)

using AlgebraicPetri

using Petri
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Theories
using Catlab.Programs
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.WiringDiagrams
using Catlab.Graphics

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# #### Step 1: Define the building block Petri nets needed to construct the model

birth_petri = PetriCospan([1], PetriNet(1, (1, (1,1))), [1]);
Graph(Petri.Model(decoration(birth_petri)))
#-
predation_petri = PetriCospan([1,2], PetriNet(2, ((1,2), (2,2))), [2]);
Graph(Petri.Model(decoration(predation_petri)))
#-
death_petri = PetriCospan([1], PetriNet(1, (1, ())), [1]);
Graph(Petri.Model(decoration(death_petri)))

# #### Step 2: Define a presentation of the free biproduct category
# that encodes the domain specific information

@present Predation(FreeBiproductCategory) begin
    prey::Ob
    predator::Ob
    birth::Hom(prey,prey)
    predation::Hom(prey⊗predator,predator)
    death::Hom(predator,predator)
end;

rabbits,wolves,birth,predation,death = generators(Predation);

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
                 rabbits=>PetriCospanOb(1),wolves=>PetriCospanOb(1),
                 birth=>birth_petri, predation=>predation_petri, death=>death_petri));

# #### Step 3: Generate models using the hom expression or program notations

lotka_volterra = (birth ⊗ id(wolves)) ⋅ predation ⋅ death
lotka_petri = Petri.Model(decoration(F(lotka_volterra)))
display_wd(lotka_volterra)
#-
Graph(lotka_petri)

# Generate appropriate vector fields, define parameters, and visualize solution

u0 = [100, 10];
p = [.3, .015, .7];
prob = ODEProblem(lotka_petri,u0,(0.0,100.0),p);
sol = solve(prob,Tsit5(),abstol=1e-8);
plot(sol)

# There is also a second syntax that is easier to write for programmers
# than the hom expression syntax. Here is an example of the same model
# as before along with a test of equivalency

lotka_volterra2 = @program Predation (r::prey, w::predator) begin
  r_2 = birth(r)
  w_2 = predation(r_2, w)
  return death(w_2)
end
lotka_petri2 = Petri.Model(decoration(F(to_hom_expr(FreeBiproductCategory, lotka_volterra2))))
lotka_petri == lotka_petri2

# #### Step 4: Extend your presentation to handle more complex phenomena
# such as a small food chain

@present DualPredation <: Predation begin
    Predator::Ob
    Predation::Hom(predator⊗Predator,Predator)
    Death::Hom(Predator,Predator)
end;

fish,Fish,Shark,birth,predation,death,Predation,Death = generators(DualPredation);

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
                 fish=>PetriCospanOb(1),Fish=>PetriCospanOb(1),
                 birth=>birth_petri, predation=>predation_petri, death=>death_petri,
                 Shark=>PetriCospanOb(1),Predation=>predation_petri, Death=>death_petri));

# Define a new model where fish are eaten by Fish which are then eaten by Sharks

dual_lv = @program DualPredation (fish::prey, Fish::predator, Shark::Predator) begin
  f_2 = birth(fish)
  F_2 = predation(f_2, Fish)
  F_3 = death(F_2)
  S_2 = Predation(F_3, Shark)
  S_3 = Death(S_2)
end
display_wd(dual_lv)
#-
dual_lv_petri = Petri.Model(decoration(F(to_hom_expr(FreeBiproductCategory, dual_lv))))
Graph(dual_lv_petri)

# Generate a new solver, provide parameters, and analyze results

u0 = [100, 10, 2];
p = [.3, .015, .7, .017, .35];
prob = ODEProblem(dual_lv_petri,u0,(0.0,100.0),p);
sol = solve(prob,Tsit5(),abstol=1e-6);
plot(sol)
