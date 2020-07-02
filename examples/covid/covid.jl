using Petri
using OrdinaryDiffEq
using Plots
using AlgebraicPetri
using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.Graphics.Graphviz: Graph

ob = PetriCospanOb(1)

spontaneous_petri = PetriCospan(
        Cospan(FinOrdFunction([1], 2),
               FinOrdFunction([2], 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1), Dict(2=>1))]))

transmission_petri = PetriCospan(
        Cospan(FinOrdFunction([1], 2),
               FinOrdFunction([2], 2)
        ), id(PetriFunctor), Petri.Model([1, 2], [(Dict(1=>1, 2=>1), Dict(2=>2))]))

exposure_petri = PetriCospan(
        Cospan(FinOrdFunction([1, 2], 3),
               FinOrdFunction([3, 2], 3)
        ), id(PetriFunctor), Petri.Model([1, 2, 3], [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]))
travel_petri = PetriCospan(
        Cospan(FinOrdFunction([1,2,3,4,5], 8),
               FinOrdFunction([6,7,8], 8)
        ), id(PetriFunctor), Petri.Model(collect(1:8), [(Dict(1=>1, 2=>1, 3=>1), Dict(6=>1, 7=>1, 8=>1))]))

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
    travel::Hom(S⊗E⊗I⊗R⊗D,S⊗E⊗I)
end

S,E,I,R,D,transmission,exposure,illness,recovery,death,travel = generators(Epidemiology)

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true)

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
        S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri,travel=>travel_petri))

# define model
sir = transmission ⋅ recovery
# get resulting petri net
p_sir = decoration(F(sir))

# display wiring diagram and petri net visualization
display_wd(sir)
Graph(p_sir)

# define initial states and transition rates
u0 = [10.0, 1, 0]
p = [0.4, 0.4]
# create and solve ODE problem
prob = ODEProblem(toODE(p_sir),u0,(0.0,7.5),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)

sei = exposure ⋅ (illness ⊗ id(I)) ⋅ ∇(I)

seir = sei ⋅ recovery
p_seir = decoration(F(seir))

display_wd(seir)
Graph(p_seir)

# define initial states and transition rates
u0 = [10.0, 1, 0, 0]
p = [0.9, 0.2, 0.5]
# create and solve ODE problem
prob = ODEProblem(toODE(p_seir),u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)

seird = sei ⋅ Δ(I) ⋅ (death ⊗ recovery)
p_seird = decoration(F(seird))

display_wd(seird)
Graph(p_seird)

# define initial states and transition rates
u0 = [10.0, 1, 0, 0, 0]
p = [0.9, 0.2, 0.5, 0.1]
# create and solve ODE problem
prob = ODEProblem(toODE(p_seird),u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
plot(sol)

# TODO: Add support for types so we can simplify to this
# seir = exposure ⋅ (illness ⊗ recovery)
# seird = seir ⋅ (death ⊗ id(R))
# display_wd(seird)

se  = id(S) ⊗ id(E)
sei = id(S) ⊗ id(E) ⊗ id(I)
seird_city = (((Δ(S) ⊗ id(E)) ⋅ (id(S) ⊗ σ(S,E))) ⊗ id(I)) ⋅ (se ⊗ exposure) ⋅ (id(S) ⊗ (∇(E) ⋅ Δ(E)) ⊗ id(I)) ⋅ (se ⊗ ((illness ⊗ id(I)) ⋅ (∇(I) ⋅ Δ(I)) ⋅ (id(I) ⊗ (Δ(I) ⋅(recovery ⊗ death))))) ⋅ travel

p_seird_city = decoration(F(seird_city))

display_wd(seird_city)
Graph(p_seird_city)

ncities(city,n::Int) = compose([city for i in 1:n]...)
seird_2 = ncities(seird_city, 2)
display_wd(seird_2)

seird_3 = ncities(seird_city, 3)
p_seird_3 = decoration(F(seird_3))
display_wd(seird_3)
Graph(p_seird_3)