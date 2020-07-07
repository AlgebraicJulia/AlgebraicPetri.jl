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
using Catlab.Graphics.Graphviz: Graph, run_graphviz

import Catlab.Theories: id

# A few helper functions
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true)
save_wd(ex, fname::AbstractString) = begin
    g = display_wd(ex)
    open(fname, "w") do io
        run_graphviz(io, g, format="svg")
    end
end
save_graph(g, fname::AbstractString) = begin
    open(fname, "w") do io
        run_graphviz(io, g, format="svg")
    end
end

splot(soln, fname) = begin
    p = plot(soln, linewidth=5, legend=false)
    savefig(p, fname)
end

splot(soln, t, fname) = begin
    p = plot(t, soln, linewidth=5, legend=false)
    savefig(p, fname)
end

splotchannels(sol, dir) = begin
    mkpath(dir)
    jp(p) = joinpath(dir, p)
    splot(sol, jp("solution.svg"))
    pltdims = [3,8,13]
    splot(sol[pltdims,:]', sol.t, jp("infecteds.svg"))
    splot(sum(sol(sol.t, idxs=pltdims), dims=1)', sol.t, jp("infecteds_sum.svg"))
    pltdims = [3,8,13] .+ 1
    splot(sol[pltdims,:]', sol.t, jp( "recovered.svg" ))
    splot(sum(sol(sol.t, idxs=pltdims), dims=1)', sol.t, jp( "recovered_sum.svg" ))
end

id(args...) = foldl((x,y)->id(x) ⊗ id(y), args)

# Step 1: Define building block Petri Net models
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
        Cospan(FinOrdFunction([1,2,3], 6),
               FinOrdFunction([4,5,6], 6)
        ), id(PetriFunctor), Petri.Model(collect(1:6), [(Dict(1=>1), Dict(4=>1)),
                                                        (Dict(2=>1), Dict(5=>1)),
                                                        (Dict(3=>1), Dict(6=>1))]))

[save_graph(Graph(decoration(p)), pname) for (p,pname) in [(spontaneous_petri, "spont.svg"),
                                               (transmission_petri, "transmission.svg"),
                                               (exposure_petri, "exposure.svg"),
                                               (travel_petri, "travel.svg")]]

# Step 2: Define a strongly type presentation of the
#         Free Biproduct Category for the desired domain
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
    travel::Hom(S⊗E⊗I,S⊗E⊗I)
end

# Create the generators
S,E,I,R,D,transmission,exposure,illness,recovery,death,travel = generators(Epidemiology)

# Define a functor from the generators to the building block Petri Nets
F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
        S=>ob, E=>ob, I=>ob, R=>ob, D=>ob,
        transmission=>transmission_petri, exposure=>exposure_petri,
        illness=>spontaneous_petri, recovery=>spontaneous_petri, death=>spontaneous_petri,travel=>travel_petri))

# Step 3: Create, visualize, and solve possible models

# SIR, SEIR, SEIRD Basic Epidemiology Models:

# define model
sir = transmission ⋅ recovery
# get resulting petri net
p_sir = decoration(F(sir))

# display wiring diagram and petri net visualization
display_wd(sir)
save_wd(sir, "sir.svg")
g = Graph(p_sir)
save_graph(g, "sir_graph.svg")

# define initial states and transition rates
u0 = [10.0, 1, 0]
p = [0.4, 0.4]
# create and solve ODE problem
prob = ODEProblem(toODE(p_sir),u0,(0.0,7.5),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
splot(sol, "sir_soln.svg")

# define model
sei = exposure ⋅ (illness ⊗ id(I)) ⋅ ∇(I)

seir = sei ⋅ recovery
# get resulting petri net
p_seir = decoration(F(seir))

# display wiring diagram and petri net visualization
display_wd(seir)
save_wd(seir, "seir.svg")
g = Graph(p_seir)
save_graph(g, "seir_graph.svg")

# define initial states and transition rates
u0 = [10.0, 1, 0, 0]
p = [0.9, 0.2, 0.5]
# create and solve ODE problem
prob = ODEProblem(toODE(p_seir),u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
splot(sol, "seir_soln.svg")

# define model
seird = sei ⋅ Δ(I) ⋅ (death ⊗ recovery)
# get resulting petri net
p_seird = decoration(F(seird))

# display wiring diagram and petri net visualization
display_wd(seird)
save_wd(seird, "seird.svg")
g = Graph(p_seird)
save_graph(g, "seird_graph.svg")

# define initial states and transition rates
u0 = [10.0, 1, 0, 0, 0]
p = [0.9, 0.2, 0.5, 0.1]
# create and solve ODE problem
prob = ODEProblem(toODE(p_seird),u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
splot(sol, "seird_soln.svg")

# TODO: Add support for types so we can simplify to this
# seir = exposure ⋅ (illness ⊗ recovery)
# seird = seir ⋅ (death ⊗ id(R))
# display_wd(seird)


# COVID-19 TRAVEL MODEL:
# SEIRD City Model with travel as S ⊗ E ⊗ I → S ⊗ E ⊗ I
# Manually defined Hom:
#       seird_city = (((Δ(S) ⊗ id(E)) ⋅ (id(S) ⊗ σ(S,E))) ⊗ id(I)) ⋅ (id(S, E) ⊗ exposure) ⋅ (id(S) ⊗ (∇(E) ⋅ Δ(E)) ⊗ id(I)) ⋅ (id(S, E) ⊗ ((illness ⊗ id(I)) ⋅ (∇(I) ⋅ Δ(I)) ⋅ (id(I) ⊗ (Δ(I) ⋅(recovery ⊗ death))))) ⋅ (travel ⊗ ◊(R) ⊗ ◊(D))
# use the program interface for easier model definition
seird_city = @program Epidemiology (s::S, e::E, i::I) begin
    e2, i2 = exposure(s, i)
    i3 = illness(e2)
    d = death(i3)
    r = recovery(i3)
    return travel(s, [e, e2], [i2, i3])
end
seird_city = to_hom_expr(FreeBiproductCategory, seird_city)

display_wd(seird_city)
save_wd(seird_city, "seird_city.svg")

g = Graph(decoration(F(seird_city)))
save_graph(g, "seird_graph.svg")

# function to compose n city models together
ncities(city,n::Int) = compose([city for i in 1:n]...)

# create a 3 city SEIRD models
seird_3 = ncities(seird_city, 3)
pc_seird_3 = F(seird_3)
p_seird_3 = decoration(pc_seird_3)
display_wd(seird_3)
save_wd(seird_3, "seird_3.svg")
g = Graph(p_seird_3)
save_graph(g, "seird_3_graph.svg")

pltdims = [4,9, 14]
# Define time frame, 3 months
tspan = (0.0,90.0)
# Define initial states
u0 = zeros(Float64, base(pc_seird_3).n)
u0[1]  = 10000
u0[6]  = 10000
u0[11] = 10000
u0[2]  = 1
# Define transition rates
seirdparams(n::Int, k::Number) = begin
    βseird = [10/sum(u0), 1/2, 1/5, 1/16]
    βtravel = [1/20, 1/200, 1/20]/100k
    β = vcat(βseird, βtravel)
    return foldl(vcat, [(1-(0/(2n)))*β for i in 1:n])
end
params = seirdparams(3, 5)
# Generate and solve resulting ODE
prob = ODEProblem(toODE(p_seird_3),u0,tspan,params)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution


    
splotchannels(sol, "seird")

include("dynamic_rates.jl")

asymptotic(a,b,k=1) = t->(b + (a-b)/(k*t+1))
triangleasm(a,b,k=1) = t->max(b, (1-t/k)*(a-b)+b)
trianglewave(a,b,k=1) = t->((a-b)/π)*asin(sin(t*2π/k))+2b
coswave(a,b,k=1) = t->(a-b)*(cos(t*2π/k)/2)+( a+b )/2
sincwave(a,b,k=1) = t->(a-b)*(sinc(t*2π/k)/2)+( a+b )/2
modsincwave(a,b,k) = t-> max(b, (a-b)*(sinc(t*2π/k))+(b))
# sqwave(a,b,k=1) = t->(a-b)*(cos(t*2π/k)/2)+( a+b )/2
dynseirdparams(f, a,b,period, n::Int, scale::Number) = begin
    βseird = [f(a,b,period), 1/2, 1/5, 1/16]
    βtravel = [1/20, 1/200, 1/20]/100scale
    β = vcat(βseird, βtravel)
    return foldl(vcat, [β for i in 1:n])
end

waveparams(f,a,b,p) = begin
    dynseirdparams(f,a,b,p,3,5)
end

function drawsol(f, a, b, p, dir)
    tspan = (0, 240.0)
    println("INFO: Processing $dir")
    dynparams = waveparams(f, a/sum(u0), b/sum(u0),p)
    mkpath(dir)
    prob = ODEProblem(DynamicRates.toODE(p_seird_3, nothing),u0,tspan, dynparams)
    sol = OrdinaryDiffEq.solve(prob,Tsit5(), saveat=1:1:tspan[2])
    splotchannels(sol, dir)
    savefig(plot(sol.t, f(a,b,p), ylim=(0,12), legend=false, linewidth=5), "$dir/transmissability.svg")
end

a,b = 10, 1
dir = "dynamic/mitigation/"
drawsol(asymptotic,  a,b,1/6, "$dir/asymptotic")
drawsol(triangleasm, a,b,25 , "$dir/trianglesm")
drawsol(coswave,     a,b,25 , "$dir/coswave")
drawsol(sincwave,    a,b,100, "$dir/sincwave")
drawsol(modsincwave, a,b,100, "$dir/modsincwave")

a,b = 10, 1/2
dir = "dynamic/containment"
drawsol(asymptotic,  a,b,1/6, "$dir/asymptotic")
drawsol(triangleasm, a,b, 25, "$dir/trianglesm")
drawsol(coswave,     a,b, 25, "$dir/coswave")
drawsol(sincwave,    a,b,100, "$dir/sincwave")
drawsol(modsincwave, a,b,100, "$dir/modsincwave")
