# # [Multi-City COVID-19 Model](@id covid_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/covid.ipynb)

using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using Petri
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Theories
using Catlab.Programs
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.Graphics.Graphviz: run_graphviz

import Catlab.Theories: id

# helper functions
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);
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
dictreplace(f::Dict{K,N}, d::Dict{K,V}) where {K,N,V} = Dict{N,V}(Base.map(x->Pair(f[x[1]], x[2]), collect(d))) # hide
dictreplace(f::Dict, ts::Vector{Tuple{S,T}}) where {S<:Dict, T<:Dict} = [(dictreplace(f, t[1]), dictreplace(f, t[2])) for t in ts] # hide
statereplace(p::PetriWithRates, f::Dict{K,N}) where {K,N} = begin
    S = map(s->f[s], p.m.S)
    u0 = LVector(;zip(S, zeros(length(S)))...)
    map(k->u0[f[k]] = p.u0[k], keys(p.u0))
    PetriWithRates(S, dictreplace(f, p.m.Δ), p.rates, u0) # hide
end

pob = PetriCospanOb(1)
spontaneous_petri(rate, u0) = PetriCospan([1], PetriWithRates(1:2, [(Dict(1=>1), Dict(2=>1))], [rate], u0), [2])
transmission_petri(rate, u0) = PetriCospan([1,2], PetriWithRates(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))], [rate], u0), [2])
exposure_petri(rate, u0) = PetriCospan([1, 2], PetriWithRates(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))], [rate], u0), [3])

# Extend the Infectious Disease presentation,
# get the new generators, and update the functor

@present EpiCoexist <: InfectiousDiseases begin
    I2::Ob
    A::Ob
    R2::Ob
    exposure_e::Hom(S⊗E,E)
    exposure_i2::Hom(S⊗I2,E)
    exposure_a::Hom(S⊗A,E)
    progression::Hom(I,I2)
    asymptomatic_infection::Hom(E,A)
    recover_late::Hom(R,R2)
    asymptomatic_recovery::Hom(A,R)
    recovery2::Hom(I2,R)
    death2::Hom(I2,D)
end;

S,E,I,R,D,I2,A,R2,transmission,exposure,illness,recovery,death,exposure_e,exposure_i2,exposure_a,progression,asymptomatic_infection,recover_late,asymptomatic_recovery,recovery2, death2 = generators(EpiCoexist);

pop = [33540260 - 4*1000, 30895128 - 4*1000, 0]
N = sum(pop) + 2*4*1000
social_mixing_rate = [2.3276/N, 1.2015/N, .93828/N]
fatality_rate = [0.0281, 0.1868, 0]

F(ex, n) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(
    S=>pob,
    E=>pob,
    A=>pob,
    I=>pob,
    I2=>pob,
    R=>pob,
    R2=>pob,
    D=>pob,
    transmission=>transmission_petri(0, [0,0]),
    exposure=>exposure_petri(.1*social_mixing_rate[n], [pop[n],1000, 0]),
    exposure_e=>exposure_petri(.001*social_mixing_rate[n], [0,1000,0]),
    exposure_i2=>exposure_petri(.6*social_mixing_rate[n], [0,1000,0]),
    exposure_a=>exposure_petri(.5*social_mixing_rate[n], [0,1000,0]),
    progression=>spontaneous_petri(.25, [0,0]),
    asymptomatic_infection=>spontaneous_petri(.86/.14*.2, [0,0]),
    illness=>spontaneous_petri(.2, [0,0]),
    asymptomatic_recovery=>spontaneous_petri(1/15,[0,0]),
    recovery=>spontaneous_petri(0,[0,0]),
    recovery2=>spontaneous_petri(1/6,[0,0]),
    recover_late=>spontaneous_petri(1/15,[0,0]),
    death=>spontaneous_petri(0,[0,0]),
    death2=>spontaneous_petri((1/15)*(fatality_rate[n]/(1-fatality_rate[n])),[0,0])))

coexist = @program EpiCoexist (s::S, e::E, i::I, i2::I2, a::A, r::R, r2::R2, d::D) begin
    e_2 = exposure(s, i)
    e_3 = exposure_i2(s, i2)
    e_4 = exposure_a(s, a)
    e_5 = exposure_e(s, e)
    e_all = [e, e_2, e_3, e_4, e_5]
    a_2 = asymptomatic_infection(e_all)
    a_all = [a, a_2]
    r_2 = asymptomatic_recovery(a_all)
    i_2 = illness(e_all)
    i_all = [i, i_2]
    i2_2 = progression(i)
    i2_all = [i2, i2_2]
    d_2 = death2(i2_all)
    r_3 = recovery2(i2_all)
    r_all = [r, r_2, r_3]
    r2_2 = recover_late(r_all)
    r2_all = [r2, r2_2]
    d_all = [d, d_2]
    return s, e_all, i_all, i2_all, a_all, r_all, r2_all, d_all
end
# display_wd(coexist)
# save_wd(coexist, "coexist_wd.svg")
coexist = to_hom_expr(FreeBiproductCategory, coexist)

coexist_petri = decoration(F(coexist, 1))
coexist_petri = statereplace(coexist_petri, Dict(1=>:S,2=>:E,3=>:R2,4=>:D,5=>:I2,6=>:A,7=>:R,8=>:I))
tspan = (0.0,150.0)
prob = ODEProblem(coexist_petri.m,coexist_petri.u0,tspan,coexist_petri.rates)
sol = solve(prob,Tsit5())

plot(sol)

crossexposure = @program EpiCoexist (s::S, e::E, i::I, i2::I2, a::A, r::R, r2::R2, d::D,
                               s′::S, e′::E, i′::I, i2′::I2, a′::A, r′::R, r2′::R2, d′::D) begin
    e_2 = exposure(s, i′)
    e_3 = exposure_i2(s, i2′)
    e_4 = exposure_a(s, a′)
    e_5 = exposure_e(s, e′)
    e_all = [e, e_2, e_3, e_4, e_5]
    return s, e_all, i, i2, a, r, r2, d,
           s′, e′, i′, i2′, a′, r′, r2′, d′
end
crossexposure = to_hom_expr(FreeBiproductCategory, crossexposure)

#display_wd(crossexposure)
# save_wd(crossexposure, "crossexposure_wd.svg")

# 2 generation cross exposure + coexist model
population = otimes(S, E, I, I2, A, R, R2, D)
pop_hom = F(population, 1)
co_1 = F(coexist, 1)
co_2 = F(coexist, 2)
cross = F(crossexposure, 3)
twogen = cross ⋅ σ(pop_hom, pop_hom) ⋅ cross ⋅ (co_1 ⊗ co_2)

twogen_petri = decoration(twogen)
twogen_petri = statereplace(twogen_petri, Dict(1=>:S,2=>:E,3=>:I,4=>:I2,5=>:A,6=>:R,7=>:R2,8=>:D,
                                               9=>:S′,10=>:E′,11=>:I′,12=>:I2′,13=>:A′,14=>:R′,15=>:R2′,16=>:D′))
tspan = (0.0,150.0)
prob = ODEProblem(twogen_petri.m,twogen_petri.u0,tspan,twogen_petri.rates)
sol = solve(prob,Tsit5())

plot(sol)

sol_total = transpose(hcat([[s[:S] + s[:S′], s[:E] + s[:E′], s[:A] + s[:A′], s[:I] + s[:I′], s[:I2] + s[:I2′], s[:R] + s[:R′], s[:R2] + s[:R2′], s[:D] + s[:D′]] for s in sol.u]...))
plot(sol.t, sol_total, label=["S" "E" "A" "I" "I2" "R" "R2" "D"])
# display_wd(twogen)
# Graph(decoration(F(twogen)))

# n generation cross exposure + coexist model
n = 9
population = otimes(S, E, I, I2, A, R, R2, D)
pops = [population for i in 1:n]

single_exposure = foldl(⋅, [otimes(map(id, pops[1:i])...,crossexposure⋅σ(pops[i+1:i+2]...),map(id, pops[i+3:end])...) for i in 0:(n-2)])
# display_wd(single_exposure)
# save_wd(single_exposure, "$(n)gen_single_crossexposure.svg")

ngen_exposure = foldl(⋅, [single_exposure for i in 1:n])
# display_wd(ngen_exposure)
# save_wd(ngen_exposure, "$(n)gen_crossexposure.svg")

ngen_coexist = ngen_exposure ⋅ foldl(⊗, [coexist for i in 1:n])
# display_wd(ngen_coexist)
#save_wd(ngen_coexist, "$(n)gen_coexist_wd.svg")

ngen_coexist_petri = decoration(F(ngen_coexist, 1))
#Graph(ngen_coexist_petri)
#save_graph(Graph(ngen_coexist_petri), "$(n)gen_coexist_petri.svg")

# Generate some simpler diagrams to expose the hierarchical structure
@present CoexistOverview(FreeBiproductCategory) begin
    Pop::Ob
    coexist::Hom(Pop,Pop)
    crossexposure::Hom(Pop⊗Pop,Pop⊗Pop)
end
pop′,coexist′,crossexposure′ = generators(CoexistOverview)
pops = [pop′ for i in 1:n]
single_exposure = foldl(⋅, [otimes(map(id, pops[1:i])...,crossexposure′⋅σ(pops[i+1:i+2]...),map(id, pops[i+3:end])...) for i in 0:(n-2)])
#display_wd(single_exposure)
#save_wd(single_exposure, "$(n)gen_single_crossexposure_overview.svg")
ngen_exposure = foldl(⋅, [single_exposure for i in 1:n])
#display_wd(ngen_exposure)
#save_wd(ngen_exposure, "$(n)gen_crossexposure_overview.svg")
ngen_coexist = ngen_exposure ⋅ foldl(⊗, [coexist′ for i in 1:n])
#display_wd(ngen_coexist)
#save_wd(ngen_coexist, "$(n)gen_coexist_overview.svg")

# Even more generalized diagram
@present AllCoexist(FreeBiproductCategory) begin
    TotalPop::Ob
    AllCoexist::Hom(TotalPop,TotalPop)
    AllCrossExposure::Hom(TotalPop,TotalPop)
end
allpop,allcoexist,allcrossexposure = generators(AllCoexist)
all_coexist = allcrossexposure ⋅ allcoexist
#display_wd(all_coexist)
#save_wd(all_coexist, "all_coexist_wd.svg")


# Attempt to compute solution to single generation COEXIST model


# Define time frame and initial parameters
coexist_pc = F(coexist)
coexist_petri = decoration(coexist_pc)
#Graph(coexist_petri)

tspan = (0.0,150.0)
u0 = zeros(Float64, base(coexist_pc).n)
u0[2]  = 13000000

fatality_rate = 0.146
β = [0,0,0,0,1/4,.86/.14*.2,1/(10-4),1/15,1/5,1/15,(1/15)*(fatality_rate/(1-fatality_rate))]

# Generate, solve, and visualize resulting ODE
prob = ODEProblem(coexist_petri,u0,tspan,β);
sol = solve(prob,Tsit5());

plot(sol)

map(x->x[4], sol.u)

coexist_petri2 = statereplace(coexist_petri, Dict(1=>:S,2=>:E,3=>:R2,4=>:D,5=>:I2,6=>:A,7=>:R,8=>:I))
Graph(coexist_petri2)
u0 = LVector([s=0 for s in coexist_petri2.S]...)
u0[:S]  = 57345080
u0[:E]  = 1000
u0[:I]  = 1000
u0[:I2] = 1000
u0[:A]  = 1000
N = sum(u0)

β = [.001*social_mixing_rate,0.1*social_mixing_rate,0.6*social_mixing_rate,0.5*social_mixing_rate,1/4,.86/.14*.2,1/(10-4),1/15,1/5,1/15,(1/15)*(fatality_rate/(1-fatality_rate))] # hide

tspan = (0.0,150.0)
prob = ODEProblem(coexist_petri2,u0,tspan,β); # hide
sol = solve(prob,Tsit5()); # hide

coexist_petri2.Δ
u0