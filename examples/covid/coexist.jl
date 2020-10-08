# # [COEXIST Multi-Generational COVID Model](@id coexist_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/coexist.ipynb)

using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using Petri
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Programs
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

const EpiRxnNet = LabelledReactionNet{Number,Int}
const OpenEpiRxnNet = OpenLabelledReactionNet{Number,Int}
const OpenEpiRxnNetOb = OpenLabelledReactionNetOb{Number,Int}

# UPDATE THESE DEFINITIONS
ob(x::Symbol,n::Int) = codom(Open([x], EpiRxnNet(x=>n), [x]))
spontaneous_petri(x::Symbol, xn::Int, y::Symbol, yn::Int, z::Symbol, zr::Number) = Open([x], EpiRxnNet((x=>xn,y=>yn), (z,zr)=>(x=>y)), [y])
transmission_petri(name::Symbol, S::Symbol, s::Int, I::Symbol, i::Int, inf::Number) = Open([S], EpiRxnNet((S=>s,I=>i), (name,inf)=>((S,I)=>(I,I))), [I])
exposure_petri(name::Symbol, S::Symbol, s::Int, X::Symbol,x::Int,E::Symbol,e::Int,exp::Number) = Open([S, X], EpiRxnNet((S=>s,X=>x,E=>e), (name,exp)=>((S,X)=>(E,X))), [E])

# Extend the Infectious Disease presentation,
# get the new generators, and update the functor

pop = [33540260 - 4*1000, 30895128 - 4*1000, 0]
N = sum(pop) + 2*4*1000
social_mixing_rate = [2.3276/N, 1.2015/N, .93828/N]
fatality_rate = [0.0281, 0.1868, 0]

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

F(ex, n) = functor((OpenEpiRxnNetOb, OpenEpiRxnNet), ex, generators=Dict(
    S=>ob(Symbol(:S, n), pop[n]),
    E=>ob(Symbol(:E, n), 1000),
    A=>ob(Symbol(:A, n), 1000),
    I=>ob(Symbol(:I, n), 1000),
    I2=>ob(Symbol(:I2, n), 1000),
    R=>ob(Symbol(:R, n), 0),
    R2=>ob(Symbol(:R2, n), 0),
    D=>ob(Symbol(:D, n), 0),
    transmission=>transmission_petri(:transmission, :S, 0, :I, 0, 0),
    exposure=>exposure_petri(Symbol(:exp_, n), Symbol(:S,n), pop[n], Symbol(:I,n), 1000, Symbol(:E,n), 1000, .1*social_mixing_rate[n]),
    exposure_e=>exposure_petri(Symbol(:exp_e, n), Symbol(:S,n), pop[n], Symbol(:E,n),1000, Symbol(:E,n),1000, .001*social_mixing_rate[n]),
    exposure_i2=>exposure_petri(Symbol(:exp_i2, n), Symbol(:S,n), pop[n], Symbol(:I2,n), 1000, Symbol(:E,n),1000, .6*social_mixing_rate[n]),
    exposure_a=>exposure_petri(Symbol(:exp_a, n), Symbol(:S,n), pop[n], Symbol(:A,n),1000, Symbol(:E,n),1000, .5*social_mixing_rate[n]),
    progression=>spontaneous_petri(Symbol(:I,n), 1000, Symbol(:I2,n), 1000, Symbol(:prog_, n), .25),
    asymptomatic_infection=>spontaneous_petri(Symbol(:E,n), 1000, Symbol(:A,n), 1000, Symbol(:asymp_, n), .86/.14*.2),
    illness=>spontaneous_petri(Symbol(:E,n), 1000, Symbol(:I,n), 1000, Symbol(:ill_, n), .2),
    asymptomatic_recovery=>spontaneous_petri(Symbol(:A,n), 1000, Symbol(:R,n), 0, Symbol(:arec_, n), 1/15),
    recovery=>spontaneous_petri(Symbol(:I,n), 0, Symbol(:R,n), 0, Symbol(:rec_, n), 0),
    recovery2=>spontaneous_petri(Symbol(:I2,n), 1000, Symbol(:R,n), 0, Symbol(:rec_, n), 1/6),
    recover_late=>spontaneous_petri(Symbol(:R,n), 0, Symbol(:R2,n), 0, Symbol(:rec2_, n), 1/15),
    death=>spontaneous_petri(Symbol(:I,n), 0, Symbol(:D,n), 0, Symbol(:death_, n), 0),
    death2=>spontaneous_petri(Symbol(:I2,n), 1000, Symbol(:D,n), 0, Symbol(:death2_, n), (1/15)*(fatality_rate[n]/(1-fatality_rate[n])))))

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
coexist = to_hom_expr(FreeBiproductCategory, coexist)

@present EpiCrossExposure(FreeBiproductCategory) begin
    S::Ob
    E::Ob
    A::Ob
    I::Ob
    I2::Ob
    R::Ob
    R2::Ob
    D::Ob
    S′::Ob
    E′::Ob
    A′::Ob
    I′::Ob
    I2′::Ob
    R′::Ob
    R2′::Ob
    D′::Ob
    exposure::Hom(S⊗I′,E)
    exposure_e::Hom(S⊗E′,E)
    exposure_a::Hom(S⊗A′,E)
    exposure_i2::Hom(S⊗I2′,E)
    exposure′::Hom(S′⊗I,E′)
    exposure_e′::Hom(S′⊗E,E′)
    exposure_a′::Hom(S′⊗A,E′)
    exposure_i2′::Hom(S′⊗I2,E′)
end;

ce_S,ce_E,ce_A,ce_I,ce_I2,ce_R,ce_R2,ce_D,
ce_S′,ce_E′,ce_A′,ce_I′,ce_I2′,ce_R′,ce_R2′,ce_D′,
ce_exposure, ce_exposure_e, ce_exposure_a, ce_exposure_i2,
ce_exposure′, ce_exposure_e′, ce_exposure_a′, ce_exposure_i2′ = generators(EpiCrossExposure);

F_cx(ex, x,y) = functor((OpenEpiRxnNetOb, OpenEpiRxnNet), ex, generators=Dict(
    ce_S=>ob(Symbol(:S, x), pop[x]),
    ce_E=>ob(Symbol(:E, x), 1000),
    ce_A=>ob(Symbol(:A, x), 1000),
    ce_I=>ob(Symbol(:I, x), 1000),
    ce_I2=>ob(Symbol(:I2, x), 1000),
    ce_R=>ob(Symbol(:R, x), 0),
    ce_R2=>ob(Symbol(:R2, x), 0),
    ce_D=>ob(Symbol(:D, x), 0),
    ce_S′=>ob(Symbol(:S, y), pop[y]),
    ce_E′=>ob(Symbol(:E, y), 1000),
    ce_A′=>ob(Symbol(:A, y), 1000),
    ce_I′=>ob(Symbol(:I, y), 1000),
    ce_I2′=>ob(Symbol(:I2, y), 1000),
    ce_R′=>ob(Symbol(:R, y), 0),
    ce_R2′=>ob(Symbol(:R2, y), 0),
    ce_D′=>ob(Symbol(:D, y), 0),
    ce_exposure=>exposure_petri(Symbol(:exp_, x,y), Symbol(:S,x), pop[x], Symbol(:I,y), 1000, Symbol(:E,x), 1000, .1*social_mixing_rate[x]),
    ce_exposure_e=>exposure_petri(Symbol(:exp_e, x,y), Symbol(:S,x), pop[x], Symbol(:E,y),1000, Symbol(:E,x),1000, .001*social_mixing_rate[x]),
    ce_exposure_a=>exposure_petri(Symbol(:exp_a, x,y), Symbol(:S,x), pop[x], Symbol(:A,y),1000, Symbol(:E,x),1000, .5*social_mixing_rate[x]),
    ce_exposure_i2=>exposure_petri(Symbol(:exp_i2, x,y), Symbol(:S,x), pop[x], Symbol(:I2,y), 1000, Symbol(:E,x),1000, .6*social_mixing_rate[x]),
    ce_exposure′=>exposure_petri(Symbol(:exp_, y,x), Symbol(:S,y), pop[y], Symbol(:I,x), 1000, Symbol(:E,y), 1000, .1*social_mixing_rate[y]),
    ce_exposure_e′=>exposure_petri(Symbol(:exp_e, y,x), Symbol(:S,y), pop[y], Symbol(:E,x),1000, Symbol(:E,y),1000, .001*social_mixing_rate[y]),
    ce_exposure_a′=>exposure_petri(Symbol(:exp_a, y,x), Symbol(:S,y), pop[y], Symbol(:A,x),1000, Symbol(:E,y),1000, .5*social_mixing_rate[y]),
    ce_exposure_i2′=>exposure_petri(Symbol(:exp_i2, y,x), Symbol(:S,y), pop[y], Symbol(:I2,x), 1000, Symbol(:E,y),1000, .6*social_mixing_rate[y])
    ))

crossexposure = @program EpiCrossExposure (s::S, e::E, i::I, i2::I2, a::A, r::R, r2::R2, d::D,
                                           s′::S′, e′::E′, i′::I′, i2′::I2′, a′::A′, r′::R′, r2′::R2′, d′::D′) begin
    e_2 = exposure(s, i′)
    e_3 = exposure_i2(s, i2′)
    e_4 = exposure_a(s, a′)
    e_5 = exposure_e(s, e′)
    e_all = [e, e_2, e_3, e_4, e_5]
    e′_2 = exposure′(s′, i)
    e′_3 = exposure_i2′(s′, i2)
    e′_4 = exposure_a′(s′, a)
    e′_5 = exposure_e′(s′, e_all)
    e′_all = [e′, e′_2, e′_3, e′_4, e′_5]
    return s, e_all, i, i2, a, r, r2, d,
           s′, e′_all, i′, i2′, a′, r′, r2′, d′
end
crossexposure = to_hom_expr(FreeBiproductCategory, crossexposure)

@present ThreeCoexist(FreeBiproductCategory) begin
    Pop1::Ob
    Pop2::Ob
    Pop3::Ob
    crossexp12::Hom(Pop1⊗Pop2,Pop1⊗Pop2)
    crossexp13::Hom(Pop1⊗Pop3,Pop1⊗Pop3)
    crossexp23::Hom(Pop2⊗Pop3,Pop2⊗Pop3)
    coex1::Hom(Pop1,Pop1)
    coex2::Hom(Pop2,Pop2)
    coex3::Hom(Pop3,Pop3)
end;


Pop1, Pop2, Pop3, crossexp12, crossexp13, crossexp23, coex1, coex2, coex3 = generators(ThreeCoexist);

F_tcx(ex) = functor((OpenEpiRxnNetOb, OpenEpiRxnNet), ex, generators=Dict(
    Pop1=>F(otimes(S,E,I,I2,A,R,R2,D),1),
    Pop2=>F(otimes(S,E,I,I2,A,R,R2,D),2),
    Pop3=>F(otimes(S,E,I,I2,A,R,R2,D),3),
    crossexp12=>F_cx(crossexposure,1,2),
    crossexp13=>F_cx(crossexposure,1,3),
    crossexp23=>F_cx(crossexposure,2,3),
    coex1=>F(coexist,1),
    coex2=>F(coexist,2),
    coex3=>F(coexist,3)
    ))

threeNCoexist = @program ThreeCoexist (pop1::Pop1, pop2::Pop2, pop3::Pop3) begin
    pop1′, pop2′ = crossexp12(pop1, pop2)
    pop1′′, pop3′ = crossexp13(pop1′, pop3)
    pop2′′, pop3′′ = crossexp23(pop2′, pop3′)
    return coex1(pop1′′), coex2(pop2′′), coex3(pop3′′)
end
threeNCoexist = to_hom_expr(FreeBiproductCategory, threeNCoexist)
threeNCoexist_petri = apex(F_tcx(threeNCoexist))

Graph(Petri.Model(threeNCoexist_petri))
display_wd(threeNCoexist)
display_wd(crossexposure)
display_wd(coexist)

save_wd(threeNCoexist, "3gen_coexist_wd.svg")
save_wd(crossexposure, "crossexposure_wd.svg")
save_wd(coexist, "coexist_wd.svg")
save_graph(Graph(Petri.Model(threeNCoexist_petri)), "3gen_coexist_petri.svg")

# save_wd(coexist, "coexist_wd.svg")
tspan = (0.0,150.0)
prob = ODEProblem(vectorfield(threeNCoexist_petri),concentrations(threeNCoexist_petri),tspan,rates(threeNCoexist_petri))
sol = solve(prob,Tsit5())

plot(sol)

# SAVE OLD CODE WHILE WORKING ON IT
# #display_wd(crossexposure)
# # save_wd(crossexposure, "crossexposure_wd.svg")

# # 2 generation cross exposure + coexist model
# population = otimes(S, E, I, I2, A, R, R2, D)
# pop_hom_1 = F(population, 1)
# pop_hom_2 = F(population, 2)
# co_1 = F(coexist, 1)
# co_2 = F(coexist, 2)
# cross_1 = F(crossexposure, 1)
# cross_2 = F(crossexposure, 2)
# twogen = cross ⋅ braid(pop_hom, pop_hom) ⋅ cross ⋅ (co_1 ⊗ co_2)
# test = F(coexist, 1)

# twogen_petri = decoration(twogen)
# twogen_petri = statereplace(twogen_petri, Dict(1=>:S,2=>:E,3=>:I,4=>:I2,5=>:A,6=>:R,7=>:R2,8=>:D,
#                                                9=>:S′,10=>:E′,11=>:I′,12=>:I2′,13=>:A′,14=>:R′,15=>:R2′,16=>:D′))
# tspan = (0.0,150.0)
# prob = ODEProblem(twogen_petri.m,twogen_petri.u0,tspan,twogen_petri.rates)
# sol = solve(prob,Tsit5())

# plot(sol)

# sol_total = transpose(hcat([[s[:S] + s[:S′], s[:E] + s[:E′], s[:A] + s[:A′], s[:I] + s[:I′], s[:I2] + s[:I2′], s[:R] + s[:R′], s[:R2] + s[:R2′], s[:D] + s[:D′]] for s in sol.u]...))
# plot(sol.t, sol_total, label=["S" "E" "A" "I" "I2" "R" "R2" "D"])
# # display_wd(twogen)
# # Graph(decoration(F(twogen)))

# # n generation cross exposure + coexist model
# n = 9
# population = otimes(S, E, I, I2, A, R, R2, D)
# pops = [population for i in 1:n]

# single_exposure = foldl(⋅, [otimes(map(id, pops[1:i])...,crossexposure⋅σ(pops[i+1:i+2]...),map(id, pops[i+3:end])...) for i in 0:(n-2)])
# # display_wd(single_exposure)
# # save_wd(single_exposure, "$(n)gen_single_crossexposure.svg")

# ngen_exposure = foldl(⋅, [single_exposure for i in 1:n])
# # display_wd(ngen_exposure)
# # save_wd(ngen_exposure, "$(n)gen_crossexposure.svg")

# ngen_coexist = ngen_exposure ⋅ foldl(⊗, [coexist for i in 1:n])
# # display_wd(ngen_coexist)
# #save_wd(ngen_coexist, "$(n)gen_coexist_wd.svg")

# ngen_coexist_petri = decoration(F(ngen_coexist, 1))
# #Graph(ngen_coexist_petri)
# #save_graph(Graph(ngen_coexist_petri), "$(n)gen_coexist_petri.svg")

# # Generate some simpler diagrams to expose the hierarchical structure
# @present CoexistOverview(FreeBiproductCategory) begin
#     Pop::Ob
#     coexist::Hom(Pop,Pop)
#     crossexposure::Hom(Pop⊗Pop,Pop⊗Pop)
# end
# pop′,coexist′,crossexposure′ = generators(CoexistOverview)
# pops = [pop′ for i in 1:n]
# single_exposure = foldl(⋅, [otimes(map(id, pops[1:i])...,crossexposure′⋅σ(pops[i+1:i+2]...),map(id, pops[i+3:end])...) for i in 0:(n-2)])
# #display_wd(single_exposure)
# #save_wd(single_exposure, "$(n)gen_single_crossexposure_overview.svg")
# ngen_exposure = foldl(⋅, [single_exposure for i in 1:n])
# #display_wd(ngen_exposure)
# #save_wd(ngen_exposure, "$(n)gen_crossexposure_overview.svg")
# ngen_coexist = ngen_exposure ⋅ foldl(⊗, [coexist′ for i in 1:n])
# #display_wd(ngen_coexist)
# #save_wd(ngen_coexist, "$(n)gen_coexist_overview.svg")

# # Even more generalized diagram
# @present AllCoexist(FreeBiproductCategory) begin
#     TotalPop::Ob
#     AllCoexist::Hom(TotalPop,TotalPop)
#     AllCrossExposure::Hom(TotalPop,TotalPop)
# end
# allpop,allcoexist,allcrossexposure = generators(AllCoexist)
# all_coexist = allcrossexposure ⋅ allcoexist
# #display_wd(all_coexist)
# #save_wd(all_coexist, "all_coexist_wd.svg")


# # Attempt to compute solution to single generation COEXIST model


# # Define time frame and initial parameters
# coexist_pc = F(coexist)
# coexist_petri = decoration(coexist_pc)
# #Graph(coexist_petri)

# tspan = (0.0,150.0)
# u0 = zeros(Float64, base(coexist_pc).n)
# u0[2]  = 13000000

# fatality_rate = 0.146
# β = [0,0,0,0,1/4,.86/.14*.2,1/(10-4),1/15,1/5,1/15,(1/15)*(fatality_rate/(1-fatality_rate))]

# # Generate, solve, and visualize resulting ODE
# prob = ODEProblem(coexist_petri,u0,tspan,β);
# sol = solve(prob,Tsit5());

# plot(sol)

# map(x->x[4], sol.u)

# coexist_petri2 = statereplace(coexist_petri, Dict(1=>:S,2=>:E,3=>:R2,4=>:D,5=>:I2,6=>:A,7=>:R,8=>:I))
# Graph(coexist_petri2)
# u0 = LVector([s=0 for s in coexist_petri2.S]...)
# u0[:S]  = 57345080
# u0[:E]  = 1000
# u0[:I]  = 1000
# u0[:I2] = 1000
# u0[:A]  = 1000
# N = sum(u0)

# β = [.001*social_mixing_rate,0.1*social_mixing_rate,0.6*social_mixing_rate,0.5*social_mixing_rate,1/4,.86/.14*.2,1/(10-4),1/15,1/5,1/15,(1/15)*(fatality_rate/(1-fatality_rate))] # hide

# tspan = (0.0,150.0)
# prob = ODEProblem(coexist_petri2,u0,tspan,β); # hide
# sol = solve(prob,Tsit5()); # hide

# coexist_petri2.Δ
# u0
