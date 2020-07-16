# # [Multi-City COVID-19 Model](@id covid_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/covid.ipynb)

using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using Petri
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

generators(EpiCoexist)

S,E,I,R,D,I2,A,R2,transmission,exposure,illness,recovery,death,exposure_e,exposure_i2,exposure_a,progression,asymptomatic_infection,recover_late,asymptomatic_recovery,recovery2, death2 = generators(EpiCoexist);
new_functor = copy(FunctorGenerators)
new_functor[I2] = new_functor[I]
new_functor[A] = new_functor[E]
new_functor[R2] = new_functor[R]
new_functor[exposure_e] = new_functor[exposure]
new_functor[exposure_i2] = new_functor[exposure]
new_functor[exposure_a] = new_functor[exposure]
new_functor[progression] = new_functor[illness]
new_functor[asymptomatic_infection] = new_functor[illness]
new_functor[asymptomatic_recovery] = new_functor[recovery]
new_functor[recover_late] = new_functor[recovery]
new_functor[recovery2] = new_functor[recovery]
new_functor[death2] = new_functor[death]

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=new_functor);

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
display_wd(coexist)
save_wd(coexist, "coexist_wd.svg")
coexist = to_hom_expr(FreeBiproductCategory, coexist)

crossexposure = @program EpiCoexist (s::S, e::E, i::I, i2::I2, a::A, r::R, r2::R2, d::D,
                               s′::S, e′::E, i′::I, i2′::I2, a′::A, r′::R, r2′::R2, d′::D) begin
    e_2 = exposure(s, i′)
    e_3 = exposure_i2(s, i2′)
    e_4 = exposure_a(s, a′)
    e_5 = exposure_e(s, e′)
    e_all = [e, e_2, e_3, e_4, e_5]
    return s, e_all, i, i2, a, r, r2, d,
           s′, e′, i′, i2′, a′, r′, r2′, d′
    # braid the output to easily recurse
end
display_wd(crossexposure)
save_wd(crossexposure, "crossexposure_wd.svg")
crossexposure = to_hom_expr(FreeBiproductCategory, crossexposure)

# 2 generation cross exposure + coexist model
twogen = (coexist ⊗ coexist) ⋅ crossexposure
display_wd(twogen)
Graph(decoration(F(twogen)))

# n generation cross exposure + coexist model
n = 5
population = otimes(S, E, I, I2, A, R, R2, D)
pops = [population for i in 1:n]

single_exposure = foldl(⋅, [otimes(map(id, pops[1:i])...,crossexposure⋅σ(pops[i+1:i+2]...),map(id, pops[i+3:end])...) for i in 0:(n-2)])
display_wd(single_exposure)
save_wd(single_exposure, "single_crossexposure.svg")

ngen_exposure = foldl(⋅, [single_exposure for i in 1:n])
display_wd(ngen_exposure)
save_wd(ngen_exposure, "ngen_crossexposure.svg")

ngen_coexist = ngen_exposure ⋅ foldl(⊗, [coexist for i in 1:n])
display_wd(ngen_coexist)
save_wd(ngen_coexist, "ngen_coexist_wd.svg")

save_graph(Graph(decoration(F(ngen_coexist))), "ngen_coexist_petri.svg")

# Generate some simpler diagrams to expose the hierarchical structure
@present CoexistOverview(FreeBiproductCategory) begin
    Pop::Ob
    coexist::Hom(Pop,Pop)
    crossexposure::Hom(Pop⊗Pop,Pop⊗Pop)
end
pop′,coexist′,crossexposure′ = generators(CoexistOverview)
n = 5
pops = [pop′ for i in 1:n]
single_exposure = foldl(⋅, [otimes(map(id, pops[1:i])...,crossexposure′⋅σ(pops[i+1:i+2]...),map(id, pops[i+3:end])...) for i in 0:(n-2)])
display_wd(single_exposure)
save_wd(single_exposure, "single_crossexposure_overview.svg")
ngen_exposure = foldl(⋅, [single_exposure for i in 1:n])
display_wd(ngen_exposure)
save_wd(ngen_exposure, "ngen_crossexposure_overview.svg")
ngen_coexist = ngen_exposure ⋅ foldl(⊗, [coexist′ for i in 1:n])
display_wd(ngen_coexist)
save_wd(ngen_coexist, "ngen_coexist_overview.svg")

# Even more generalized diagram
@present AllCoexist(FreeBiproductCategory) begin
    TotalPop::Ob
    AllCoexist::Hom(TotalPop,TotalPop)
    AllCrossExposure::Hom(TotalPop,TotalPop)
end
allpop,allcoexist,allcrossexposure = generators(AllCoexist)
all_coexist = allcrossexposure ⋅ allcoexist
display_wd(all_coexist)
save_wd(all_coexist, "all_coexist_wd.svg")


# Attempt to compute solution to single generation COEXIST model

# identify states and transitions
# S_1 = S
# S_2 = E
# S_3 = R2
# S_4 = D
# S_5 = I2
# S_6 = A
# S_7 = R1
# S_8 = I1

# T_1 = exposure_e
# T_2 = exposure_i1
# T_3 = exposure_i2
# T_4 = exposure_a
# T_5 = progression = 1/4
# T_6 = asymptomatic_infection = 0.86/0.14 * .2
# T_7 = recovery2 = 1/(10-4)
# T_8 = asymptomatic_recovery = 1/15
# T_9 = illness = 1/5
# T_10 = recover_late = 1/15
# T_11 = death2 = recovery2 * fatality_hospital_ratio / 1 - fatality_hospital_ratio

# Define time frame and initial parameters
coexist_pc = F(coexist)
coexist_petri = decoration(coexist_pc)
Graph(coexist_petri)

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