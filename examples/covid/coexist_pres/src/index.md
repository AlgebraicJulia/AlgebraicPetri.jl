class: title, no-footer

# AlgebraicPetri: COEXIST COVID-19 Model

## Micah Halter <br/> Evan Patterson <br/> *James Fairbanks* <br/> Georgia Tech Research Institute (GTRI)

### July 23, 2020

---

class: roomy

# Agenda

1. Defining the theory of epidemiology
1. Basic epidemiology models
1. Extend basic epidemiology models
1. Build and simulate COEXIST COVID-19 model

---

```@setup coexist
using LabelledArrays;
using AlgebraicPetri;
using AlgebraicPetri.Epidemiology;
using Petri;
using DifferentialEquations;
using Plots;
using Catlab;
using Catlab.Theories;
using Catlab.Programs;
using Catlab.CategoricalAlgebra.ShapeDiagrams;
using Catlab.WiringDiagrams;
using Catlab.Graphics;
using Catlab.Graphics.Graphviz: run_graphviz;
⊗ = otimes;
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);
mkdir("figs")
save_wd(ex, fname::AbstractString) = begin
    g = display_wd(ex)
    open(joinpath("figs", fname), "w") do io
        run_graphviz(io, g, format="svg")
    end
end;
save_graph(g, fname::AbstractString) = begin
    open(joinpath("figs", fname), "w") do io
        run_graphviz(io, g, format="svg")
    end
end;
dictreplace(f::Dict{K,N}, d::Dict{K,V}) where {K,N,V} = Dict{N,V}(Base.map(x->Pair(f[x[1]], x[2]), collect(d)))
dictreplace(f::Dict, ts::Vector{Tuple{S,T}}) where {S<:Dict, T<:Dict} = [(dictreplace(f, t[1]), dictreplace(f, t[2])) for t in ts]
statereplace(m::Model, f::Dict{K,N}) where {K,N} = Petri.Model(map(s->f[s], m.S), dictreplace(f, m.Δ))

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

save_graph(Graph(decoration(F_epi(illness))), "illness_petri.svg");
save_wd(illness, "illness_wd.svg");
save_graph(Graph(decoration(F_epi(recovery))), "recovery_petri.svg");
save_wd(recovery, "recovery_wd.svg");
save_graph(Graph(decoration(F_epi(transmission))), "transmission_petri.svg");
save_wd(transmission, "transmission_wd.svg");
save_graph(Graph(decoration(F_epi(exposure))), "exposure_petri.svg");
save_wd(exposure, "exposure_wd.svg");
```

# Defining the Theory of Epidemiology

```julia
@present InfectiousDiseases(FreeBiproductCategory) begin
    S::Ob
    E::Ob
    I::Ob
    R::Ob
    D::Ob
    transmission::Hom(S⊗I, I)
    exposure::Hom(S⊗I, E)
    illness::Hom(E,I)
    recovery::Hom(I,R)
    death::Hom(I,D)
end
```

---

class: compact

# Epidemiology Building Blocks

|    Algebraic Expression   |         Wiring Diagram         |                     Petri Net                    |
|:-------------------------:|:------------------------------:|:------------------------------------------------:|
|      $illness: E → I$     | ![](figs/illness_wd.svg# mh-5) |    ![](figs/illness_petri.svg# mh-5)   |
| $transmission: S ⊗ I → I$ |  ![](figs/transmission_wd.svg# mh-5) | ![](figs/transmission_petri.svg# mh-5) |
|   $exposure: S ⊗ I → E$   |    ![](figs/exposure_wd.svg# mh-5)   |   ![](figs/exposure_petri.svg# mh-5)   |
|     $recovery: I → R$     |    ![](figs/recovery_wd.svg# mh-5)   |   ![](figs/recovery_petri.svg# mh-5)   |

---

# Basic SIR Model

```@setup coexist
sir = transmission ⋅ recovery
sir_petri = statereplace(decoration(F_epi(sir)), Dict(1=>:S,2=>:I,3=>:R))
save_graph(Graph(sir_petri), "sir_petri.svg")
save_wd(sir, "sir_wd.svg")
```

.center[
$$sir = transmission ⋅ recovery$$

<br/>
![](figs/sir_wd.svg)

![](figs/sir_petri.svg# pt-4)
]

---

class: compact

# SEIR Model

```@example coexist
seir = @program InfectiousDiseases (s::S,i::I) begin
  e = exposure(s,i)
  i2 = illness(e)
  i_all = [i,i2]
  return recovery(i_all)
end
nothing # hide
```

```@setup coexist
seir =  to_hom_expr(FreeBiproductCategory, seir)
seir_petri = statereplace(decoration(F_epi(seir)), Dict(1=>S,2=>E,3=>I,4=>R))
save_graph(Graph(seir_petri), "seir_petri.svg")
save_wd(seir, "seir_wd.svg")
```

![](figs/seir_wd.svg# center)

![](figs/seir_petri.svg# center)

---

# COEXIST COVID-19 Model

![COEXIST Model](assets/coexist.png# maxw-80pct center pt-4)

---

# Extend Basic Epidemiology

```julia
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
end
```

---

class: compact

# Defining COEXIST

```@example coexist
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
nothing # hide
```

```@setup coexist
save_wd(coexist, "coexist_wd.svg")
coexist =  to_hom_expr(FreeBiproductCategory, coexist)
coexist_petri = decoration(F(coexist))
coexist_petri = statereplace(coexist_petri, Dict(1=>:S,2=>:E,3=>:R2,4=>:D,5=>:I2,6=>:A,7=>:R,8=>:I))
save_graph(Graph(coexist_petri), "coexist_petri.svg")
```

---

# COEXIST SEIRD Model Petri Net

![COEXIST Petri Net](figs/coexist_petri.svg# center)

---

class: compact

# Simulate the Model

```@example coexist
u0 = LVector(;zip(coexist_petri.S, zeros(length(coexist_petri.S)))...) # hide
u0.S  = 57345080
u0.E  = 1000
u0.I  = 1000
u0.I2 = 1000
u0.A  = 1000
nothing # hide
```

```@setup coexist
tspan = (0.0,200.0)
social_mixing_rate = 1.2232/sum(u0)
social_mixing_rate_with_distancing = 0.7748/sum(u0)
fatality_rate = 0.146
β = [.001*social_mixing_rate,0.1*social_mixing_rate,0.6*social_mixing_rate,0.5*social_mixing_rate,1/4,.86/.14*.2,1/(10-4),1/15,1/5,1/15,(1/15)*(fatality_rate/(1-fatality_rate))]

prob = ODEProblem(coexist_petri,u0,tspan,β)
sol = solve(prob,Tsit5())

plot(sol)
savefig(joinpath("figs","sol_plot.svg"))
```

![](figs/sol_plot.svg# center)

---

class: compact

# Define Inter-Generational Cross Exposure

```@example coexist
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
nothing # hide
```

```@setup coexist
save_wd(crossexposure, "crossexposure_wd.svg")
crossexposure = to_hom_expr(FreeBiproductCategory, crossexposure)

@present CoexistOverview(FreeBiproductCategory) begin
    Pop::Ob
    coexist::Hom(Pop,Pop)
    crossexposure::Hom(Pop⊗Pop,Pop⊗Pop)
end
pop′,coexist′,crossexposure′ = generators(CoexistOverview)
twogen′ = crossexposure′ ⋅ σ(pop′, pop′) ⋅ crossexposure′ ⋅ (coexist′ ⊗ coexist′)
save_wd(twogen′, "twogen_overview_wd.svg")
```

![](figs/twogen_overview_wd.svg# center pt-5)

---

```@setup coexist
population = otimes(S, E, I, I2, A, R, R2, D)
twogen = crossexposure ⋅ σ(population, population) ⋅ crossexposure ⋅ (coexist ⊗ coexist)
twogen_petri = decoration(F(twogen))
twogen_petri = statereplace(twogen_petri, Dict(1=>:S,2=>:E,3=>:I,4=>:I2,5=>:A,6=>:R,7=>:R2,8=>:D,
                                               9=>:S′,10=>:E′,11=>:I′,12=>:I2′,13=>:A′,14=>:R′,15=>:R2′,16=>:D′))
save_wd(twogen, "twogen_wd.svg")
save_graph(Graph(twogen_petri), "twogen_petri.svg")
```

![Two Generation Petri Net](figs/twogen_petri.svg# maxw-80pct center)

---

# Petri Nets with Rates

![Formalization of Petri Nets with Rates](assets/petri_with_rates.png# maxw-50pct center pt-5)

---

# Undirected Wiring Diagrams

![](figs/twogen_overview_wd.svg# center)

![Two Generation Undirected Wiring Diagram](assets/twogen_coexist_uwd.png# maxw-50pct center)

---

# $\mathcal{n}$-Generational COEXIST Model

![Three Generation Undirected Wiring Diagram](assets/threegen_coexist_uwd.png# maxw-50pct center)

---

**Coexist:**

![Coexist Undirected Wiring Diagram](assets/coexist_uwd.png# maxw-50pct center)

**Cross Exposure:**

![Cross Exposure Undirected Wiring Diagram](assets/crossexposure_uwd.png# maxw-50pct center)

---

class: compact

# Relation Syntax

```julia
coexist = @relation (s::S, e::E, i::I, i2::I2, a::A, r::R, r2::R2, d::D) begin
    exposure(s, i, e)
    exposure(s, i2, e)
    exposure(s, a, e)
    exposure(s, e, e)
    illness(e, a)
    illness(e, i)
    progression(i, i2)
    death(i2, d)
    recovery(a, r)
    recovery(i2, r)
    progression(r, r2)
end
```

```julia
crossexposure = @relation (s::S, e::E, i::I, i2::I2, a::A, r::R, r2::R2, d::D,
                           s′::S, e′::E, i′::I, i2′::I2, a′::A, r′::R, r2′::R2, d′::D) begin
    exposure(s, i′, e)
    exposure(s, i2′, e)
    exposure(s, a′, e)
    exposure(s, e′, e)
end
```
