```@meta
EditURL = "https://github.com/AlgebraicJulia/AlgebraicPetri.jl/blob/master/docs/examples/covid/epidemiology.jl"
```

# [Basic Epidemiology Models](@id epidemiology_example)

[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](https://nbviewer.jupyter.org/github/AlgebraicJulia/AlgebraicPetri.jl/blob/gh-pages/dev/examples/covid/epidemiology.ipynb)

````@example epidemiology
using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Programs.RelationalPrograms

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));
nothing #hide
````

#### SIR Model:

define model

````@example epidemiology
sir = @relation (s,i,r) begin
    infection(s,i)
    recovery(i,r)
end
display_uwd(sir)
````

````@example epidemiology
p_sir = apex(oapply_epi(sir))
Graph(p_sir)
````

define initial states and transition rates, then
create, solve, and visualize ODE problem

````@example epidemiology
u0 = LVector(S=10, I=1, R=0);
p = LVector(inf=0.4, rec=0.4);
nothing #hide
````

The C-Set representation has direct support for generating a DiffEq vector field

````@example epidemiology
prob = ODEProblem(vectorfield(p_sir),u0,(0.0,7.5),p);
sol = solve(prob,Tsit5())

plot(sol)
````

#### SEIR Model:

define model

````@example epidemiology
seir = @relation (s,e,i,r) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
end
display_uwd(seir)
````

````@example epidemiology
p_seir = apex(oapply_epi(seir))
Graph(p_seir)
````

define initial states and transition rates, then
create, solve, and visualize ODE problem

````@example epidemiology
u0 = LVector(S=10, E=1, I=0, R=0);
p = LVector(exp=.9, ill=.2, rec=.5);

prob = ODEProblem(vectorfield(p_seir),u0,(0.0,15.0),p);
sol = solve(prob,Tsit5())

plot(sol)
````

#### SEIRD Model:

define model

````@example epidemiology
seird = @relation (s,e,i,r,d) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
    death(i,d)
end
display_uwd(seird)
````

````@example epidemiology
p_seird = apex(oapply_epi(seird))
Graph(p_seird)
````

define initial states and transition rates, then
create, solve, and visualize ODE problem

````@example epidemiology
u0 = LVector(S=10, E=1, I=0, R=0, D=0);
p = LVector(exp=0.9, ill=0.2, rec=0.5, death=0.1);

prob = ODEProblem(vectorfield(p_seird),u0,(0.0,15.0),p);
sol = solve(prob,Tsit5())

plot(sol)
````

