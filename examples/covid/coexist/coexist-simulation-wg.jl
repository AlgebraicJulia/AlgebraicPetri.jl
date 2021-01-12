# # [COEXIST Multi-Generational COVID Model](@id coexist_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/coexist/coexist.ipynb)

using AlgebraicPetri

using OrdinaryDiffEq
using Plots
using JSON

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using Catlab.Theories
using Catlab.WiringDiagrams
using Catlab.Graphics.Graphviz: run_graphviz

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>"1"));
save_fig(g, fname::AbstractString, format::AbstractString) = begin
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, g, format=format)
    end
end

# Define helper functions for defining the two types of
# reactions in an epidemiology Model. Either a state
# spontaneously changes, or one state causes another to change

ob(x::Symbol) = codom(Open([x], LabelledPetriNet(x), [x])).ob;
function spontaneous_petri(transition::Symbol, s::Symbol, t::Symbol)
    Open(LabelledPetriNet(unique([s, t]), transition=>(s=>t)))
end;
function exposure_petri(transition::Symbol, s::Symbol, e::Symbol, t::Symbol)
    Open(LabelledPetriNet(unique([s,e,t]), transition=>((s,e)=>(t,e))))
end;

# Define an `oapply` function that connects the building block Petri nets to
# the operations we will use in the model.

F(ex, n) = oapply(ex, Dict(
    :exposure=>exposure_petri(Symbol(:exp_, n), Symbol(:S,n), Symbol(:I,n), Symbol(:E,n)),
    :exposure_e=>exposure_petri(Symbol(:exp_e, n), Symbol(:S,n), Symbol(:E,n),Symbol(:E,n)),
    :exposure_i2=>exposure_petri(Symbol(:exp_i2, n), Symbol(:S,n), Symbol(:I2,n), Symbol(:E,n)),
    :exposure_a=>exposure_petri(Symbol(:exp_a, n), Symbol(:S,n), Symbol(:A,n),Symbol(:E,n)),
    :progression=>spontaneous_petri(Symbol(:prog_, n), Symbol(:I,n), Symbol(:I2,n)),
    :asymptomatic_infection=>spontaneous_petri(Symbol(:asymp_, n), Symbol(:E,n), Symbol(:A,n)),
    :illness=>spontaneous_petri(Symbol(:ill_, n), Symbol(:E,n), Symbol(:I,n)),
    :asymptomatic_recovery=>spontaneous_petri(Symbol(:arec_, n), Symbol(:A,n), Symbol(:R,n)),
    :recovery=>spontaneous_petri(Symbol(:rec_, n), Symbol(:I2,n), Symbol(:R,n)),
    :recover_late=>spontaneous_petri(Symbol(:rec2_, n), Symbol(:R,n), Symbol(:R2,n)),
    :death=>spontaneous_petri(Symbol(:death2_, n), Symbol(:I2,n), Symbol(:D,n))));

# Define the COEXIST model using the `@relation` macro

coexist = @relation (s, e, i, i2, a, r, r2, d) begin
    exposure(s, i, e)
    exposure_i2(s, i2, e)
    exposure_a(s, a, e)
    exposure_e(s, e, e)
    asymptomatic_infection(e, a)
    asymptomatic_recovery(a, r)
    illness(e, i)
    progression(i, i2)
    death(i2, d)
    recovery(i2, r)
    recover_late(r, r2)
end;
display_uwd(coexist)
Graph(F(coexist, 1))

# Define an `oapply` function that can be used to create a model of
# cross exposure between two sets of populations

F_cx(ex, x,y) = oapply(ex, Dict(
    :exposure=>exposure_petri(Symbol(:exp_, x,y), Symbol(:S,x), Symbol(:I,y), Symbol(:E,x)),
    :exposure_e=>exposure_petri(Symbol(:exp_e, x,y), Symbol(:S,x), Symbol(:E,y),Symbol(:E,x)),
    :exposure_a=>exposure_petri(Symbol(:exp_a, x,y), Symbol(:S,x), Symbol(:A,y),Symbol(:E,x)),
    :exposure_i2=>exposure_petri(Symbol(:exp_i2, x,y), Symbol(:S,x), Symbol(:I2,y), Symbol(:E,x)),
    :exposure′=>exposure_petri(Symbol(:exp_, y,x), Symbol(:S,y), Symbol(:I,x), Symbol(:E,y)),
    :exposure_e′=>exposure_petri(Symbol(:exp_e, y,x), Symbol(:S,y), Symbol(:E,x),Symbol(:E,y)),
    :exposure_a′=>exposure_petri(Symbol(:exp_a, y,x), Symbol(:S,y), Symbol(:A,x),Symbol(:E,y)),
    :exposure_i2′=>exposure_petri(Symbol(:exp_i2, y,x), Symbol(:S,y), Symbol(:I2,x), Symbol(:E,y))
  ),
  Dict(
    :s=>ob(Symbol(:S, x)),
    :e=>ob(Symbol(:E, x)),
    :a=>ob(Symbol(:A, x)),
    :i=>ob(Symbol(:I, x)),
    :i2=>ob(Symbol(:I2, x)),
    :r=>ob(Symbol(:R, x)),
    :r2=>ob(Symbol(:R2, x)),
    :d=>ob(Symbol(:D, x)),
    :s′=>ob(Symbol(:S, y)),
    :e′=>ob(Symbol(:E, y)),
    :a′=>ob(Symbol(:A, y)),
    :i′=>ob(Symbol(:I, y)),
    :i2′=>ob(Symbol(:I2, y)),
    :r′=>ob(Symbol(:R, y)),
    :r2′=>ob(Symbol(:R2, y)),
    :d′=>ob(Symbol(:D, y))
  ));

# Use this new presentation to define a model
# of cross exposure between two populations

crossexposure = @relation (s, e, i, i2, a, r, r2, d, s′, e′, i′, i2′, a′, r′, r2′, d′) begin
    exposure(s, i′, e)
    exposure_i2(s, i2′, e)
    exposure_a(s, a′, e)
    exposure_e(s, e′, e)
    exposure′(s′, i, e′)
    exposure_i2′(s′, i2, e′)
    exposure_a′(s′, a, e′)
    exposure_e′(s′, e, e′)
end;
display_uwd(crossexposure)

# To combine these two models, we need to create a final relational model and
# use the `bundle_legs` function in our `oapply` that enables us to model 3
# population wires instead of each individual state as a wire. Each of these
# populations has their own COEXIST model, and interact through cross exposure

bundled_cross(x,y) = bundle_legs(F_cx(crossexposure, x, y), [tuple([1:8;]...), tuple([9:16;]...)])
bundled_coex(x) = bundle_legs(F(coexist, x), [tuple([1:8;]...)])
F_tcx(ex) = oapply(ex, Dict(
    :crossexp12=>bundled_cross(3,4),
    :crossexp13=>bundled_cross(3,5),
    :crossexp23=>bundled_cross(4,5),
    :coex1=>bundled_coex(3),
    :coex2=>bundled_coex(4),
    :coex3=>bundled_coex(5)));

threeNCoexist = @relation (pop1, pop2, pop3) begin
    crossexp12(pop1, pop2)
    crossexp13(pop1, pop3)
    crossexp23(pop2, pop3)
    coex1(pop1)
    coex2(pop2)
    coex3(pop3)
end;
display_uwd(threeNCoexist)
#-
threeNCoexist_algpetri = apex(F_tcx(threeNCoexist))
Graph(threeNCoexist_algpetri);
save_fig(Graph(threeNCoexist_algpetri), "3ncoexist_petri", "svg"); # hide

# ![3-generation COEXIST model petri net](3ncoexist_petri.svg)

# We can JSON to convert this Petri net into an
# easily shareable format

JSON.print(threeNCoexist_algpetri.tables)

# We can now easily generate a solver for DifferentialEquations.jl
# because we encoded the intitial parameters and rates throughout
# the construction of the model, the final result knows its
# concentrations and rates.

tspan = (0.0,100.0);
prob = ODEProblem(vectorfield(threeNCoexist_algpetri),concentrations(threeNCoexist_algpetri),tspan,rates(threeNCoexist_algpetri));
sol = solve(prob,Tsit5());
plot(sol, xlabel="Time", ylabel="Number of people")

# If we want to model other intervention methods,
# we can simply adjust the rates of exposure to
# represent stay at home orders and mask wearing.
# Because of how we have defined our rates, we can
# simply update the social mixing rates, and
# resolve the model.

for i in 1:length(social_mixing_rate)
  for j in 1:length(social_mixing_rate[1])
    social_mixing_rate[i][j] = social_mixing_rate[i][j] / (i != j ? 10 : 5);
  end
end
threeNCoexist_algpetri = apex(F_tcx(threeNCoexist));

prob = ODEProblem(vectorfield(threeNCoexist_algpetri),concentrations(threeNCoexist_algpetri),tspan,rates(threeNCoexist_algpetri));
sol = solve(prob,Tsit5());
plot(sol, xlabel="Time", ylabel="Number of people")
