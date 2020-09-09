# # [Chime Model](@id chime_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/chime.ipynb)

using AlgebraicPetri.Epidemiology

using Petri
using OrdinaryDiffEq
using StochasticDiffEq
using Plots

using Catlab.Meta
using Catlab.Theories
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.Graphics

using JSON

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# help capture JSON of defined functions
macro capture(funcname, exname, ex)
    quote
        $(esc(exname)) = $(JSON.json(strip_lines(ex, recurse=true)))
        $(esc(funcname)) = $ex
    end
end

sir = transmission ⋅ recovery

p_sir = decoration(F_epi(sir));
display_wd(sir)
#-
Graph(p_sir)

policy_days = [20,60,120]
contact_rate = 0.05

u0 = [990.0, 10, 0];
t_span = (17.0,120.0)

@capture γ γ_json 1/14
@capture β β_json t->begin
    pol = findfirst(x->t<=x, policy_days) # array of days when policy changes
    growth_rate = ifelse(pol == 1, 0.0, 2^(1/((pol-1)*5)) - 1) # growth rate depending on policy
    return (growth_rate + γ) / u0[1] * (1-contact_rate) # calculate rate of infection
end
p = [β, γ];

prob = ODEProblem(p_sir,u0,t_span,p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())

plot(sol)

prob,cb = SDEProblem(p_sir,u0,t_span,p);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)

using AlgebraicPetri.Types

cssir = LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, 1)=>((:S, :I)=>(:I,:I)), (:rec, 1)=>(:I=>:R))
csdict = JSON.parse(JSON.json(cssir.tables))
csdict["T"][1]["rate"] = JSON.parse(β_json)
csdict["T"][2]["rate"] = JSON.parse(γ_json)
final_json = JSON.json(csdict)
println(final_json)