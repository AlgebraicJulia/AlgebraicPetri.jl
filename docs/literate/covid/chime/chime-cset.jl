using AlgebraicPetri
using OrdinaryDiffEq
using Plots
using Catlab.Meta
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using JSON

import OrdinaryDiffEq: ODEProblem
ODEProblem(p::LabelledReactionNet, t) = ODEProblem(vectorfield(p), concentrations(p), t, rates(p))

# help capture JSON of defined functions
macro capture(funcname, exname, ex)
    quote
        $(esc(exname)) = $(repr(strip_lines(ex, recurse=true)))
        $(esc(funcname)) = $ex
    end
end

@capture γ γ_text 1/14
@capture β β_text t->begin
    policy_days = [20,60,120] .+ 17
    contact_rate = 0.05
    pol = findfirst(x->t<=x, policy_days) # array of days when policy changes
    growth_rate = pol == 1 ? 0.0 : (2^(1/((pol-1)*5)) - 1) # growth rate depending on policy
    return (growth_rate + γ) / 990 * (1-contact_rate) # calculate rate of infection
end

sir_cset= LabelledReactionNet{Function, Float64}((:S=>990, :I=>10, :R=>0), (:inf, β)=>((:S, :I)=>(:I,:I)), (:rec, t->γ)=>(:I=>:R))

to_graphviz(sir_cset)

prob = ODEProblem(sir_cset, (17.0, 120.0))
sol = OrdinaryDiffEq.solve(prob,Tsit5())
plot(sol)

## Getting Sharable JSON
sir_cset_string = LabelledReactionNet{String, Int}((:S=>990, :I=>10, :R=>0), (:inf, β_text)=>((:S, :I)=>(:I,:I)), (:rec, γ_text)=>(:I=>:R))
JSON.print(tables(sir_cset_string), 2)
