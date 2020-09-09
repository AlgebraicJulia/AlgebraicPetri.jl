# # [Chime Model](@id chime_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/chime.ipynb)

using AlgebraicPetri.Epidemiology

using Petri
using OrdinaryDiffEq
using StochasticDiffEq
using Plots

using Catlab.Theories
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.Graphics

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

sir = transmission ⋅ recovery

p_sir = decoration(F_epi(sir));
display_wd(sir)
#-
Graph(p_sir)

β(S,γ,contact_rate,policy_days,t) = begin
    pol = findfirst(x->t<=x, policy_days) # array of days when policy changes
    growth_rate = ifelse(pol == 1, 0.0, 2^(1/((pol-1)*5)) - 1) # growth rate depending on policy
    return (growth_rate + γ) / S * (1-contact_rate) # calculate rate of infection
end
γ = 1/14

policy_days = [20,60,120]
contact_rate = 0.05

u0 = [990.0, 10, 0];
p = [t->β(u0[1],γ,contact_rate,policy_days,t), γ];
t_span = (17.0,120.0)

prob = ODEProblem(p_sir,u0,t_span,p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())

plot(sol)

prob,cb = SDEProblem(p_sir,u0,t_span,p);
sol = solve(prob,SRA1(),callback=cb)

plot(sol)