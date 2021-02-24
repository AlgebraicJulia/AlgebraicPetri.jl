include("EnzymeReactions.jl")
using .EnzymeReactions
include("Estimators.jl")
using .Estimators

using AlgebraicPetri
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using JSON

using DifferentialEquations
using Distributions
using Turing
using Plots, StatsPlots

#######################
# EDIT CONSTANTS HERE #
#######################

# Initial Concentrations
K = :K=>33000;
S = :S=>33000;
L = :L=>33000;
Kinact = :Kinact=>0;
Sinact = :Sinact=>0;
Linact = :Linact=>0;
E = :E=>700000;
G = :G=>1300000;

# Parameter Distributions
kon = Uniform(log(6e-5), log(6e-2))
koff = Uniform(log(6e-4), log(6e3))
kcat = Uniform(log(1e-10), log(6e4))

# Parameter distributions for each type of reaction
rxns = Dict(
  :K => [inactivate(K, kcat)
         bindunbind(K, K, kon, koff)
         degrade(K, K, kcat)
         bindunbind(K, Kinact, kon, koff)
         degrade(K, Kinact, kcat)],
  :S => [inactivate(S, kcat)
         bindunbind(S, S, kon, koff)
         degrade(S, S, kcat)
         bindunbind(S, Sinact, kon, koff)
         degrade(S, Sinact, kcat)],
  :L => [inactivate(L, kcat)
         bindunbind(L, L, kon, koff)
         degrade(L, L, kcat)
         bindunbind(L, Linact, kon, koff)
         degrade(L, Linact, kcat)],
  :KE => [bindunbind(K, E, kon, koff)
          degrade(K, E, kcat)],
  :KG => [bindunbind(K, G, kon, koff)
          degrade(K, G, kcat)],
  :SE => [bindunbind(S, E, kon, koff)
          degrade(S, E, kcat)],
  :SG => [bindunbind(S, G, kon, koff)
          degrade(S, G, kcat)],
  :LE => [bindunbind(L, E, kon, koff)
          degrade(L, E, kcat)],
  :LG => [bindunbind(L, G, kon, koff)
          degrade(L, G, kcat)],
  :KS => [bindunbind(K, S, kon, koff)
          degrade(K, S, kcat)
          bindunbind(K, Sinact, kon, koff)
          degrade(K, Sinact, kcat)],
  :KL => [bindunbind(K, L, kon, koff)
          degrade(K, L, kcat)
          bindunbind(K, Linact, kon, koff)
          degrade(K, Linact, kcat)],
  :SK => [bindunbind(S, K, kon, koff)
          degrade(S, K, kcat)
          bindunbind(S, Kinact, kon, koff)
          degrade(S, Kinact, kcat)],
  :SL => [bindunbind(S, L, kon, koff)
          degrade(S, L, kcat)
          bindunbind(S, Linact, kon, koff)
          degrade(S, Linact, kcat)],
  :LK => [bindunbind(L, K, kon, koff)
          degrade(L, K, kcat)
          bindunbind(L, Kinact, kon, koff)
          degrade(L, Kinact, kcat)],
  :LS => [bindunbind(L, S, kon, koff)
          degrade(L, S, kcat)
          bindunbind(L, Sinact, kon, koff)
          degrade(L, Sinact, kcat)]
);

# define labels to reaction network mappings
functor(x) = oapply(x, Dict(
  :catK=>enz(rxns, K),
  :catS=>enz(rxns, S),
  :catL=>enz(rxns, L),
  :catKcatS=>enz_enz(rxns, K,S),
  :catKcatL=>enz_enz(rxns, K,L),
  :catScatK=>enz_enz(rxns, S,K),
  :catScatL=>enz_enz(rxns, S,L),
  :catLcatK=>enz_enz(rxns, L,K),
  :catLcatS=>enz_enz(rxns, L,S),
  :catKsubE=>enz_sub(rxns, K,E),
  :catSsubE=>enz_sub(rxns, S,E),
  :catLsubE=>enz_sub(rxns, L,E),
  :catKsubG=>enz_sub(rxns, K,G),
  :catSsubG=>enz_sub(rxns, S,G),
  :catLsubG=>enz_sub(rxns, L,G)));

# helper function to convert undirected wiring diagram to reaction network
enzyme_reaction(args...) = enzyme_uwd(args...) |> functor |> apex

######################
# DEFINE MODELS HERE #
######################

# Define the Enzyme system
KSE = enzyme_reaction([:K, :S], [:E])

# Upload data from .json file
j_data = JSON.parsefile("data/KSE_data.json");

# Perform parameter estimation using Enzyme system and .json files
pred = Estimators.estimate_rates(KSE, j_data, iter_method=PG(100, 100), sample_steps=100,
                                 error_scale=x->log(x+1e-15), param_scale=x->exp(x));
plot(pred)
savefig("pred_results.png")

tspan = (minimum(j_data["time_data"]), maximum(j_data["time_data"]))

# First show plot with the means of the priors
prior_KSE = LabelledReactionNet{Number, Number}(KSE, concentrations(KSE), exp.(mean.(rates(KSE))))
prob = ode(prior_KSE, tspan)
plot(solve(prob))
plot!(j_data["time_data"], [j_data["data"][k] for k in keys(j_data["data"])], seriestype = :scatter)
savefig("prior_model.png")

# Then show the plot with the means of the posteriors
pred_rates = meanRates(pred)
tuned_KSE = LabelledReactionNet{Number, Number}(KSE, concentrations(KSE),
                                                Dict(k=>exp(pred_rates[k]) for k in keys(pred_rates)))
tspan = (minimum(j_data["time_data"]), maximum(j_data["time_data"]))
prob = ode(tuned_KSE, tspan)
plot(solve(prob))
plot!(j_data["time_data"], [j_data["data"][k] for k in keys(j_data["data"])], seriestype = :scatter)
savefig("posterior_model.png")

# Other examples of enzyme reactions
KSLE = enzyme_reaction([:K, :S, :L], [:E])

KSLEG = enzyme_reaction([:K, :S, :L], [:E, :G])