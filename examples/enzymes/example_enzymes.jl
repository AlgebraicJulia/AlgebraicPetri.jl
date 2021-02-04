include("EnzymeReactions.jl")
using .EnzymeReactions
include("Estimators.jl")
using .Estimators

using AlgebraicPetri
using Catlab.WiringDiagrams

using DifferentialEquations
using Distributions
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
kon = Uniform(6e-5, 6e-2)
koff = Uniform(6e-4, 6e3)
kcat = Uniform(0, 6e4)

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

KSE = enzyme_reaction([:K, :S], [:E])

KSLE = enzyme_reaction([:K, :S, :L], [:E])

KSLEG = enzyme_reaction([:K, :S, :L], [:E, :G])