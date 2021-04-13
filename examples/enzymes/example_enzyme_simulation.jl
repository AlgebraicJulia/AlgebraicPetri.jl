# # [Cathepsin Enzyme Reactions](@id enzyme_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/enzymes/example_enzyme_simulation.ipynb)

include("EnzymeReactions.jl")
using .EnzymeReactions

using AlgebraicPetri
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra

using DifferentialEquations
using Plots

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

# Parameter Rates (units of pM, min)
rxns = Dict(
  :K => [inactivate(K, 7.494e-10)
         bindunbind(K, K, 7.814e-4, 3.867e-3)
         degrade(K, K, 2.265e-1)
         bindunbind(K, Kinact, 7.814e-4, 3.867e-3)
         degrade(K, Kinact, 2.265e-1)],
  :S => [inactivate(S, 1.906e-2)
         bindunbind(S, S, 3.534e-7, 1.688e2)
         degrade(S, S, 1.433e1)
         bindunbind(S, Sinact, 3.534e-7, 1.688e2)
         degrade(S, Sinact, 1.433e1)],
  :L => [inactivate(L, 7.810e-3)
         bindunbind(L, L, 1.000e-11, 7.440e3)
         degrade(L, L, 2.670e2)
         bindunbind(L, Linact, 1.000e-11, 7.440e3)
         degrade(L, Linact, 2.670e2)],
  :KE => [bindunbind(K, E, 9.668e-6,1.000e-2)
          degrade(K, E, 1.728e0)],
  :KG => [bindunbind(K, G, 2.764e-6, 8.780e-1)
          degrade(K, G, 1.502e0)],
  :SE => [bindunbind(S, E, 4.197e-7, 1.06e-3)
          degrade(S, E, 1.384e4)],
  :SG => [bindunbind(S, G, 5.152e-8, 3.894e-3)
          degrade(S, G, 8.755e-1)],
  :LE => [bindunbind(L, E, 1.977e-8, 1.000e-2)
          degrade(L, E, 1.066e2)],
  :LG => [bindunbind(L, G, 3.394e-8, 2.365e1)
          degrade(L, G, 4.352e0)],
  :KS => [bindunbind(K, S, 8.822e-4, 4.114e5)
          degrade(K, S, 9.000e-10)
          bindunbind(K, Sinact, 8.822e-4, 4.114e5)
          degrade(K, Sinact, 9.000e-10)],
  :KL => [bindunbind(K, L, 1.756e-4, 3.729e4)
          degrade(K, L, 6.505e6)
          bindunbind(K, Linact, 1.756e-4, 3.729e4)
          degrade(K, Linact, 6.505e6)],
  :SK => [bindunbind(S, K, 3.679e-4, 1.562e3)
          degrade(S, K, 4.410e2)
          bindunbind(S, Kinact, 3.679e-4, 1.562e3)
          degrade(S, Kinact, 4.410e2)],
  :SL => [bindunbind(S, L, 1.000e-3, 5.000e2)
          degrade(S, L, 1.000e-7)
          bindunbind(S, Linact, 1.000e-3, 5.000e2)
          degrade(S, Linact, 1.000e-7)],
  :LK => [bindunbind(L, K, 1.000e-3, 4.118e3)
          degrade(L, K, 3.234e1)
          bindunbind(L, Kinact, 1.000e-3, 4.118e3)
          degrade(L, Kinact, 3.234e1)],
  :LS => [bindunbind(L, S, 1.056e-12, 5.000e2)
          degrade(L, S, 5.000e-1)
          bindunbind(L, Sinact, 1.056e-12, 5.000e2)
          degrade(L, Sinact, 5.000e-1)]
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
prob = ode(KSE, (0.0,120.0))
plot(solve(prob))

KSLE = enzyme_reaction([:K, :S, :L], [:E])
prob = ode(KSLE, (0.0,120.0))
plot(solve(prob))

KSLEG = enzyme_reaction([:K, :S, :L], [:E, :G])
prob = ode(KSLEG, (0.0,120.0))
plot(solve(prob))