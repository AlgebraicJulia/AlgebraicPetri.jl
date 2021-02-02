using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using LabelledArrays
using DifferentialEquations
using Plots, StatsPlots, MCMCChains
using Turing, Distributions

import Turing: sample

Turing.setadbackend(:forwarddiff)

# API DEVELOPMENT

# ========
# = GOAL =
# ========

# # Define estimation problem and sample
# full_prob = EstimationProblem(sir, [Uniform(0,1) for i ∈ 1:nt(sir)], [:S, :I, :R]=>full_data, u0)
# full_pred = sample(full_prob, HMC(0.01, 5), 1000)

# # Define estimation problem and sample
# inf_prob = EstimationProblem(sir, [Uniform(0,1) for i ∈ 1:nt(sir)], :I=>inf_data, u0)
# inf_pred = sample(inf_prob, HMC(0.01, 5), 1000)

# ============
# = API CODE =
# ============

struct EstimationProblem
  petri::AbstractPetriNet
  tspan
  priors
  data
  u0
end

function EstimationProblem(rxn::Union{AbstractReactionNet,AbstractLabelledReactionNet}, tspan, priors, data)
  EstimationProblem(rxn, tspan, priors, data, concentrations(rxn))
end

function turing_model(prob::EstimationProblem)
  @model function fit_fun(data_idxs, data, prob′)
    σ ~ InverseGamma(2,3)

    # pick distribution with possible range of parameters
    p ~ product_distribution(prob.priors)

    prob′′ = remake(prob′,p=p)
    predicted = solve(prob′′,Tsit5(),saveat=0.1)

    # modify for data format
    for i = 1:length(predicted)
      for (j,s) in enumerate(data_idxs)
        data[j,i] ~ Normal(predicted[i][s], σ)
      end
    end
  end
  fit_fun(prob.data..., ODEProblem(vectorfield(prob.petri), prob.u0, tspan, rates(prob.petri)))
end

sample(prob::EstimationProblem, args...) = sample(turing_model(prob), args...)

# =========
# = USAGE =
# =========

# Generate some fake data

rxn = ReactionNet{Number,Number}([10.0, 1.0, 0.0], (0.3, (1,2)=>(2,2)), (0.6, 2=>3))
tspan = (0.0, 10.0)

prob = ODEProblem(vectorfield(rxn),concentrations(rxn),tspan,rates(rxn));
sol = solve(prob,Tsit5(), saveat=0.1)

σ = 0.2
dt = 0.1
measurements = Array(sol) + σ * randn(size(Array(sol)))
infected_measurement = measurements[2,:]
infected_measurement = reshape(infected_measurement, 1, length(infected_measurement))
recovered_measurement = measurements[3,:]
recovered_measurement = reshape(recovered_measurement, 1, length(recovered_measurement))
scatter(sol.t, infected_measurement', legend=false, xlabel="Time", ylabel="Population")

# Run the estimator

# Fit to just infected data

sir = ReactionNet{Number,Number}([10.0, 1.0, 0.0], (0.3, (1,2)=>(2,2)), (0.6, 2=>3))
Graph(sir)
priors = [Uniform(0,1) for i in 1:nt(sir)]
# est_prob = EstimationProblem(sir, tspan, priors, [2]=>infected_measurement)
est_prob = EstimationProblem(sir, tspan, priors, [2, 3]=>[infected_measurement; recovered_measurement])
pred = sample(est_prob, HMC(0.01, 5), 1000)
β = mean(pred).nt.mean[1:nt(sir)]

plot(pred)

scatter(sol.t, infected_measurement', legend = false)
scatter!(sol.t, recovered_measurement')
ode_prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan,β)
sol = solve(ode_prob, Tsit5(), saveat=0.1)
plot!(sol, legend = false, xlabel="Time", ylabel="Population")

# Fit to data for all states

priors = [Uniform(0,1) for i in 1:nt(sir)]
est_prob = EstimationProblem(sir, tspan, priors, snames(sir)=>measurements)
pred = sample(est_prob, HMC(0.01, 5), 1000)
β = mean(pred).nt.mean[1:nt(sir)]

scatter(sol.t, measurements', legend = false)
ode_prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan,β)
sol = solve(ode_prob, Tsit5(), saveat=0.1)
plot!(sol, legend = false)