include("Estimators.jl")
using .Estimators

using AlgebraicPetri
using JSON

using DifferentialEquations
using Distributions, Turing
using Plots, StatsPlots

# Define epidemiological system
sir = LabelledReactionNet{Distribution, Number}([:S=>10.0, :I=>1.0, :R=>0.0],
                                                (:inf=>truncated(Normal(1,1),0,Inf))=>((:S,:I)=>(:I,:I)),
                                        (:rec=>truncated(Normal(1,1),0,Inf))=>(:I=>:R))
# Parse input data
j_data = JSON.parsefile("data/sir_data.json");

# Perform parameter estimation
pred = Estimators.estimate_rates(sir, j_data, iter_method=PG(100, 100), sample_steps=100);
plot(pred)
savefig("pred_results.png")

tspan = (minimum(j_data["time_data"]), maximum(j_data["time_data"]))

# First show plot with the means of the priors
prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan, mean.(rates(sir)))
plot(solve(prob))
plot!(j_data["time_data"], [j_data["data"][k] for k in keys(j_data["data"])], seriestype = :scatter)
savefig("prior_model.png")

# Then show plot with the means of the posteriors
tuned_SIR = LabelledReactionNet{Number, Number}(sir, concentrations(sir), meanRates(pred))
tspan = (minimum(j_data["time_data"]), maximum(j_data["time_data"]))
prob = ODEProblem(vectorfield(tuned_SIR), concentrations(tuned_SIR), tspan, rates(tuned_SIR))
plot(solve(prob))
plot!(j_data["time_data"], [j_data["data"][k] for k in keys(j_data["data"])], seriestype = :scatter)
savefig("posterior_model.png")
