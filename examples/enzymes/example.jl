include("Estimators.jl")
using .Estimators
using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using DifferentialEquations
using Distributions
using Plots, StatsPlots

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
pred = estimate_rates(sir, tspan, priors, [2, 3]=>[infected_measurement; recovered_measurement])
β = mean(pred).nt.mean[1:nt(sir)]

plot(pred)

scatter(sol.t, infected_measurement', legend = false)
scatter!(sol.t, recovered_measurement')
ode_prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan,β)
sol = solve(ode_prob, Tsit5(), saveat=0.1)
plot!(sol, legend = false, xlabel="Time", ylabel="Population")

# Fit to data for all states

priors = [Uniform(0,1) for i in 1:nt(sir)]
pred = estimate_rates(sir, tspan, priors, snames(sir)=>measurements)
β = mean(pred).nt.mean[1:nt(sir)]

scatter(sol.t, measurements', legend = false)
ode_prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan,β)
sol = solve(ode_prob, Tsit5(), saveat=0.1)
plot!(sol, legend = false)