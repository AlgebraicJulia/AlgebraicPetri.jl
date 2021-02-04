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

sir = LabelledReactionNet{Distribution,Number}([:S=>10.0, :I=>1.0, :R=>0.0], (:inf, Uniform(0,1))=>((:S,:I)=>(:I,:I)), (:rec, Uniform(0,1))=>(:I=>:R));
Graph(sir)
pred = estimate_rates(sir, tspan, [:I, :R]=>[infected_measurement; recovered_measurement])
β = Dict(tname(sir,t)=>mean(pred).nt.mean[t] for t in 1:nt(sir))
plot(pred)

scatter(sol.t, infected_measurement', legend = false)
scatter!(sol.t, recovered_measurement')
ode_prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan,β)
sol = solve(ode_prob, Tsit5(), saveat=0.1)
plot!(sol, legend = false, xlabel="Time", ylabel="Population")

# Fit to data for all states

pred = estimate_rates(sir, tspan, snames(sir)=>measurements)
β = Dict(tname(sir,t)=>mean(pred).nt.mean[t] for t in 1:nt(sir))

scatter(sol.t, measurements', legend = false)
ode_prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan,β)
sol = solve(ode_prob, Tsit5(), saveat=0.1)
plot!(sol, legend = false)