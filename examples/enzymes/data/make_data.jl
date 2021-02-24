using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using DifferentialEquations
using Plots, StatsPlots, MCMCChains
using Turing, Distributions
using JSON

# define model

sir = LabelledReactionNet{Number, Number}([:S=>10.0, :I=>1.0, :R=>0.0],
                                                (:inf=>0.3)=>((:S,:I)=>(:I,:I)),
                                                (:rec=>0.6)=>(:I=>:R))

prob = ODEProblem(vectorfield(sir),concentrations(sir),(0.0,10.0),rates(sir));
sol = solve(prob,Tsit5(), saveat=0.1)

plot(sol)

# Generate test data

σ = 0.2
dt = 0.1
measurements = Array(sol) + σ * randn(size(Array(sol)))
infected_measurement = measurements[2,:]
infected_measurement = reshape(infected_measurement, 1, length(infected_measurement))

times = range(0, 10, length=100)

results = Dict("time_data"=>times, "data"=>Dict{Symbol, Array{Float64}}())
for k in keys(sol(0))
  results["data"][k] = [abs(sol(t)[k] + σ * randn()) for t in times]
end

open("sir_data.json", "w") do f
  write(f, JSON.json(results, 2))
end
