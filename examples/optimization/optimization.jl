using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using LabelledArrays
using DifferentialEquations
using Plots

using Catlab
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Programs.RelationalPrograms

using Turing, Distributions
using MCMCChains, StatsPlots
Turing.setadbackend(:forwarddiff)

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

# #### SIR Model:

# define model
sir = PetriNet(3, (1,2)=>(2,2), 2=>3)
Graph(sir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10.0, 1.0, 0.0]
β = [0.4, 0.4]

# The C-Set representation has direct support for generating a DiffEq vector field

prob = ODEProblem(vectorfield(sir),u0,(0.0,10.0),β);
sol = solve(prob,Tsit5(), saveat=0.1)

plot(sol)

σ = 0.2
dt = 0.1
measurements = Array(sol) + σ * randn(size(Array(sol)))
plot(sol, alpha=0.5, legend = false); scatter!(sol.t, measurements')

@model function fit_rates(measurements, prob)
  σ ~ InverseGamma(2,3)

  p ~ product_distribution([Uniform(0,10) for i in 1:length(prob.p)])

  prob = remake(prob,p=p)
  predicted = solve(prob,Tsit5(),saveat=0.1)

  for i = 1:length(predicted)
    measurements[:,i] ~ MvNormal(predicted[i], σ)
  end
end

model = fit_rates(measurements, prob)

sample(model,HMC(0.1,5),1000)

@model function fit_infection(measurements, prob)
  σ ~ InverseGamma(2,3)

  p ~ product_distribution([Uniform(0,10) for i in 1:length(prob.p)])

  prob = remake(prob,p=p)
  predicted = solve(prob,Tsit5(),saveat=0.1)

  for i = 1:length(predicted)
    measurements[i] ~ Normal(predicted[i][2], σ)
  end
end

infected_measurement = measurements[2,:]

model_inf = fit_infection(infected_measurement, prob)
predicted = sample(model_inf,HMC(0.1,5),1000)
pred_arr = Array(predicted)

pl = scatter(sol.t, infected_measurement)

  resol = solve(remake(prob, p=pred_arr[rand(1:size(pred_arr)[1]), 1:length(prob.p)]),Tsit5(),saveat=0.1)
  plot!(resol, alpha=0.1, color = "#BBBBBB", legend = false)
