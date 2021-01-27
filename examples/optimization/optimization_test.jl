using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using LabelledArrays
using DifferentialEquations
using Plots, StatsPlots, MCMCChains
using Turing, Distributions
Turing.setadbackend(:forwarddiff)

# define model

sir = PetriNet(3, (1,2)=>(2,2), 2=>3)
Graph(sir)

u0 = [10.0, 1.0, 0.0]
β = [.3, .6]

prob = ODEProblem(vectorfield(sir),u0,(0.0,10.0),β);
sol = solve(prob,Tsit5(), saveat=0.1)

plot(sol)

# Generate test data

σ = 0.2
dt = 0.1
measurements = Array(sol) + σ * randn(size(Array(sol)))
infected_measurement = measurements[2,:]
infected_measurement = reshape(infected_measurement, 1, length(infected_measurement))
scatter!(sol.t, measurements')

# Fit to mock data

@model function fit_rates(measurements, prob′)
  σ ~ InverseGamma(2,3)

  p ~ product_distribution([Uniform(0,1) for i in 1:length(prob′.p)])

  prob′′ = remake(prob′,p=p)
  predicted = solve(prob′′,Tsit5(),saveat=0.1)

  for i = 1:length(predicted)
    measurements[:,i] ~ MvNormal(predicted[i], σ)
  end
end

model = fit_rates(measurements, prob)

predicted = sample(model,HMC(0.01,5),1000)
pred_arr = Array(predicted)

scatter(sol.t, measurements', legend = false)
resol = solve(remake(prob,p=pred_arr[rand(1:(size(pred_arr)[1])), 1:2]), Tsit5(), saveat=0.1)
plot!(resol, legend = false)

# Fit to just infected data

@model function fit_infection(measurements, prob′)
  σ ~ InverseGamma(2,3)

  # pick distribution with possible range of parameters
  p ~ product_distribution([Uniform(0,1) for i in 1:length(prob′.p)])
  # p ~ product_distribution([truncated(Normal(.5,.5), 0, 1) for i in 1:length(prob′.p)])

  prob′′ = remake(prob′,p=p)
  predicted = solve(prob′′,Tsit5(),saveat=0.1)

  for i = 1:length(predicted)
    measurements[i] ~ Normal(predicted[i][2], σ)
  end
end

model_inf = fit_infection(infected_measurement, prob)
# Set gradient size
predicted = sample(model_inf,HMC(0.01,5),1000)
pred_arr = Array(predicted)

scatter(sol.t, infected_measurement', legend = false)
resol = solve(remake(prob,p=pred_arr[rand(1:(size(pred_arr)[1])), 1:2]), Tsit5(), saveat=0.1)
plot!(resol, legend = false)