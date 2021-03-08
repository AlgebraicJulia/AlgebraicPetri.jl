using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using DifferentialEquations
using JSON

function gen_sir_data(;sample_times=range(0,10,length=100), error=0.2, dest="sir_data.json")
  # define model

  sir = LabelledReactionNet{Number, Number}([:S=>10.0, :I=>1.0, :R=>0.0],
                                                  (:inf=>0.3)=>((:S,:I)=>(:I,:I)),
                                                  (:rec=>0.6)=>(:I=>:R))
 
  tspan = (minimum(sample_times), maximum(sample_times))
  prob = ODEProblem(vectorfield(sir),concentrations(sir),tspan,rates(sir));
  sol = solve(prob,Tsit5(), saveat=sample_times)
  
  # Generate test data
  
  σ = error
  measurements = Array(sol) + σ * randn(size(Array(sol)))
  infected_measurement = measurements[2,:]
  infected_measurement = reshape(infected_measurement, 1, length(infected_measurement)) 
  
  results = Dict("time_data"=>sample_times, "data"=>Dict{Symbol, Array{Float64}}())
  for k in keys(sol(0))
    results["data"][k] = [abs(sol(t)[k] + σ * randn()) for t in sample_times]
  end
  
  open(dest, "w") do f
    write(f, JSON.json(results, 2))
  end
end
