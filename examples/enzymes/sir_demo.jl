include("Estimators.jl")
using .Estimators

using AlgebraicPetri
using JSON
using Dates

using DifferentialEquations
using Distributions, Turing
using Plots, StatsPlots
using FileIO

function run_sir(;input="data/sir_data.json", 
                 folder="enzyme_res/$(Dates.format(Dates.now(), "yyyymmdd_HHMM"))",
                 plots=true,
                 results="results.jls")
  # Define epidemiological system
  sir = LabelledReactionNet{Distribution, Number}([:S=>10.0, :I=>1.0, :R=>0.0],
                                          (:inf=>truncated(Normal(1,1),0,Inf))=>((:S,:I)=>(:I,:I)),
                                          (:rec=>truncated(Normal(1,1),0,Inf))=>(:I=>:R))
  # Parse input data
  j_data = JSON.parsefile("data/sir_data.json");
  
  # Perform parameter estimation
  #pred = Estimators.estimate_rates(sir, j_data, iter_method=PG(1000, 100), sample_steps=100);
  pred = Estimators.estimate_rates(sir, j_data, iter_method=HMC(0.1, 5), sample_steps=1000);
  
  folder = "enzyme_res/$(Dates.format(Dates.now(), "yyyymmdd_HHMM"))"
  mkpath(folder)
  
  # In order to re-import these results, run:
  # using MCMCChains
  # read("results.jls", Chains)
  write("$folder/$results", pred)

  if(plots)
    mkpath("$folder/prior")
    mkpath("$folder/pred")
    mkpath("$folder/posterior")
    plot(pred)
    savefig("$folder/pred/pred_results.png")
    
    tspan = (minimum(j_data["time_data"]), maximum(j_data["time_data"]))
    times = j_data["time_data"]
    est_times = range(tspan[1], tspan[2], length=100)
    
    # First show plot with the means of the priors
    prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan, mean.(rates(sir)))
    sol = solve(prob)
    for k in keys(j_data["data"])
      plot(est_times, [sol(t)[Symbol(k)] for t in est_times], label="Estimation")
      plot!(j_data["time_data"], j_data["data"][k], label="Measured data", seriestype = :scatter)
      savefig("$folder/prior/$k.png")
    end
    
    # Then show plot with the means of the posteriors
    tuned_SIR = LabelledReactionNet{Number, Number}(sir, concentrations(sir), meanRates(pred))
    tspan = (minimum(j_data["time_data"]), maximum(j_data["time_data"]))
    prob = ODEProblem(vectorfield(tuned_SIR), concentrations(tuned_SIR), tspan, rates(tuned_SIR))
    sol = solve(prob)
    for k in keys(j_data["data"])
      plot(est_times, [sol(t)[Symbol(k)] for t in est_times], label="Estimation")
      plot!(j_data["time_data"], j_data["data"][k], label="Measured data", seriestype = :scatter)
      savefig("$folder/posterior/$k.png")
    end
  end
  pred
end
