using Dates
using Statistics
using JSON

# Define helper function for calculating the root mean squared error
function rmse(sol, data, times, key)
  # Calculate squared differences
  se = sum([(sol(times[i])[Symbol(key)]-data[key][i])^2 for i in 1:length(times)])

  mse = se / length(times)
  sqrt(mse)
end

include("sir_demo.jl")
include("data/make_data.jl")

# Make the parameter estimation directory for the date of param run
folder="param_res/$(Dates.format(Dates.now(), "yyyymmdd_HHMM"))"
mkpath(folder)

# Initialize this simulation (rates don't matter, since we'll use the rates
# from the prediction results
sir = LabelledReactionNet{Number, Number}([:S=>10.0, :I=>1.0, :R=>0.0],
                                          (:inf=>1.0)=>((:S,:I)=>(:I,:I)),
                                          (:rec=>1.0)=>(:I=>:R))
results = Array{Dict{String, Number}, 1}()

# Iterate through possible error to add to generated data
for σ_ind ∈ 1:10
  σ = range(0,1,length=10)[σ_ind]
  println("Starting Sample $σ")
  # Iterate through the number of samples taken for the simulated data
  for n_samps ∈ [5, 10, 50, 100, 500, 1000]
    println("\tStarting $n_samps samples")
    error_S = Array{Float64,1}(undef, 5)
    error_I = Array{Float64,1}(undef, 5)
    error_R = Array{Float64,1}(undef, 5)

    # Iterate through the number of runs make for these parameters
    for run ∈ 1:5
      println("\t\tRun $run")
      sub_folder = "$folder/$(σ)_$n_samps"
      mkpath(sub_folder)

      # Generate data file for this run
      data_loc = "$sub_folder/data_$run.json"
      gen_sir_data(sample_times=range(0,10,length=n_samps), error=σ, dest=data_loc)

      # Get prediction
      pred = run_sir(input=data_loc, folder=sub_folder, plots=false, results="results_$run.json")

      # Read in the generated data file
      j_data = JSON.parsefile(data_loc)
      times = j_data["time_data"]
      tspan = (minimum(times), maximum(times))

      # Simulate using the predicted rates, and calculate the rmse
      prob = ODEProblem(vectorfield(sir), concentrations(sir), tspan, meanRates(pred))
      sol = solve(prob)
      error_S[run] = rmse(sol, j_data["data"], times, "S")
      error_I[run] = rmse(sol, j_data["data"], times, "I")
      error_R[run] = rmse(sol, j_data["data"], times, "R")
    end
    # Add errors to results
    cur_row = Dict{String, Number}("σ"=>σ, "samples"=>n_samps)
    cur_row["error_S"] = mean(error_S)
    cur_row["error_S_std"] = std(error_S)
    cur_row["error_I"] = mean(error_I)
    cur_row["error_I_std"] = std(error_I)
    cur_row["error_R"] = mean(error_R)
    cur_row["error_R_std"] = std(error_R)
    push!(results, cur_row)
  end

  # Write results to json
  open("$folder/summary_pt$σ_ind.json", "w") do f
    write(f, JSON.json(results, 2))
  end
end

open("$folder/summary.json", "w") do f
  write(f, JSON.json(results, 2))
end
