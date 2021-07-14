""" Parameter Estimation with AlgebraicPetri
This script provides an example of how AlgebraicPetri can work well with
other Julia libraries to perform parameter estimation. This specific example
fits a complex enzyme Petri net to data collected from a laboratory
environment.
"""

include("EnzymeReactions.jl")
include("AffinityViz.jl")
include("ParameterEstimation.jl")

using .EnzymeReactions, .ParameterEstimation
using Catlab.CategoricalAlgebra, Catlab.Graphics
using JSON: parsefile, print
using AlgebraicPetri
using Catalyst: ReactionSystem
using OrdinaryDiffEq
using Plots

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger());

gen = enzyme_generators([:K,:S,:L],[:G,:E,:P])
lfunctor(x) = oapply(x, gen);
Graph(gen[:catKsubG])

# Generate the two models we fit to

uwd_kgp = enzyme_uwd([:K], [:G, :P])
model_kgp = uwd_kgp |> lfunctor |> apex;

uwd_kg = enzyme_uwd([:K], [:G])
model_kg = uwd_kg |> lfunctor |> apex;

to_graphviz(uwd_kgp, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"))

# Convert Petri net to ReactionSystem
enzyme_kg = ReactionSystem(model_kg)

""" Define optimization data
We use rates from [doi: 10.1073/pnas.1912207117] as an initial estimate for
fitting to the data.

Additionally, since we are using luminosity data to compute the concentration
of degraded gelatin, we add a constant to scale from luminosity to molarity
"""

init_rates = parsefile("data/KSLGE_rates.json");

G_LUM_CONST = 7.5/64207/39e3/150*1e12/7.1 #pM/LU
kg_data = ModelParams("K_G.json"; t_offset=-7, scale=Dict("G_deg"=>G_LUM_CONST));

res, p1 = optimize(model_kg, enzyme_kg, init_rates, [kg_data]; max_iters = 200, get_md=true, progress=true);

# Plot the results
p = p1
rate_i = 1:nt(model_kg)
t_shift = kg_data.t_offset
t_steps = kg_data.time_data
u0 =  ["$c" ∈ keys(kg_data.concs) ? kg_data.concs["$c"] : 0.0 for c in snames(model_kg)]
prob = ODEProblem(enzyme_kg, u0, (minimum(t_steps) + t_shift,maximum(t_steps)), 10 .^ p)
sol_estimate = solve(prob, Tsit5(), tstops=t_steps .+ t_shift)

plot(t_steps,kg_data.data["G_deg"],seriestype=:scatter,alpha=0.6,
     xlimit=(-7,90), ylabel="Concentration (pM)", markersize=6,
     labels=["G_deg"])
plot!(sol_estimate, linestyle=:dash,lw=4,
      labels=permutedims(String.(snames(model_kg))),
      title="Parameter Estimation Results", xlabel="Time (min)",
      ylimit=(0,1.2e5), legend=:outertopright)
savefig("kg_fit.png")

""" Perform fit on data which includes the Covid-19 spike protien

We will use the rates from the previous fit and use them as assumed values for
fitting to experiments which include the spike protien.
"""

data_arr = Array{ModelParams,1}()
t_offsets = Dict( "K_GS_5.json"=>-6.0,
                  "K_GS_10.json"=>-6.0)

for f in ["K_GS_5.json", "K_GS_10.json"]
  mp = ModelParams(f; t_offset=t_offsets[f], scale=Dict("G_deg"=>G_LUM_CONST))
  push!(data_arr, mp)
end

enzyme_kgp = ReactionSystem(model_kgp; fixed_rates=Dict(tnames(model_kg)[t]=>10 ^ p1[t] for t in 1:nt(model_kg)))

spike_guesses = Dict("deg_KP"=>1e-4,
                     "bind_KP"=>1e-8,
                     "unbind_KP"=>1e-8)
res, p2 = optimize(model_kgp, enzyme_kgp, spike_guesses, data_arr; max_iters = 200, get_md=true, progress=true);

p = p2
ind = 1
data = data_arr[ind]
rate_i = 1:nt(model_kgp)
t_shift = data.t_offset
t_steps = data.time_data
u0 =  ["$c" ∈ keys(data.concs) ? data.concs["$c"] : 0.0 for c in snames(model_kgp)]
prob = ODEProblem(enzyme_kgp, u0, (minimum(t_steps) + t_shift,maximum(t_steps)), 10 .^ p)
sol_estimate = solve(prob,Tsit5(), tstops=t_steps .+ t_shift)
plot(t_steps,data.data["G_deg"],seriestype=:scatter,alpha=0.6, xlimit=(-7,90),
     ylabel="Concentration (pM)", markersize=6, labels=["G_deg"])
plot!(sol_estimate, linestyle=:dash,lw=4,
      labels=permutedims(String.(snames(model_kgp))),
      title="Parameter Estimation Results", xlabel="Time (min)",
      ylimit=(0,1.2e5), legend=:outertopright)

savefig("kgs_fit_res1.png")


ind = 2
data = data_arr[ind]
rate_i = 1:nt(model_kgp)
t_shift = data.t_offset
t_steps = data.time_data
u0 =  ["$c" ∈ keys(data.concs) ? data.concs["$c"] : 0.0 for c in snames(model_kgp)]
prob = ODEProblem(enzyme_kgp, u0, (minimum(t_steps) + t_shift,maximum(t_steps)), 10 .^ p)
sol_estimate = solve(prob, Tsit5(), tstops=t_steps .+ t_shift)
plot(t_steps,data.data["G_deg"],seriestype=:scatter,alpha=0.6, xlimit=(-7,90),
     ylabel="Concentration (pM)", markersize=6, labels=["G_deg"])
plot!(sol_estimate, linestyle=:dash,lw=4,
      labels=permutedims(String.(snames(model_kgp))),
      title="Parameter Estimation Results", xlabel="Time (min)",
      ylimit=(0,1.2e5), legend=:outertopright)

savefig("kgs_fit_res2.png")

open("pred_param.json", "w") do f
  print(f, res, 2)
end
