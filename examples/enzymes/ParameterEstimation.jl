module ParameterEstimation
  using JSON: parsefile
  using Catalyst: reactions, Sym, ReactionSystem, params
  using DifferentialEquations
  using AlgebraicPetri
  using DiffEqFlux, Flux, OrdinaryDiffEq

  export ModelParams, optimize

  struct ModelParams
    concs::Dict
    data::Dict
    time_data::Array{<:Real,1}
    t_offset::Real
  end

  function ModelParams(fname::String; t_offset=0, scale=Dict{Symbol, Float64}())
    filedata = parsefile(fname)
    samples = filedata["data"]
    t_steps = float.(filedata["time_data"])
    concs = filedata["init_conc"]
    for k in keys(scale)
      samples[k] .*= scale[k]
    end
    ModelParams(concs, samples, t_steps, t_offset)
  end

	""" TODO:
	Add tooling for estimating concentrations
	Add tooling for estimating start times
	Add tooling for combining multiple experiments with multiple intersecting parameters

	Allow for a range to be specified for each value, providing a linear ramp-up of loss if the parameter is outside of that range.
	"""
	function optimize(model, rs::ReactionSystem, p_init, data::Array{ModelParams};
                          max_iters=1000,
                          p_ranges = fill((-Inf, Inf), length(p_init)),
                          l_factor=10, get_md=false, kw...)

    prob = ODEProblem(rs, zeros(ns(model)), (0.0,100.0), zeros(length(params(rs))))
    sym_rates = filter(r->(reactions(rs)[r].rate isa Sym), 1:nt(model))
    rate_var = tnames(model)[sym_rates]
    rate_ind = 1:length(rate_var)
    p_init_arr = [log10(p_init["$r"]) for r in rate_var]

	  function loss(p)
		loss = 0
		for param in data
      u0 = ["$c" ∈ keys(param.concs) ? param.concs["$c"] : 0.0 for c in snames(model)]
			t_steps = param.time_data
			t_off = param.t_offset
			samples = param.data
			sol = solve(remake(prob, u0 = u0,
						tspan=(minimum(t_steps) + t_off,maximum(t_steps)),
						p=  10.0 .^ (p[rate_ind])),
					Tsit5(), saveat=t_steps)
			for k in keys(samples)
				vals = [sol.u[i][findfirst(==(Symbol(k)), snames(model))] for i in 1:length(t_steps)]
				loss += sum((vals .- samples[k]).^2 ./(samples[k].^2))
			end
			for i in 1:length(p)
				cur_p = p[i]
				p_range = p_ranges[i]
				if cur_p > p_range[2]
					loss += l_factor * (cur_p - p_range[2])
				elseif cur_p < p_range[1]
					loss += l_factor * (p_range[1] - cur_p)
				end
			end
	    end
	    return loss#, sol
	  end
	  res = DiffEqFlux.sciml_train(loss,p_init_arr,ADAM(0.1);maxiters = max_iters, kw...)
    real_p = Dict(rate_var[r] => 10.0 ^ res[r] for r in rate_ind)
    for t in 1:nt(model)
      if !(t ∈ sym_rates)
        real_p[tname(model, t)] = reactions(rs)[t].rate
      end
    end

    if get_md
      return real_p, res
    else
      return real_p
    end
	end
end
