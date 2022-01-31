module Sensitivity
using ForwardDiff
using AlgebraicPetri
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using DifferentialEquations
using PerceptualColourMaps
using Colors
using StatsBase
using LabelledArrays
using Base.Iterators: flatten
using JSON

export max_metric, int_metric, final_metric, sensitivity, GraphHeatmap

function max_metric(model, concs; t_range=(0.0,100.0))
  prob = ODEProblem(vectorfield(model), concentrations(model), t_range, rates(model))
  metric(p) = begin
    prob = remake(prob; p=p)
    sol = solve(prob)
    maximum(sum(sol(t)[k] for k in concs) for t in sol.t)
  end
end


function final_metric(model, concs; t_range=(0.0,100.0))
  prob = ODEProblem(vectorfield(model), concentrations(model), t_range, rates(model))
  metric(p) = begin
    prob = remake(prob; p=p)
    sol = solve(prob)
    sum(sol(t_range[2])[k] for k in concs)
  end
end

function int_metric(model, concs; t_range=(0,100.0), kw...)
  rs = rates(model)
  prob = ODEProblem(vectorfield(model; kw...), concentrations(model), t_range, rs)
  metric(p) = begin
    prob = remake(prob; p=p)
    sol = solve(prob)
    prob2 = ODEProblem((du,u,p,t)->(du .= sol(t)), sol(0), t_range, [])
    sol2 = solve(prob2)
    sum(sol2(t_range[2])[k] for k in concs)
  end
end

function sensitivity(metric, rates)
  ForwardDiff.gradient(metric, rates)
end

function f2range(v, clims, N)
    v = (max(clims[1], min(clims[2], v)) - clims[1]) / (clims[2] - clims[1])
    round(Int64, (N-1)*v + 1)
end

function to_position(val)
  isnothing(val) && return ""
  "$(val[1]),$(val[2])!"
end

function GraphHeatmap(p::AbstractPetriNet, rates; clims=nothing,positions=Dict(:T=>fill(nothing, nt(p)), :S=>fill(nothing, ns(p))))

  colors = string.(hex.(cmap("D01")))
  clims = isnothing(clims) ? (minimum(rates), maximum(rates)) : clims

  statenodes = [Node("s$s", Attributes(:label=>"$(sname(p, s))",:shape=>"circle", :color=>"#6C9AC3",
                                       :pos=>to_position(positions[:S][s]))) for s in 1:ns(p)]

  transnodes = [Node("t$k", Attributes(:label=>"$(tname(p, k))", :shape=>"square", :color=>"#$(colors[f2range(rates[tname(p,k)], clims, 256)][3:end])",
                                       :pos=>to_position(positions[:T][k]))) for k in 1:nt(p)]

  graph_attrs = Attributes(:rankdir=>"LR")
  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  stmts = vcat(statenodes, transnodes)

  edges = map(1:nt(p)) do k
    vcat(edgify(countmap_wrap(inputs(p, k)), k, false),
         edgify(countmap_wrap(outputs(p, k)), k, true))
  end |> flatten |> collect

  stmts = vcat(stmts, edges)
  g = Graphviz.Graph(;name="G", stmts=stmts, directed=true,
                     prog= all(isnothing.(vcat(positions[:T], positions[:S]))) ? "dot" : "fdp",
                     graph_attrs=graph_attrs, node_attrs=node_attrs,
                     edge_attrs=edge_attrs)
  return g
end
countmap_wrap(a) = isempty(a) ? Dict{Int, Int}() : countmap(a)

function edgify(δ::Dict{Int64, Int64}, transition, reverse::Bool; pre="")
  return [Edge(reverse ? ["\"$(pre)t$transition\"", "\"$(pre)s$k\""] :
                         ["\"$(pre)s$k\"", "\"$(pre)t$transition\""],
               Attributes(:label=>"$(δ[k])", :labelfontsize=>"6"))
           for k in collect(keys(δ)) if δ[k] != 0]
end

function edgify(δ::Dict{Tuple{Int64, Bool}, Int64}, transition, reverse::Bool; pre="", lw=3.0)
  return [Edge(reverse ? ["\"$(pre)t$transition\"", "\"$(pre)s$(k[1])\""] :
                         ["\"$(pre)s$(k[1])\"", "\"$(pre)t$transition\""],
               Attributes(:label=>"$(δ[k])", :labelfontsize=>"6",
                          :penwidth=>(k[2] ? "$lw" : "1.0")))
           for k in collect(keys(δ)) if δ[k] != 0]
end

function to_json(model, sens::Dict)

	sensDict = Dict(keys(sens).=>values(sens))
	rateDict = Dict(keys(rates(model)).=>values(rates(model)))
	concDict = Dict(keys(concentrations(model)).=>values(concentrations(model)))
	results = Dict(:model => "SimpleSIRDPetriNetClassic",
	               :rates => rateDict,
	               :init_concentrations => concDict,
	    :tspan => Dict(:start => 0.0, :end => 100.0),
	               :sensDict => sensDict)
	open("SIRD_int_SR.json", "w") do f
	    JSON.print(f, results, 2)
	end
end
end
