# +
include("epidemiology.jl")

using OrdinaryDiffEq
using Catlab.Graphics.Graphviz: run_graphviz

# +
save_wd(ex, fname::AbstractString) = begin
    g = display_wd(ex)
    open(fname, "w") do io
        run_graphviz(io, g, format="svg")
    end
end
save_graph(g, fname::AbstractString) = begin
    open(fname, "w") do io
        run_graphviz(io, g, format="svg")
    end
end

splot(soln, fname) = begin
    p = plot(soln, linewidth=5, legend=false)
    savefig(p, fname)
end

splot(soln, t, fname) = begin
    p = plot(t, soln, linewidth=5, legend=false)
    savefig(p, fname)
end

splotchannels(sol, dir) = begin
    mkpath(dir)
    jp(p) = joinpath(dir, p)
    splot(sol, jp("solution.svg"))
    pltdims = [3,8,13]
    splot(sol[pltdims,:]', sol.t, jp("infecteds.svg"))
    splot(sum(sol(sol.t, idxs=pltdims), dims=1)', sol.t, jp("infecteds_sum.svg"))
    pltdims = [3,8,13] .+ 1
    splot(sol[pltdims,:]', sol.t, jp( "recovered.svg" ))
    splot(sum(sol(sol.t, idxs=pltdims), dims=1)', sol.t, jp( "recovered_sum.svg" ))
end
# -

save_wd(sir, "sir.svg")
save_graph(Graph(p_sir), "sir_graph.svg")

# define initial states and transition rates
u0 = [10.0, 1, 0]
p = [0.4, 0.4]
# create and solve ODE problem
prob = ODEProblem(p_sir,u0,(0.0,7.5),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
splot(sol, "sir_soln.svg")

save_wd(seir, "seir.svg")
save_graph(Graph(p_seir), "seir_graph.svg")

# define initial states and transition rates
u0 = [10.0, 1, 0, 0]
p = [0.9, 0.2, 0.5]
# create and solve ODE problem
prob = ODEProblem(p_seir,u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
splot(sol, "seir_soln.svg")

save_wd(seird, "seird.svg")
save_graph(Graph(p_seird), "seird_graph.svg")

# define initial states and transition rates
u0 = [10.0, 1, 0, 0, 0]
p = [0.9, 0.2, 0.5, 0.1]
# create and solve ODE problem
prob = ODEProblem(p_seird,u0,(0.0,15.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
# visualize the solution
splot(sol, "seird_soln.svg")

include("covid.jl")

[save_graph(Graph(decoration(p)), pname) for (p,pname) in [(spontaneous_petri, "spont.svg"),
                                               (transmission_petri, "transmission.svg"),
                                               (exposure_petri, "exposure.svg"),
                                               (travel_petri, "travel.svg")]]

save_wd(seird_city, "seird_city.svg")
save_graph(Graph(decoration(F(seird_city))), "seird_graph.svg")

save_wd(seird_3, "seird_3.svg")
save_graph(Graph(p_seird_3), "seird_3_graph.svg")

pltdims = [4,9, 14]
# Define time frame, 3 months
tspan = (0.0,90.0)
# Define initial states
u0 = zeros(Float64, length(base(pc_seird_3)))
u0[1]  = 10000
u0[6]  = 10000
u0[11] = 10000
u0[2]  = 1
# Define transition rates
params = seirdparams(3, 5)
# Generate and solve resulting ODE
prob = ODEProblem(p_seird_3,u0,tspan,params)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
splotchannels(sol, "seird")

function drawsol(f, a, b, p, dir)
    tspan = (0, 240.0)
    println("INFO: Processing $dir")
    dynparams = waveparams(f, a/sum(u0), b/sum(u0),p)
    mkpath(dir)
    prob = ODEProblem(p_seird_3,u0,tspan, dynparams)
    sol = OrdinaryDiffEq.solve(prob,Tsit5(), saveat=1:1:tspan[2])
    splotchannels(sol, dir)
    savefig(plot(sol.t, f(a,b,p), ylim=(0,12), legend=false, linewidth=5), "$dir/transmissability.svg")
end

a,b = 10, 1
dir = "dynamic/mitigation/"
drawsol(asymptotic,  a,b,1/6, "$dir/asymptotic")
drawsol(triangleasm, a,b,25 , "$dir/trianglesm")
drawsol(coswave,     a,b,25 , "$dir/coswave")
drawsol(sincwave,    a,b,100, "$dir/sincwave")
drawsol(modsincwave, a,b,100, "$dir/modsincwave")

a,b = 10, 1/2
dir = "dynamic/containment"
drawsol(asymptotic,  a,b,1/6, "$dir/asymptotic")
drawsol(triangleasm, a,b, 25, "$dir/trianglesm")
drawsol(coswave,     a,b, 25, "$dir/coswave")
drawsol(sincwave,    a,b,100, "$dir/sincwave")
drawsol(modsincwave, a,b,100, "$dir/modsincwave")
