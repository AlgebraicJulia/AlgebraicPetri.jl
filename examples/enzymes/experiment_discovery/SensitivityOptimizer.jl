include("EnzymeReactions.jl")
using .EnzymeReactions: enzyme_generators, enzyme_uwd
using GlobalSensitivity
using QuasiMonteCarlo
using Statistics
using OrdinaryDiffEq
using JSON: parsefile, print
using AlgebraicPetri
using Catlab.Programs
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catalyst: reactions, Sym, ReactionSystem, params
using Zygote
# using Catlab.CategoricalAlgebra, Catlab.Graphics
using Plots
using GalacticOptim, Optim
using ForwardDiff
using LinearAlgebra
using CSV
using DataFrames

function formatStrArr(arr)
    #Change symbol array to string
    B = Array{String}(undef, length(arr))
    for i in 1:length(arr)
         B[i] = string(arr[i])
    end
    return B
end


function generatorFunction()
    gen = enzyme_generators([:K,:S,:L],[:G,:E,:P])
    lfunctor(x) = oapply(x, gen);
    Graph(gen[:catKsubG])
    # Generate the two models we fit to
    uwd_kgp = enzyme_uwd([:S], [:G, :P])
    model_kgp = uwd_kgp |> lfunctor |> apex;
    uwd_kg = enzyme_uwd([:S], [:G, :P])
    model_kg = uwd_kgp |> lfunctor |> apex;
    to_graphviz(uwd_kgp, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"))
    # Convert Petri net to ReactionSystem
    enzyme_kg = ReactionSystem(model_kg)


    """making GSA parameters here"""
    init_rates = parsefile("pred_param_KSSpike.json")

    # rateArr = collect(values(init_rates))
    #

    """Load in the data"""
    G_LUM_CONST = 7.5/64207/39e3/150*1e12/7.1 #pM
    home_dir = "Path"
    fileDir = string(home_dir,"\\K_GS_10.json")
    filedata = parsefile(fileDir)
    samples = filedata["data"]
    t_steps = filedata["time_data"]
    concs = filedata["init_conc"]
    weights = filedata["weights"]

    u0 = ["$c" ∈ keys(concs) ? concs["$c"] : 0.0 for c in snames(model_kgp)]
    # println(u0)
    prob = ODEProblem(enzyme_kg, zeros(ns(model_kgp)), (0.0,1.0), zeros(length(params(enzyme_kg))))
    # println(prob.u0)
    sym_rates = filter(r->(reactions(enzyme_kg)[r].rate isa Sym), 1:nt(model_kgp))
    rate_var = tnames(model_kgp)[sym_rates]
    rate_ind = 1:length(rate_var)
    rate_var_str = formatStrArr(rate_var)
    rate_arr = getindex.(Ref(init_rates),rate_var_str)
    p_init_arr = [log10(init_rates["$r"]) for r in formatStrArr(rate_var)]
    # println(p_init_arr)
    #Isolate the data needed
    times = t_steps[1]["G_deg"]
    G_data = samples[1]["G_deg"]

    loc_G_deg = findall(x->x==:G_deg, snames(model_kgp))
    println(loc_G_deg)
    sol = solve(prob;saveat = 1.0)
    # println(sol.u[100][8])
    # init_sol = solve(prob,u0 = u0, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    #     Tsit5(); saveat = 1.0)

    function runSimulation(parameters)
        # init_sol = solve(prob,u0 = u0, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
        #     Tsit5(); saveat = 1.0)
        # u = [parameters[1], 0.0, 0.0, 0.0, 0.0, parameters[2], 0.0, 0.0, parameters[3], 0.0, 0.0]
        u = [parameters[1], 0.0, 0.0, 0.0, 0.0, 1.28e6, 0.0, 0.0, 6.67e4, 0.0, 0.0]
        # u = [parameters[1], 0.0, 0.0, 0.0, 0.0, 1.28e6, 0.0, 0.0, 0.0, 0.0, 0.0]
        # u = [parameters[1], 0.0, 0.0, 0.0, 0.0, parameters[2], 0.0, 0.0, 6.67e4, 0.0, 0.0, 0.0, 0.0, 0.0, 1.28e6, 0.0, 0.0, 6.67e4, 0.0, 0.0, 0.0, 0.0]
        prob1 = remake(prob;u0 = u)
        sol = solve(prob1,p= 10.0 .^ (p_init_arr[rate_ind]),Tsit5();saveat = 1.0)
        # error = 0
        #Making cost function
        # ival = init_sol.u
        val = sol.u
        return val[end][loc_G_deg]
    # return error

    end
    return runSimulation
end



function createCostFunc(f)
    """cost function to run teh GSA in. GSA_cf is the cost function for the GSA
    It will take in a specific value of the parameters and then perform a GSA on
    those by creating an upper and lower bound of those variables
     """

     function costF(x)
         jcbn = ForwardDiff.jacobian(f,x)
         return -LinearAlgebra.norm(jcbn)
     end

end

test = generatorFunction()

test2 = createCostFunc(test)

u0 = [3.0e4, 0.0, 0.0, 0.0, 0.0, 1.28e6, 0.0, 0.0, 6.67e4, 0.0, 0.0]

ux = [3.0e4]
lb = [2.1e1]
ub = [3.0e4]

KRange = LinRange(0.0,3.0e4,500)
# t = map(test,GRange)
t2 = reduce(vcat,test.(KRange))
plot(KRange,t2,title = "Starting S vs end Gelatin")
savefig("Kresponse.png")
t3 = zeros(length(t2)-1)
for i = 1:length(t2)-1
    d = t2[i+1]-t2[i]
    t3[i] = d
end
plot(KRange[1:length(KRange)-1], t3,ylabel ="Change in degraded Gelatin (nM)", xlabel = "Cathepsin S Starting Concentration (nM)",title = "Change in Degraded Gelatin concentration")
savefig("SresponseDf.png")

a = ForwardDiff.jacobian(test,ux)
b = test2(ux)

optTest = optimize(test2,lb,ub,ux)

println(Optim.minimizer(optTest))
optVal = Optim.minimizer(optTest)[1]

#Now I should plot the results to check
uplot = [optVal, 0.0, 0.0, 0.0, 0.0, 1.28e6, 0.0, 0.0, 6.67e4, 0.0, 0.0]

gen = enzyme_generators([:K,:S,:L],[:G,:E,:P])
lfunctor(x) = oapply(x, gen);
Graph(gen[:catKsubG])
# Generate the two models we fit to
uwd_kgp = enzyme_uwd([ :S], [:G, :P])
model_kgp = uwd_kgp |> lfunctor |> apex;
uwd_kg = enzyme_uwd([ :S], [:G, :P])
model_kg = uwd_kgp |> lfunctor |> apex;
to_graphviz(uwd_kgp, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"))
# Convert Petri net to ReactionSystem
enzyme_kg = ReactionSystem(model_kg)


"""making GSA parameters here"""
init_rates = parsefile("pred_param_KSSpike.json")

"""Load in the data"""
G_LUM_CONST = 7.5/64207/39e3/150*1e12/7.1 #pM
home_dir = "Directory"
fileDir = string(home_dir,"K_GS_10.json")
filedata = parsefile(fileDir)
samples = filedata["data"]
t_steps = filedata["time_data"]
concs = filedata["init_conc"]
weights = filedata["weights"]

uOpt = ["$c" ∈ keys(concs) ? concs["$c"] : 0.0 for c in snames(model_kgp)]
uOpt[1] = optVal
# println(u0)
Optval1 = .25*optVal
Optval2 = .5*optVal
Optval3 = .75*optVal
Optval4 = 1.25*optVal
Optval5 = 1.5*optVal
Optval6 = 1.75*optVal
probOpt = ODEProblem(enzyme_kg, zeros(ns(model_kgp)), (0.0,100.0), zeros(length(params(enzyme_kg))))

sym_rates = filter(r->(reactions(enzyme_kg)[r].rate isa Sym), 1:nt(model_kgp))
rate_var = tnames(model_kgp)[sym_rates]
rate_ind = 1:length(rate_var)
rate_var_str = formatStrArr(rate_var)
rate_arr = getindex.(Ref(init_rates),rate_var_str)
p_init_arr = [log10(init_rates["$r"]) for r in formatStrArr(rate_var)]

solOpt = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)

plot(solOpt,ylims = [0,7.5e4], title = "Simulation with optimal Cathepsin concentration",ylabel ="Concentration (nM)", xlabel = "Time (seconds)", labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot.png")
names = snames(model_kgp)
# prepend!(names, [:timestamp])
names = formatStrArr(names)

dfOpt = DataFrame(solOpt)
# rename!(dfOpt,names)
CSV.write("optimal_conc.csv",dfOpt)

uOpt[1] = Optval1
solOpt1 = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)
dfOpt1 = DataFrame(solOpt1)
# rename!(dfOpt1,names)
CSV.write("optimal_conc25.csv",dfOpt1)

uOpt[1] = Optval2
solOpt2 = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)
dfOpt2 = DataFrame(solOpt2)
# rename!(dfOpt2,names)
CSV.write("optimal_conc5.csv",dfOpt2)

uOpt[1] = Optval3
solOpt3 = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)
dfOpt3 = DataFrame(solOpt3)
# rename!(dfOpt3,names)
CSV.write("optimal_conc75.csv",dfOpt3)


uOpt[1] = Optval4
solOpt4 = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)
dfOpt4 = DataFrame(solOpt4)
# rename!(dfOpt4,names)
CSV.write("optimal_conc125.csv",dfOpt4)


uOpt[1] = Optval5
solOpt5 = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)
dfOpt5 = DataFrame(solOpt5)
# rename!(dfOpt5,names)
CSV.write("optimal_conc15.csv",dfOpt5)

uOpt[1] = Optval6
solOpt6 = solve(probOpt,u0 = uOpt, tspan=(0.0 ,100.0),p=  10.0 .^ (p_init_arr[rate_ind]),
    Tsit5(); saveat = 1.0)
dfOpt6 = DataFrame(solOpt6)
# rename!(dfOpt6,names)
CSV.write("optimal_conc175.csv",dfOpt6)

plot(solOpt1, ylims = [0,7.5e4], ylabel ="Concentration (nM)", xlabel = "Time (seconds)", title = "Cathepsin concentration .25 optimal", label = labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot1.png")

plot(solOpt2,ylims = [0,7.5e4], ylabel ="Concentration (nM)", xlabel = "Time (seconds)", title = "Cathepsin concentration .5 optimal",label = labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot2.png")
plot(solOpt3,ylims = [0,7.5e4], ylabel ="Concentration (nM)", xlabel = "Time (seconds)", title = "Cathepsin concentration .75 optimal",label = labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot3.png")
plot(solOpt4,ylims = [0,7.5e4],ylabel ="Concentration (nM)", xlabel = "Time (seconds)",  title = "Cathepsin concentration 1.25 optimal",label = labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot4.png")
plot(solOpt5,ylims = [0,7.5e4],ylabel ="Concentration (nM)", xlabel = "Time (seconds)", title = "Cathepsin concentration 1.5 optimal",label = labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot5.png")

plot(solOpt5,ylims = [0,7.5e4], ylabel ="Concentration (nM)", xlabel = "Time (seconds)",  title = "Cathepsin concentration 1.75 optimal",label = labels=permutedims(String.(snames(model_kgp))))
savefig("optimalPlot6.png")
