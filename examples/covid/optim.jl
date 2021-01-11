using AlgebraicPetri
using AlgebraicPetri.Epidemiology

using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Programs.RelationalPrograms

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

# #### SIR Model:

# define model
sir = PetriNet(3, (1,2)=>(2,2), 2=>3)
Graph(sir)

# define initial states and transition rates, then
# create, solve, and visualize ODE problem

u0 = [10, 1, 0]
β = [0.4, 0.4]

# The C-Set representation has direct support for generating a DiffEq vector field

prob = ODEProblem(vectorfield(sir),u0,(0.0,7.5),β);
sol = solve(prob,Tsit5())

# Create mock data to optimize to
t = collect(range(0,stop=7.5,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(3)) for i in 1:length(t)])
data = convert(Array,randomized)

using DiffEqParamEstim
using Optim
using LeastSquaresOptim

# Create loss function
cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data), maxiters=10000,verbose=false)
# Optimize with initial parameter guesses of 0, 0
result_bfgs = optimize(cost_function, [0.0, 0.0], BFGS())
result_bfgs.minimizer

# Create loss function
cost_function = build_lsoptim_objective(prob,t,data,Tsit5())
# Initial parameter guesses of 0, 0
x = [0.0, 0.0]
# Optimize
res = optimize!(LeastSquaresProblem(x = x, f! = cost_function,
                output_length = length(t)*length(prob.u0)),
                Dogleg(LeastSquaresOptim.LSMR()))
res.minimizer

# Run ODE Solver on optimized parameters
prob = ODEProblem(vectorfield(sir),u0,(0.0,7.5),res.minimizer);
sol = solve(prob,Tsit5())
plot(sol)