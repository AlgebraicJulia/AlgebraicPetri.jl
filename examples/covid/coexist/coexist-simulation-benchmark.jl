# # [COEXIST Multi-Generational COVID Model](@id coexist_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/coexist/coexist.ipynb)
using Dates
using Profile
println(now(), " Starting compilation")
flush(stdout)
include("stratification.jl")

println(now(), " Starting coexist construction")
flush(stdout)


ob(x::Symbol) = codom(Open([x], LabelledPetriNet(x), [x])).ob;
function spontaneous_petri(transition::Symbol, s::Symbol, t::Symbol)
    Open(LabelledPetriNet(unique([s, t]), transition=>(s=>t)))
end;
function exposure_petri(transition::Symbol, s::Symbol, e::Symbol, t::Symbol)
    Open(LabelledPetriNet(unique([s,e,t]), transition=>((s,e)=>(t,e))))
end;

# Define an `oapply` function that connects the building block Petri nets to
# the operations we will use in the model.

F(ex) = oapply(ex, Dict(
    :exposure=>exposure_petri(Symbol(:exp_), Symbol(:S), Symbol(:I), Symbol(:E)),
    :exposure_e=>exposure_petri(Symbol(:exp_e), Symbol(:S), Symbol(:E),Symbol(:E)),
    :exposure_i2=>exposure_petri(Symbol(:exp_i2), Symbol(:S), Symbol(:I2), Symbol(:E)),
    :exposure_a=>exposure_petri(Symbol(:exp_a), Symbol(:S), Symbol(:A),Symbol(:E)),
    :progression=>spontaneous_petri(Symbol(:prog_), Symbol(:I), Symbol(:I2)),
    :asymptomatic_infection=>spontaneous_petri(Symbol(:asymp_), Symbol(:E), Symbol(:A)),
    :illness=>spontaneous_petri(Symbol(:ill_), Symbol(:E), Symbol(:I)),
    :asymptomatic_recovery=>spontaneous_petri(Symbol(:arec_), Symbol(:A), Symbol(:R)),
    :recovery=>spontaneous_petri(Symbol(:rec_), Symbol(:I2), Symbol(:R)),
    :recover_late=>spontaneous_petri(Symbol(:rec2_), Symbol(:R), Symbol(:R2)),
    :death=>spontaneous_petri(Symbol(:death2_), Symbol(:I2), Symbol(:D))));

# Define the COEXIST model using the `@relation` macro

coexist = @relation (s, e, i, i2, a, r, r2, d) begin
    exposure(s, i, e)
    exposure_i2(s, i2, e)
    exposure_a(s, a, e)
    exposure_e(s, e, e)
    asymptomatic_infection(e, a)
    asymptomatic_recovery(a, r)
    illness(e, i)
    progression(i, i2)
    death(i2, d)
    recovery(i2, r)
    recover_late(r, r2)
end;
#display_uwd(coexist)
#AlgebraicPetri.Graph(F(coexist))


# To combine these two models, we need to create a final relational model and
# use the `bundle_legs` function in our `oapply` that enables us to model 3
# population wires instead of each individual state as a wire. Each of these
# populations has their own COEXIST model, and interact through cross exposure


# ![3-generation COEXIST model petri net](3ncoexist_petri.svg)

# We can JSON to convert this Petri net into an
# easily shareable format


""" diff_strat(epi_model, connection_graph)

  This function takes in a LabelledPetriNet and a graph which describes
  geographical connections. It returns a LabelledPetriNet which models
  diffusion between geographic populations described by the given graph.
"""

function diff_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph)
  diff_conn = diff_petri(epi_model)
  stratify(epi_model, (diff_conn, connection_graph))
end

""" dem_strat(epi_model, connection_graph, sus_state, exp_state, inf_states)

  This function takes in a LabelledPetriNet and a graph which describes
  infection connections between populations. It also takes in the symbol used
  for susceptible states, the symbol used for the exposed state, and an array
  of symbols for states that can expose susceptible states. It returns a
  LabelledPetriNet which models diffusion between populations described by the
  given graph.
"""

function dem_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph, sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol})
  dem_conn = dem_petri(epi_model, sus_state, exp_state, inf_states)
  stratify(epi_model, (dem_conn, connection_graph))
end

coex = apex(F(coexist));
println(now(), " Finished generating coexist")
flush(stdout)

benchmark(demo, cities) = 
begin
  println(now(), " Generating dem strat with $(demo)-clique")
  flush(stdout)
  coex_dem = dem_strat(coex, clique(demo), :S, :E, [:I, :E, :I2, :A]);
  println(now(), " Finished generating dem with $(demo)-clique")
  flush(stdout)
  println(now(), " Generating city strat with $(cities)-cycle")
  flush(stdout)
  coex_diff = diff_strat(coex_dem, cycle(cities))
  println(now(), " Finished city strat with $(cities)-cycle")
  flush(stdout)
end

println("\nBenchmark with 1 city with 1 population")
benchmark(1,1)
println("\nBenchmark with 10 cities with 10 populations")
benchmark(10,10)
println("\nBenchmark with 100 cities with 100 populations")
benchmark(100,100)

dem_conn = dem_petri(coex, :S, :E, [:I, :E, :I2, :A])
coex_dem = stratify(coex, (dem_conn, cycle(3)));
diff_conn = diff_petri(coex_dem)
coex_diff = stratify(coex_dem, (diff_conn, cycle(3)))


save_fig(AlgebraicPetri.Graph(coex_dem), "3cycle_dem", "svg");
save_fig(AlgebraicPetri.Graph(coex_diff), "3cycle_diff", "svg");

dem_conn = dem_petri(coex, :S, :E, [:I, :E, :I2, :A])
diff_conn = diff_petri(coex)
coex_diff = stratify(coex, (diff_conn, cycle(5)), (dem_conn, clique(5)))
save_fig(AlgebraicPetri.Graph(coex_diff), "dem_diff", "svg");
#
#Profile.clear()
#@profile benchmark(100, 100)
#
#Profile.print()
