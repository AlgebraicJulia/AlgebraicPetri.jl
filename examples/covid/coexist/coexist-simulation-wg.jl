# # [COEXIST Multi-Generational COVID Model](@id coexist_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/coexist/coexist.ipynb)

using AlgebraicPetri

using OrdinaryDiffEq
using Plots
using JSON

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using Catlab.Theories
using Catlab.WiringDiagrams
using Catlab.Graphics.Graphviz: run_graphviz

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>"1"));
save_fig(g, fname::AbstractString, format::AbstractString) = begin
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, g, format=format)
    end
end

# Define helper functions for defining the two types of
# reactions in an epidemiology Model. Either a state
# spontaneously changes, or one state causes another to change

using Catlab.Graphs.BasicGraphs
function cross_uwd(g::Catlab.Graphs.BasicGraphs.Graph)
    rel = RelationDiagram{Symbol}(0)
    
    # Add populations
    juncs = add_junctions!(rel, nv(g), variable=[Symbol("pop$i") for i in 1:nv(g)])
    srcs = subpart(g, :src)
    tgts = subpart(g, :tgt)
    # Add cross boxes
    add_parts!(rel, :Box, ne(g), name=[Symbol("cross_$(srcs[i])_$(tgts[i])") for i in 1:ne(g)])
    add_parts!(rel, :Port, ne(g), junction=srcs, box=1:ne(g))
    add_parts!(rel, :Port, ne(g), junction=tgts, box=1:ne(g))
    
    # Add epidemiology model boxes
    boxes = add_parts!(rel, :Box, nv(g), name=[Symbol("ep$i") for i in 1:nv(g)])
    add_parts!(rel, :Port, nv(g), junction=juncs, box=boxes)
    rel 
end

function stratify(epi_petri::Function, connection_graph::Catlab.Graphs.BasicGraphs.Graph, diffusion_petri::Function)
    conn_uwd = cross_uwd(connection_graph)
    
    # Calls diffusion_petri for each edge as (src, tgt)
    ep_map = Dict([Symbol("ep$i")=>epi_petri(i) for i in 1:nv(connection_graph)])
    srcs = subpart(connection_graph, :src)
    tgts = subpart(connection_graph, :tgt)
    for i in 1:ne(connection_graph)
        ep_map[Symbol("cross_$(srcs[i])_$(tgts[i])")] = diffusion_petri(srcs[i], tgts[i])
    end
    
    oapply(conn_uwd, ep_map)
end

dem_connection(epi_model::Function, sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol}, x::Int, y::Int) = begin
    append_ind(x::Symbol, ind::Int) = Symbol("$(x)_$ind")
    
    ep1 = apex(epi_model(x))
    ep2 = apex(epi_model(y))
    
    sus1 = append_ind(sus_state, x)
    sus2 = append_ind(sus_state, y)
    exp1 = append_ind(exp_state, x)
    exp2 = append_ind(exp_state, y)
    inf1 = [append_ind(inf, x) for inf in inf_states]
    inf2 = [append_ind(inf, y) for inf in inf_states]
    
    LabelledPetriNet(vcat(subpart(ep1, :sname), subpart(ep2, :sname)),
                     vcat([Symbol("crx_$(sus1)_$(inf)")=>((sus1, inf)=>(inf, exp1)) for inf in inf2],
                          [Symbol("crx_$(sus2)_$(inf)")=>((sus2, inf)=>(inf, exp2)) for inf in inf1])...)
    
end

diff_connection(epi_model::Function, x::Int, y::Int) = begin
    ep1 = apex(epi_model(x))
    ep2 = apex(epi_model(y))
    states1 = subpart(ep1, :sname)
    states2 = subpart(ep2, :sname)
    LabelledPetriNet(vcat(states1, states2),
                     vcat([Symbol("diff_$(states1[i])_$(states2[i])")=>(states1[i]=>states2[i]) for i in 1:ns(ep1)],
                          [Symbol("diff_$(states2[i])_$(states1[i])")=>(states2[i]=>states1[i]) for i in 1:ns(ep2)])...)
end

add_index(epi_model::LabelledPetriNet, ind::Int) = begin
    new_petri = deepcopy(epi_model)
    snames = subpart(epi_model, :sname)
    tnames = subpart(epi_model, :tname)
    set_subpart!(new_petri, :sname, [Symbol("$(name)_$ind") for name in snames])
    set_subpart!(new_petri, :tname, [Symbol("$(name)_$ind") for name in tnames])
    new_petri
end


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
display_uwd(coexist)
AlgebraicPetri.Graph(F(coexist))

# Define an `oapply` function that can be used to create a model of
# cross exposure between two sets of populations

F_cx(ex) = oapply(ex, Dict(
    :exposure=>exposure_petri(Symbol(:exp_), Symbol(:S), Symbol(:I), Symbol(:E)),
    :exposure_e=>exposure_petri(Symbol(:exp_e), Symbol(:S), Symbol(:E),Symbol(:E)),
    :exposure_a=>exposure_petri(Symbol(:exp_a), Symbol(:S), Symbol(:A),Symbol(:E)),
    :exposure_i2=>exposure_petri(Symbol(:exp_i2), Symbol(:S), Symbol(:I2), Symbol(:E)),
    :exposure′=>exposure_petri(Symbol(:exp_), Symbol(:S), Symbol(:I), Symbol(:E)),
    :exposure_e′=>exposure_petri(Symbol(:exp_e), Symbol(:S), Symbol(:E),Symbol(:E)),
    :exposure_a′=>exposure_petri(Symbol(:exp_a), Symbol(:S), Symbol(:A),Symbol(:E)),
    :exposure_i2′=>exposure_petri(Symbol(:exp_i2), Symbol(:S), Symbol(:I2), Symbol(:E))
  ),
  Dict(
    :s=>ob(Symbol(:S)),
    :e=>ob(Symbol(:E)),
    :a=>ob(Symbol(:A)),
    :i=>ob(Symbol(:I)),
   :i2=>ob(Symbol(:I2)),
    :r=>ob(Symbol(:R)),
    :r2=>ob(Symbol(:R2)),
    :d=>ob(Symbol(:D)),
    :s′=>ob(Symbol(:S)),
    :e′=>ob(Symbol(:E)),
    :a′=>ob(Symbol(:A)),
    :i′=>ob(Symbol(:I)),
    :i2′=>ob(Symbol(:I2)),
    :r′=>ob(Symbol(:R)),
    :r2′=>ob(Symbol(:R2)),
    :d′=>ob(Symbol(:D))
  ));

# Use this new presentation to define a model
# of cross exposure between two populations

crossexposure = @relation (s, e, i, i2, a, r, r2, d, s′, e′, i′, i2′, a′, r′, r2′, d′) begin
    exposure(s, i′, e)
    exposure_i2(s, i2′, e)
    exposure_a(s, a′, e)
    exposure_e(s, e′, e)
    exposure′(s′, i, e′)
    exposure_i2′(s′, i2, e′)
    exposure_a′(s′, a, e′)
    exposure_e′(s′, e, e′)
end;
display_uwd(crossexposure)

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
    epi_func(ind) = bundle_legs(Open(add_index(epi_model, ind)), [tuple([1:ns(epi_model);]...)])
    conn(x, y) = bundle_legs(Open(diff_connection(epi_func, x, y)), [tuple([1:ns(epi_model);]...), tuple([(ns(epi_model)+1):(2*ns(epi_model));]...)])
    stratify(epi_func, connection_graph, conn)
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
    epi_func(ind) = bundle_legs(Open(add_index(epi_model, ind)), [tuple([1:ns(epi_model);]...)])
    conn(x, y) = bundle_legs(Open(
                    dem_connection(epi_func, sus_state::Symbol, exp_state::Symbol,
                                    inf_states::Array{Symbol}, x, y)),
                    [tuple([1:ns(epi_model);]...), tuple([(ns(epi_model)+1):(2*ns(epi_model));]...)])
    stratify(epi_func, connection_graph, conn)
end

cycle_3 = Catlab.Graphs.BasicGraphs.Graph(3)
add_edges!(cycle_3, [1,2,3], [2,3,1])

coex = apex(F(coexist));

coex_dem = apex(dem_strat(coex, cycle_3, :S, :E, [:I, :E, :I2, :A]));
coex_diff = apex(diff_strat(coex, cycle_3))

save_fig(AlgebraicPetri.Graph(coex_dem), "coex_crx", "svg");
save_fig(AlgebraicPetri.Graph(coex_diff), "coex_diff", "svg");
