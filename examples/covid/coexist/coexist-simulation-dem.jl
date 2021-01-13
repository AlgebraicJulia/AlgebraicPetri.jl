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

add_index(epi_model::LabelledPetriNet, ind::Int) = begin
    new_petri = deepcopy(epi_model)
    snames = subpart(epi_model, :sname)
    tnames = subpart(epi_model, :tname)
    set_subpart!(new_petri, :sname, [Symbol("$(name)_$ind") for name in snames])
    set_subpart!(new_petri, :tname, [Symbol("$(name)_$ind") for name in tnames])
    new_petri
end

function dem_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph, sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol})
    epi_func(ind) = bundle_legs(Open(add_index(epi_model, ind)), [tuple([1:ns(epi_model);]...)])
    conn(x, y) = bundle_legs(Open(
                    dem_connection(epi_func, sus_state::Symbol, exp_state::Symbol,
                                    inf_states::Array{Symbol}, x, y)),
                    [tuple([1:ns(epi_model);]...), tuple([(ns(epi_model)+1):(2*ns(epi_model));]...)])
    stratify(epi_func, connection_graph, conn)
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

F(ex, n) = oapply(ex, Dict(
    :exposure=>exposure_petri(Symbol(:exp_, n), Symbol(:S,n), Symbol(:I,n), Symbol(:E,n)),
    :exposure_e=>exposure_petri(Symbol(:exp_e, n), Symbol(:S,n), Symbol(:E,n),Symbol(:E,n)),
    :exposure_i2=>exposure_petri(Symbol(:exp_i2, n), Symbol(:S,n), Symbol(:I2,n), Symbol(:E,n)),
    :exposure_a=>exposure_petri(Symbol(:exp_a, n), Symbol(:S,n), Symbol(:A,n),Symbol(:E,n)),
    :progression=>spontaneous_petri(Symbol(:prog_, n), Symbol(:I,n), Symbol(:I2,n)),
    :asymptomatic_infection=>spontaneous_petri(Symbol(:asymp_, n), Symbol(:E,n), Symbol(:A,n)),
    :illness=>spontaneous_petri(Symbol(:ill_, n), Symbol(:E,n), Symbol(:I,n)),
    :asymptomatic_recovery=>spontaneous_petri(Symbol(:arec_, n), Symbol(:A,n), Symbol(:R,n)),
    :recovery=>spontaneous_petri(Symbol(:rec_, n), Symbol(:I2,n), Symbol(:R,n)),
    :recover_late=>spontaneous_petri(Symbol(:rec2_, n), Symbol(:R,n), Symbol(:R2,n)),
    :death=>spontaneous_petri(Symbol(:death2_, n), Symbol(:I2,n), Symbol(:D,n))));

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
AlgebraicPetri.Graph(F(coexist, 1))

# Define an `oapply` function that can be used to create a model of
# cross exposure between two sets of populations

F_cx(ex, x,y) = oapply(ex, Dict(
    :exposure=>exposure_petri(Symbol(:exp_, x,y), Symbol(:S,x), Symbol(:I,y), Symbol(:E,x)),
    :exposure_e=>exposure_petri(Symbol(:exp_e, x,y), Symbol(:S,x), Symbol(:E,y),Symbol(:E,x)),
    :exposure_a=>exposure_petri(Symbol(:exp_a, x,y), Symbol(:S,x), Symbol(:A,y),Symbol(:E,x)),
    :exposure_i2=>exposure_petri(Symbol(:exp_i2, x,y), Symbol(:S,x), Symbol(:I2,y), Symbol(:E,x)),
    :exposure′=>exposure_petri(Symbol(:exp_, y,x), Symbol(:S,y), Symbol(:I,x), Symbol(:E,y)),
    :exposure_e′=>exposure_petri(Symbol(:exp_e, y,x), Symbol(:S,y), Symbol(:E,x),Symbol(:E,y)),
    :exposure_a′=>exposure_petri(Symbol(:exp_a, y,x), Symbol(:S,y), Symbol(:A,x),Symbol(:E,y)),
    :exposure_i2′=>exposure_petri(Symbol(:exp_i2, y,x), Symbol(:S,y), Symbol(:I2,x), Symbol(:E,y))
  ),
  Dict(
    :s=>ob(Symbol(:S, x)),
    :e=>ob(Symbol(:E, x)),
    :a=>ob(Symbol(:A, x)),
    :i=>ob(Symbol(:I, x)),
    :i2=>ob(Symbol(:I2, x)),
    :r=>ob(Symbol(:R, x)),
    :r2=>ob(Symbol(:R2, x)),
    :d=>ob(Symbol(:D, x)),
    :s′=>ob(Symbol(:S, y)),
    :e′=>ob(Symbol(:E, y)),
    :a′=>ob(Symbol(:A, y)),
    :i′=>ob(Symbol(:I, y)),
    :i2′=>ob(Symbol(:I2, y)),
    :r′=>ob(Symbol(:R, y)),
    :r2′=>ob(Symbol(:R2, y)),
    :d′=>ob(Symbol(:D, y))
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

bundled_cross(x,y) = bundle_legs(F_cx(crossexposure, x, y), [tuple([1:8;]...), tuple([9:16;]...)])
bundled_coex(x) = bundle_legs(F(coexist, x), [tuple([1:8;]...)])
epi_coex(x) = bundled_coex(x+2)
epi_cross(x,y) = bundled_cross(x+2,y+2)

cross_3 = Catlab.Graphs.BasicGraphs.Graph(3)
add_edges!(cross_3, [1,2,3], [2,3,1])

#-
threeNCoexist_algpetri = apex(dem_strat(apex(epi_coex(2)), cross_3, :S4, :E4, [:I4, :E4, :I24, :A4]))
AlgebraicPetri.Graph(threeNCoexist_algpetri);
save_fig(AlgebraicPetri.Graph(threeNCoexist_algpetri), "3ncoexist_petri", "svg"); # hide

# ![3-generation COEXIST model petri net](3ncoexist_petri.svg)

# We can JSON to convert this Petri net into an
# easily shareable format

JSON.print(threeNCoexist_algpetri.tables)
