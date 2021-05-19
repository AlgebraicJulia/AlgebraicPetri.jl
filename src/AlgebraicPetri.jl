""" Computing in the category of finite sets and Petri cospans
"""
module AlgebraicPetri

export TheoryPetriNet, PetriNet, OpenPetriNetOb, AbstractPetriNet, ns, nt, ni, no,
  add_species!, add_transition!, add_transitions!,
  add_input!, add_inputs!, add_output!, add_outputs!, inputs, outputs,
  TransitionMatrices, vectorfield,
  TheoryLabelledPetriNet, LabelledPetriNet, AbstractLabelledPetriNet, sname, tname, snames, tnames,
  TheoryReactionNet, ReactionNet, AbstractReactionNet, concentration, concentrations, rate, rates,
  TheoryLabelledReactionNet, LabelledReactionNet, AbstractLabelledReactionNet,
  Open, OpenPetriNet, OpenLabelledPetriNet, OpenReactionNet, OpenLabelledReactionNet,
  OpenPetriNetOb, OpenLabelledPetriNetOb, OpenReactionNetOb, OpenLabelledReactionNetOb

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Present
using Catlab.Theories
using LabelledArrays
using LinearAlgebra: mul!

vectorify(n::Vector) = n
vectorify(n::Tuple) = length(n) == 1 ? [n] : n
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

# Petri Nets
############

""" ACSet definition for a Petri net.

See Catlab.jl documentation for description of the @present syntax.
"""
@present TheoryPetriNet(FreeSchema) begin
  T::Ob
  S::Ob
  I::Ob
  O::Ob

  it::Hom(I,T)
  is::Hom(I,S)
  ot::Hom(O,T)
  os::Hom(O,S)
end

const AbstractPetriNet = AbstractACSetType(TheoryPetriNet)
const PetriNet = CSetType(TheoryPetriNet,index=[:it,:is,:ot,:os])
const OpenPetriNetOb, OpenPetriNet = OpenCSetTypes(PetriNet,:S)

""" Open(p::AbstractPetriNet)

Converts a PetriNet to an OpenPetriNet where each state is exposed as a leg of
the cospan. The OpenPetriNet can be composed over an undirected wiring diagram
(see this
[blog post](https://www.algebraicjulia.org/blog/post/2020/11/structured-cospans-2/)
for a description of this compositional tooling)
"""
Open(p::AbstractPetriNet) = OpenPetriNet(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)

""" Open(p::AbstractPetriNet, legs...)

Generates on OpenPetriNet with legs bundled as described by `legs`
"""
Open(p::AbstractPetriNet, legs...) = OpenPetriNet(p, map(l->FinFunction(l, ns(p)), legs)...)

""" Open(n, p::AbstractPetriNet, m)

Generates on OpenPetriNet with two legs, `n` and `m`
"""
Open(n, p::AbstractPetriNet, m) = Open(p, n, m)

""" PetriNet(n::Int, ts::Vararg{Union{Pair,Tuple}})

Constructs a PetriNet object with `n` states and transitions described by `ts`.
Transitions are given as `(input_states)=>(output_states)`.

A PetriNet modelling the SIR model with 3 states and 2 transitions can be
constructed as follows:
```@example
PetriNet(3, (1,2)=>(2,2), 2=>3)
```
"""
PetriNet(n::Int, ts::Vararg{Union{Pair,Tuple}}) = begin
  p = PetriNet()
  add_species!(p, n)
  add_transitions!(p, length(ts))
  for (i,(ins,outs)) in enumerate(ts)
    ins = vectorify(ins)
    outs = vectorify(outs)
    add_inputs!(p, length(ins), repeat([i], length(ins)), collect(ins))
    add_outputs!(p, length(outs), repeat([i], length(outs)), collect(outs))
  end
  p
end

""" Number of states in a Petri net
"""
ns(p::AbstractPetriNet) = nparts(p,:S)

""" Number of transitions in a Petri net
"""
nt(p::AbstractPetriNet) = nparts(p,:T)

""" Number of input relationships in a Petri net
"""
ni(p::AbstractPetriNet) = nparts(p,:I)

""" Number of output relationships in a Petri net
"""
no(p::AbstractPetriNet) = nparts(p,:O)

""" Add a species to the Petri net. Label and concentration can be provided
depending on the kind of Petri net.

Returns the ID of the species
"""
add_species!(p::AbstractPetriNet;kw...) = add_part!(p,:S;kw...)

""" Add `n` species to the Petri net. Label and concentration can be provided
depending on the kind of Petri net.

Returns the ID of the species
"""
add_species!(p::AbstractPetriNet,n;kw...) = add_parts!(p,:S,n;kw...)

""" Add a transition to the Petri net. Label and rate can be provided
depending on the kind of Petri net.

Returns the ID of the transition
"""
add_transition!(p::AbstractPetriNet;kw...) = add_part!(p,:T;kw...)

""" Add `n` transitions to the Petri net. Label and rate can be provided
depending on the kind of Petri net.

Returns the ID of the transition
"""
add_transitions!(p::AbstractPetriNet,n;kw...) = add_parts!(p,:T,n;kw...)

""" add_input!(p::AbstractPetriNet,t,s;kw...)

Add an input relationship to the Petri net between the transition `t` and species `s`.

Returns the ID of the input relationship
"""
add_input!(p::AbstractPetriNet,t,s;kw...) = add_part!(p,:I;it=t,is=s,kw...)

""" add_inputs!(p::AbstractPetriNet,n,t,s;kw...)

Add input relationships to the Petri net between the transitions `t` and species `s`.

Returns the ID of the input relationship
"""
add_inputs!(p::AbstractPetriNet,n,t,s;kw...) = add_parts!(p,:I,n;it=t,is=s,kw...)

""" add_output!(p::AbstractPetriNet,t,s;kw...)

Add an output relationship to the Petri net between the transition `t` and species `s`.

Returns the ID of the input relationship
"""
add_output!(p::AbstractPetriNet,t,s;kw...) = add_part!(p,:O;ot=t,os=s,kw...)

""" add_outputs!(p::AbstractPetriNet,n,t,s;kw...)

Add output relationships to the Petri net between the transitions `t` and species `s`.

Returns the ID of the input relationship
"""
add_outputs!(p::AbstractPetriNet,n,t,s;kw...) = add_parts!(p,:O,n;ot=t,os=s,kw...)

""" Name of species

Note that this returns an index if labels are not present in the PetriNet
"""
sname(p::AbstractPetriNet,s) = (1:ns(p))[s]

""" Name of transition

Note that this returns an index if labels are not present in the PetriNet
"""
tname(p::AbstractPetriNet,t) = (1:nt(p))[t]

""" Names of species in  a Petri net

Note that this returns indices if labels are not present in the PetriNet
"""
snames(p::AbstractPetriNet) = map(s->sname(p, s), 1:ns(p))

""" Names of transitions in  a Petri net

Note that this returns indices if labels are not present in the PetriNet
"""
tnames(p::AbstractPetriNet) = map(t->tname(p, t), 1:nt(p))

# Note: although indexing makes this pretty fast, it is often faster to bulk-convert
# the PetriNet net into a transition matrix, if you are working with all of the transitions
""" Input relationships for a transition
"""
inputs(p::AbstractPetriNet,t) = subpart(p,incident(p,t,:it),:is)

""" Output relationships for a transition
"""
outputs(p::AbstractPetriNet,t) = subpart(p,incident(p,t,:ot),:os)

""" TransitionMatrices

This data structure stores the transition matrix of an AbstractPetriNet object.
This is primarily used for constructing the vectorfield representation of the
Petri net.
"""
struct TransitionMatrices
  input::Matrix{Int}
  output::Matrix{Int}
  TransitionMatrices(p::AbstractPetriNet) = begin
    input, output = zeros(Int,(nt(p),ns(p))), zeros(Int,(nt(p),ns(p)))
    for i in 1:ni(p)
      input[subpart(p,i,:it),subpart(p,i,:is)] += 1
    end
    for o in 1:no(p)
      output[subpart(p,o,:ot),subpart(p,o,:os)] += 1
    end
    new(input,output)
  end
end

""" PetriNet(tm::TransitionMatrices)

Constructs a PetriNet from its TransitionMatrices representation.
"""
PetriNet(tm::TransitionMatrices) = begin
  (m,n) = size(tm.input)
  p = PetriNet()
  add_species!(p,n)
  add_transitions!(p,m)
  for i in 1:m
    for j in 1:n
      add_inputs!(p,tm.input[i,j],i,j)
      add_outputs!(p,tm.output[i,j],i,j)
    end
  end
  p
end

valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end

""" vectorfield(pn::AbstractPetriNet)

Generates a Julia function which calculates the vectorfield of the Petri net
being simulated under the law of mass action.

The resulting function has a signature of the form `f!(du, u, p, t)` and can be
passed to the DifferentialEquations.jl solver package.
"""
vectorfield(pn::AbstractPetriNet) = begin
  tm = TransitionMatrices(pn)
  dt = tm.output - tm.input
  f(du,u,p,t) = begin
    rates = zeros(eltype(du),nt(pn))
    u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
    p_m = [p[tname(pn, i)] for i in 1:nt(pn)]
    for i in 1:nt(pn)
      rates[i] = valueat(p_m[i],u,t) * prod(u_m[j] ^ tm.input[i,j] for j in 1:ns(pn))
    end
    for j in 1:ns(pn)
      du[sname(pn, j)] = sum(rates[i] * dt[i,j] for i in 1:nt(pn))
    end
    return du
  end
  return f
end

""" ACSet definition for a Petri net with labels on transitions and states.

See Catlab.jl documentation for description of the @present syntax.
"""
@present TheoryLabelledPetriNet <: TheoryPetriNet begin
  Name::Data

  tname::Attr(T, Name)
  sname::Attr(S, Name)
end

const AbstractLabelledPetriNet = AbstractACSetType(TheoryLabelledPetriNet)
const LabelledPetriNetUntyped = ACSetType(TheoryLabelledPetriNet, index=[:it,:is,:ot,:os])
const LabelledPetriNet = LabelledPetriNetUntyped{Symbol}
const OpenLabelledPetriNetObUntyped, OpenLabelledPetriNetUntyped = OpenACSetTypes(LabelledPetriNetUntyped,:S)
const OpenLabelledPetriNetOb, OpenLabelledPetriNet = OpenLabelledPetriNetObUntyped{Symbol}, OpenLabelledPetriNetUntyped{Symbol}


Open(p::AbstractLabelledPetriNet) = OpenLabelledPetriNet(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open(p::AbstractLabelledPetriNet, legs...) = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledPetriNet(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)),legs)...)
end
Open(n, p::AbstractLabelledPetriNet, m) = Open(p, n, m)

""" LabelledPetriNet(n, ts::Vararg{Union{Pair,Tuple}})

Constructs a LabelledPetriNet object with state names as elements of `n` and
labelled transitions described by `ts`.
Transitions are given as `transition_name=>((input_states)=>(output_states))`.

A LabelledPetriNet modelling the SIR model with 3 states and 2 transitions can be
constructed as follows:
```@example
LabelledPetriNet([:S, :I, :R], :inf=>((:S,:I)=>(:I,:I)), :rec=>(:I=>:R))
```
"""
LabelledPetriNet(n, ts::Vararg{Union{Pair,Tuple}}) = begin
  p = LabelledPetriNet()
  n = vectorify(n)
  state_idx = state_dict(n)
  add_species!(p, length(n), sname=n)
  for (name,(ins,outs)) in ts
    i = add_transition!(p, tname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    add_inputs!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    add_outputs!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
  end
  p
end

# Reaction Nets
###############

""" ACSet definition for a Petri net with rates on transitions and
concentrations on states.

See Catlab.jl documentation for description of the @present syntax.
"""
@present TheoryReactionNet <: TheoryPetriNet begin
  Rate::Data
  Concentration::Data

  rate::Attr(T, Rate)
  concentration::Attr(S, Concentration)
end

const AbstractReactionNet = AbstractACSetType(TheoryReactionNet)
const ReactionNet = ACSetType(TheoryReactionNet, index=[:it,:is,:ot,:os])
const OpenReactionNetOb, OpenReactionNet = OpenACSetTypes(ReactionNet,:S)

Open(p::AbstractReactionNet{R,C}) where {R,C} = OpenReactionNet{R,C}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open(p::AbstractReactionNet{R,C}, legs...) where {R,C} = OpenReactionNet{R,C}(p, map(l->FinFunction(l, ns(p)), legs)...)
Open(n, p::AbstractReactionNet, m) = Open(p, n, m)

""" ReactionNet{R,C}(n, ts::Vararg{Union{Pair,Tuple}}) where {R,C}

Constructs a ReactionNet object with state concentrations as elements of `n`
and transitions described by `ts`. `R` is the data type used to store rates and
`C` is the data type used to store concentrations.

Transitions are given as `transition_rate=>((input_states)=>(output_states))`.

A ReactionNet modelling the SIR model with 3 states and 2 transitions, an
initial population of 10 susceptible, 1 infected, 0 recovered and an infection
rate of 0.5 and recovery rate of 0.1 can be
constructed as follows:
```@example
ReactionNet{Float64, Float64}([10,1,0], 0.5=>((1,2)=>(2,2)), 0.1=>(2=>3))
```
"""
ReactionNet{R,C}(n, ts::Vararg{Union{Pair,Tuple}}) where {R,C} = begin
  p = ReactionNet{R,C}()
  add_species!(p, length(n), concentration=n)
  for (i, (rate,(ins,outs))) in enumerate(ts)
    i = add_transition!(p, rate=rate)
    ins = vectorify(ins)
    outs = vectorify(outs)
    add_inputs!(p, length(ins), repeat([i], length(ins)), collect(ins))
    add_outputs!(p, length(outs), repeat([i], length(outs)), collect(outs))
  end
  p
end

concentration(p::AbstractReactionNet,s) = subpart(p,s,:concentration)
rate(p::AbstractReactionNet,t) = subpart(p,t,:rate)

concentrations(p::AbstractReactionNet) = map(s->concentration(p, s), 1:ns(p))
rates(p::AbstractReactionNet) = map(t->rate(p, t), 1:nt(p))

""" ACSet definition for a ReactionNet with labels on transitions and states.

See Catlab.jl documentation for description of the @present syntax.
"""
@present TheoryLabelledReactionNet <: TheoryReactionNet begin
  Name::Data

  tname::Attr(T, Name)
  sname::Attr(S, Name)
end

const AbstractLabelledReactionNet = AbstractACSetType(TheoryLabelledReactionNet)
const LabelledReactionNetUntyped = ACSetType(TheoryLabelledReactionNet, index=[:it,:is,:ot,:os])
const LabelledReactionNet{R,C} = LabelledReactionNetUntyped{R,C,Symbol}
const OpenLabelledReactionNetObUntyped, OpenLabelledReactionNetUntyped = OpenACSetTypes(LabelledReactionNetUntyped,:S)
const OpenLabelledReactionNetOb{R,C} = OpenLabelledReactionNetObUntyped{R,C,Symbol}
const OpenLabelledReactionNet{R,C} = OpenLabelledReactionNetUntyped{R,C,Symbol}

Open(p::AbstractLabelledReactionNet{R,C}) where {R,C} = OpenLabelledReactionNet{R,C}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open(p::AbstractLabelledReactionNet{R,C}, legs...) where {R,C} = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledReactionNet{R,C}(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)), legs)...)
end
Open(n, p::AbstractLabelledReactionNet, m) = Open(p, n, m)

# Ex. LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .3/1000)=>((:S, :I)=>(:I,:I)), (:rec, .2)=>(:I=>:R))
""" LabelledReactionNet{R,C}(n, ts::Vararg{Union{Pair,Tuple}}) where {R,C}

Constructs a LabelledReactionNet object with labelled state concentrations as elements of `n`
and labelled transitions described by `ts`. `R` is the data type used to store rates and
`C` is the data type used to store concentrations.

Transitions are given as `(t_name=>t_rate)=>((input_states)=>(output_states))`.

A LabelledReactionNet modelling the SIR model with 3 states and 2 transitions, an
initial population of 10 susceptible, 1 infected, 0 recovered and an infection
rate of 0.5 and recovery rate of 0.1 can be
constructed as follows:
```@example
ReactionNet{Float64, Float64}([:S=>10,:I=>1,:R=>0], (:inf=>0.5)=>((1,2)=>(2,2)), (:rec=>0.1)=>(2=>3))
```
"""
LabelledReactionNet{R,C}(n, ts::Vararg{Union{Pair,Tuple}}) where {R,C} = begin
  p = LabelledReactionNet{R,C}()
  n = vectorify(n)
  states = map(first, collect(n))
  concentrations = map(last, collect(n))
  state_idx = state_dict(states)
  add_species!(p, length(states), concentration=concentrations, sname=states)
  for (i, ((name,rate),(ins,outs))) in enumerate(ts)
    i = add_transition!(p,rate=rate, tname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    add_inputs!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    add_outputs!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
  end
  p
end

sname(p::Union{AbstractLabelledPetriNet, AbstractLabelledReactionNet},s) = subpart(p,s,:sname)
tname(p::Union{AbstractLabelledPetriNet, AbstractLabelledReactionNet},t) = subpart(p,t,:tname)

# Interoperability between different types

PetriNet(pn::AbstractPetriNet) = begin
  pn′ = PetriNet()
  copy_parts!(pn′, pn)
  pn′
end

LabelledPetriNet(pn::Union{AbstractLabelledPetriNet, AbstractLabelledReactionNet}) = begin
  pn′ = LabelledPetriNet()
  copy_parts!(pn′, pn)
  pn′
end

LabelledPetriNet(pn::AbstractPetriNet, snames, tnames) = begin
  pn′ = LabelledPetriNet()
  copy_parts!(pn′, pn)
  map(k->set_subpart!(pn′, k, :sname, snames[k]), keys(snames))
  map(k->set_subpart!(pn′, k, :tname, tnames[k]), keys(tnames))
  pn′
end

ReactionNet{R,C}(pn::Union{AbstractReactionNet, AbstractLabelledReactionNet}) where {R, C} = begin
  pn′ = ReactionNet{R,C}()
  copy_parts!(pn′, pn)
  pn′
end

ReactionNet{R,C}(pn::AbstractPetriNet, concentrations, rates) where {R, C} = begin
  pn′ = ReactionNet{R,C}()
  copy_parts!(pn′, pn)
  map(k->set_subpart!(pn′, k, :concentration, concentrations[k]), keys(concentrations))
  map(k->set_subpart!(pn′, k, :rate, rates[k]), keys(rates))
  pn′
end

LabelledReactionNet{R,C}(pn::AbstractLabelledReactionNet) where {R, C} = begin
  pn′ = LabelledReactionNet{R,C}()
  copy_parts!(pn′, pn)
  pn′
end

LabelledReactionNet{R,C}(pn::AbstractPetriNet, s_labels, t_labels, concentrations, rates) where {R, C} = begin
  pn′ = LabelledReactionNet{R,C}()
  copy_parts!(pn′, pn)
  map(k->set_subpart!(pn′, k, :sname, s_labels[k]), keys(s_labels))
  map(k->set_subpart!(pn′, k, :tname, t_labels[k]), keys(t_labels))
  map(k->set_subpart!(pn′, k, :concentration, concentrations[k]), keys(concentrations))
  map(k->set_subpart!(pn′, k, :rate, rates[k]), keys(rates))
  pn′
end

LabelledReactionNet{R,C}(pn::Union{AbstractPetriNet}, states, transitions) where {R, C} = begin
  pn′ = LabelledReactionNet{R,C}()
  copy_parts!(pn′, pn)
  for (i, (k, v)) in enumerate(states)
    set_subpart!(pn′, i, :sname, k)
    set_subpart!(pn′, i, :concentration, v)
  end
  for (i, (k, v)) in enumerate(transitions)
    set_subpart!(pn′, i, :tname, k)
    set_subpart!(pn′, i, :rate, v)
  end
  pn′
end

""" Concentration of a ReactionNet
"""
concentration(p::AbstractLabelledReactionNet,s) = subpart(p,s,:concentration)

""" Rate of a RectionNet
"""
rate(p::AbstractLabelledReactionNet,t) = subpart(p,t,:rate)

""" All concentrations of a ReactionNet
"""
concentrations(p::AbstractLabelledReactionNet) = begin
  snames = [sname(p, s) for s in 1:ns(p)]
  LVector(;[(snames[s]=>concentration(p, s)) for s in 1:ns(p)]...)
end

""" All rates of a ReactionNet
"""
rates(p::AbstractLabelledReactionNet) = begin
  tnames = [tname(p, s) for s in 1:nt(p)]
  LVector(;[(tnames[t]=>rate(p, t)) for t in 1:nt(p)]...)
end

include("interoperability.jl")
include("visualization.jl")
include("Epidemiology.jl")


end
