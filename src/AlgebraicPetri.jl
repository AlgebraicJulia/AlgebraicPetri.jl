""" Computing in the category of finite sets and Petri cospans
"""
module AlgebraicPetri

export TheoryPetriNet, PetriNet, OpenPetriNetOb, AbstractPetriNet, ns, nt, ni, no,
  add_species!, add_transition!, add_transitions!,
  add_input!, add_inputs!, add_output!, add_outputs!, inputs, outputs,
  TransitionMatrices, vectorfield,
  TheoryLabelledPetriNet, LabelledPetriNet, AbstractLabelledPetriNet, sname, tname,
  TheoryReactionNet, ReactionNet, AbstractReactionNet, concentration, concentrations, rate, rates,
  TheoryLabelledReactionNet, LabelledReactionNet, AbstractLabelledReactionNet,
  Open, OpenPetriNet, OpenLabelledPetriNet, OpenReactionNet, OpenLabelledReactionNet,
  OpenPetriNetOb, OpenLabelledPetriNetOb, OpenReactionNetOb, OpenLabelledReactionNetOb

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Present
using Catlab.Theories
using Petri
using LabelledArrays
using LinearAlgebra: mul!
import Petri: Model, Graph, vectorfield 

vectorify(n::Vector) = n
vectorify(n::Tuple) = length(n) == 1 ? [n] : n
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

# Petri Nets
############

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

Open(p::AbstractPetriNet) = OpenPetriNet(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open(p::AbstractPetriNet, legs...) = OpenPetriNet(p, map(l->FinFunction(l, ns(p)), legs)...)
Open(n, p::AbstractPetriNet, m) = Open(p, n, m)

# PetriNet([:S, :I, :R], :infection=>((1, 2), 3))

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

ns(p::AbstractPetriNet) = nparts(p,:S)
nt(p::AbstractPetriNet) = nparts(p,:T)
ni(p::AbstractPetriNet) = nparts(p,:I)
no(p::AbstractPetriNet) = nparts(p,:O)

add_species!(p::AbstractPetriNet;kw...) = add_part!(p,:S;kw...)
add_species!(p::AbstractPetriNet,n;kw...) = add_parts!(p,:S,n;kw...)

add_transition!(p::AbstractPetriNet;kw...) = add_part!(p,:T;kw...)
add_transitions!(p::AbstractPetriNet,n;kw...) = add_parts!(p,:T,n;kw...)

add_input!(p::AbstractPetriNet,t,s;kw...) = add_part!(p,:I;it=t,is=s,kw...)
add_inputs!(p::AbstractPetriNet,n,t,s;kw...) = add_parts!(p,:I,n;it=t,is=s,kw...)

add_output!(p::AbstractPetriNet,t,s;kw...) = add_part!(p,:O;ot=t,os=s,kw...)
add_outputs!(p::AbstractPetriNet,n,t,s;kw...) = add_parts!(p,:O,n;ot=t,os=s,kw...)

sname(p::AbstractPetriNet,s) = (1:ns(p))[s]
tname(p::AbstractPetriNet,t) = (1:nt(p))[t]

# Note: although indexing makes this pretty fast, it is often faster to bulk-convert
# the PetriNet net into a transition matrix, if you are working with all of the transitions
inputs(p::AbstractPetriNet,t) = subpart(p,incident(p,t,:it),:is)
outputs(p::AbstractPetriNet,t) = subpart(p,incident(p,t,:ot),:os)

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

# Interoperability with Petri.jl
Petri.Model(p::AbstractPetriNet) = begin
  ts = TransitionMatrices(p)
  t_in = map(i->Dict(k=>v for (k,v) in enumerate(ts.input[i,:]) if v != 0), 1:nt(p))
  t_out = map(i->Dict(k=>v for (k,v) in enumerate(ts.output[i,:]) if v != 0), 1:nt(p))

  Δ = Dict(i=>t for (i,t) in enumerate(zip(t_in, t_out)))
  return Petri.Model(ns(p), Δ)
end

Petri.Model(p::Union{AbstractLabelledPetriNet, AbstractLabelledReactionNet}) = begin
  snames = [sname(p, s) for s in 1:ns(p)]
  tnames = [tname(p, t) for t in 1:nt(p)]
  ts = TransitionMatrices(p)
  t_in = map(i->LVector(;[(snames[k]=>v) for (k,v) in enumerate(ts.input[i,:]) if v != 0]...), 1:nt(p))
  t_out = map(i->LVector(;[(snames[k]=>v) for (k,v) in enumerate(ts.output[i,:]) if v != 0]...), 1:nt(p))

  Δ = LVector(;[(tnames[i]=>t) for (i,t) in enumerate(zip(t_in, t_out))]...)
  return Petri.Model(collect(values(snames)), Δ)
end

concentration(p::AbstractLabelledReactionNet,s) = subpart(p,s,:concentration)
rate(p::AbstractLabelledReactionNet,t) = subpart(p,t,:rate)

concentrations(p::AbstractLabelledReactionNet) = begin
  snames = [sname(p, s) for s in 1:ns(p)]
  LVector(;[(snames[s]=>concentration(p, s)) for s in 1:ns(p)]...)
end

rates(p::AbstractLabelledReactionNet) = begin
  tnames = [tname(p, s) for s in 1:nt(p)]
  LVector(;[(tnames[t]=>rate(p, t)) for t in 1:nt(p)]...)
end

include("visualization.jl")
include("Epidemiology.jl")


end
