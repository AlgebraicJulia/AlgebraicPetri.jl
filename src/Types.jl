module Types
export TheoryPetriNet, PetriNet, AbstractPetriNet, ns, nt, ni, no,
  add_species!, add_transition!, add_transitions!,
  add_input!, add_inputs!, add_output!, add_outputs!, inputs, outputs,
  TransitionMatrices, rate_eq,
  TheoryLabelledPetriNet, LabelledPetriNet, AbstractLabelledPetriNet, sname, tname,
  TheoryReactionNet, ReactionNet, AbstractReactionNet, concentration, rate,
  TheoryLabelledReactionNet, LabelledReactionNet, AbstractLabelledReactionNet

using Catlab
using Catlab.CategoricalAlgebra.CSets
using Catlab.Present
using Catlab.Theories
using Petri
import Petri: Model

vectorify(n) = begin
  if !(typeof(n) <: Union{Vector,Tuple}) || (typeof(n) <: Tuple && length(n) == 1)
    [n]
  else
    n
  end
end

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

# PetriNet([:S, :I, :R], :infection=>((1, 2), 3))

PetriNet(n,ts...) = begin
  p = PetriNet()
  add_species!(p, n)
  add_transitions!(p, length(ts))
  for (i,(ins,outs)) in enumerate(ts)
    if !(typeof(ins) <: Union{Vector,Tuple}) || length(ins) == 1 ins = [ins] end
    if !(typeof(outs) <: Union{Vector,Tuple}) || length(outs) == 1 outs = [outs] end
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

rate_eq(p::AbstractPetriNet) = begin
  tm = TransitionMatrices(p)
  f(du,u,p,t) = begin
    log_rates = log.(p) .+ tm.input * log.(u)
    du = (tm.output - tm.input) * exp.(log_rates)
  end
  f
end

@present TheoryLabelledPetriNet <: TheoryPetriNet begin
  Name::Data

  tname::Attr(T, Name)
  sname::Attr(S, Name)
end

const AbstractLabelledPetriNet = AbstractACSetType(TheoryLabelledPetriNet)
const LabelledPetriNet = ACSetType(TheoryLabelledPetriNet, index=[:it,:is,:ot,:os]){Symbol}

LabelledPetriNet(n,ts...) = begin
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

sname(p::AbstractLabelledPetriNet,s) = subpart(p,s,:sname)
tname(p::AbstractLabelledPetriNet,t) = subpart(p,t,:tname)

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

ReactionNet{R,C}(n,ts...) where {R,C} = begin
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

@present TheoryLabelledReactionNet <: TheoryReactionNet begin
  Name::Data

  tname::Attr(T, Name)
  sname::Attr(S, Name)
end

const AbstractLabelledReactionNet = AbstractACSetType(TheoryLabelledReactionNet)
const LabelledReactionNet{R,C} = ACSetType(TheoryLabelledReactionNet, index=[:it,:is,:ot,:os]){R,C,Symbol}

# Ex. LabelledReactionNet{Number, Int}((:S=>990, :I=>10, :R=>0), (:inf, .3/1000)=>((:S, :I)=>(:I,:I)), (:rec, .2)=>(:I=>:R))
LabelledReactionNet{R,C}(n,ts...) where {R,C} = begin
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

# Interoperability with Petri.jl
Petri.Model(p::Union{PetriNet,ReactionNet}) = begin
  ts = TransitionMatrices(p)
  return Petri.Model(1:ns(p), collect(zip(eachrow(ts.input), eachrow(ts.output))))
end

Petri.Model(p::Union{LabelledPetriNet,LabelledReactionNet}) = begin
  snames = Dict(s=>sname(p,s) for s in 1:ns(p))
  tnames = Dict(t=>tname(p,t) for t in 1:nt(p))
  ts = TransitionMatrices(p)
  t_in = map(i->Dict(snames[k]=>v for (k,v) in enumerate(i) if v != 0), eachrow(ts.input))
  t_out = map(i->Dict(snames[k]=>v for (k,v) in enumerate(i) if v != 0), eachrow(ts.output))
  Δ = Dict(tnames[k]=>v for (k,v) in enumerate(zip(t_in, t_out)))
  return Petri.Model(collect(values(snames)), Δ)
end

end