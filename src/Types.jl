module Types
export TheoryPetri, Petri, AbstractPetri, ns, nt, ni, no,
  add_species!, add_transition!, add_transitions!,
  add_input!, add_inputs!, add_output!, add_outputs!, inputs, outputs,
  TransitionMatrices, rate_eq,
  TheoryReactionNet, ReactionNet, AbstractReactionNet

using Catlab
using Catlab.CategoricalAlgebra.CSets
using Catlab.Present
using Catlab.Theories
using Petri

# Petri Nets
############

@present TheoryPetri(FreeSchema) begin
  T::Ob
  S::Ob
  I::Ob
  O::Ob

  it::Hom(I,T)
  is::Hom(I,S)
  ot::Hom(O,T)
  os::Hom(O,S)
end

const AbstractPetri = AbstractACSetType(TheoryPetri)
const Petri = CSetType(TheoryPetri,index=[:it,:is,:ot,:os])

(pt::typeof(Petri))(n,ts...) = begin
  p = pt()
  add_species!(p, n)
  add_transitions!(p, length(ts))
  for (i,(ins,outs)) in enumerate(ts)
    if length(ins) == 1 ins = [ins] end
    if length(outs) == 1 outs = [outs] end
    add_inputs!(p, length(ins), repeat([i], length(ins)), collect(ins))
    add_outputs!(p, length(outs), repeat([i], length(outs)), collect(outs))
  end
  p
end

ns(p::AbstractPetri) = nparts(p,:S)
nt(p::AbstractPetri) = nparts(p,:T)
ni(p::AbstractPetri) = nparts(p,:I)
no(p::AbstractPetri) = nparts(p,:O)

add_species!(p::AbstractPetri;kw...) = add_part!(p,:S;kw...)
add_species!(p::AbstractPetri,n;kw...) = add_parts!(p,:S,n;kw...)

add_transition!(p::AbstractPetri;kw...) = add_part!(p,:T;kw...)
add_transitions!(p::AbstractPetri,n;kw...) = add_parts!(p,:T,n;kw...)

add_input!(p::AbstractPetri,t,s;kw...) = add_part!(p,:I;it=t,is=s,kw...)
add_inputs!(p::AbstractPetri,n,t,s;kw...) = add_parts!(p,:I,n;it=t,is=s,kw...)

add_output!(p::AbstractPetri,t,s;kw...) = add_part!(p,:O;ot=t,os=s,kw...)
add_outputs!(p::AbstractPetri,n,t,s;kw...) = add_parts!(p,:O,n;ot=t,os=s,kw...)

# Note: although indexing makes this pretty fast, it is often faster to bulk-convert
# the Petri net into a transition matrix, if you are working with all of the transitions
inputs(p::AbstractPetri,t) = subpart(p,incident(p,t,:it),:is)
outputs(p::AbstractPetri,t) = subpart(p,incident(p,t,:ot),:os)

struct TransitionMatrices
  input::Matrix{Int}
  output::Matrix{Int}
  TransitionMatrices(p::AbstractPetri) = begin
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

Petri(tm::TransitionMatrices) = begin
  (m,n) = size(tm.input)
  p = Petri()
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

rate_eq(p::AbstractPetri) = begin
  tm = TransitionMatrices(p)
  f(du,u,p,t) = begin
    log_rates = log.(p) .+ tm.input * log.(u)
    du = (tm.output - tm.input) * exp.(log_rates)
  end
  f
end

# Reaction Nets
###############

@present TheoryReactionNet <: TheoryPetri begin
  Rate::Data
  Concentration::Data

  rate::Attr(T, Rate)
  concentration::Attr(S, Concentration)
end

const AbstractReactionNet = AbstractACSetType(TheoryReactionNet)
const ReactionNet = ACSetType(TheoryReactionNet, index=[:it,:is,:ot,:os])

end
