module Petri

using Catlab
using Catlab.CategoricalAlgebra.CSets
using Catlab.Present
using Catlab.Theories

# Petri Nets
############

@present TheoryPetri(FreeCategory) begin
    T::Ob
    S::Ob
    I::Ob
    O::Ob

    it::Hom(I,T)
    is::Hom(I,S)
    ot::Hom(O,T)
    os::Hom(O,S)
end

const AbstractPetri = AbstractCSetType(TheoryPetri)
const Petri = CSetType(TheoryPetri, index=[:is,:it,:os,:ot])

ns(p::AbstractPetri) = nparts(p,:S)
nt(p::AbstractPetri) = nparts(p,:T)
ni(p::AbstractPetri) = nparts(p,:I)
no(p::AbstractPetri) = nparts(p,:O)

add_species!(p::AbstractPetri) = add_part!(p,:S)
add_species!(p::AbstractPetri,n) = add_parts!(p,:S,n)

add_transition!(p::AbstractPetri) = add_part!(p,:T)
add_transitions!(p::AbstractPetri,n) = add_parts!(p,:T,n)

add_input!(p::AbstractPetri,t,s) = add_part!(p,:I,(is=s,it=t))
add_inputs!(p::AbstractPetri,n,t,s) = add_parts!(p,:I,n,(is=s,it=t))

add_output!(p::AbstractPetri,t,s) = add_part!(p,:O,(os=s,ot=t))
add_output!(p::AbstractPetri,n,t,s) = add_parts!(p,:O,n,(os=s,ot=t))

# Note: although indexing makes this pretty fast, it is often faster to bulk-convert
# the Petri net into a transition matrix, if you are working with all of the transitions
inputs(p::AbstractPetri,t) = subpart(p,incident(p,t,:it),:is)
outputs(p::AbstractPetri,t) = subpart(p,incident(p,t,:ot),os)

totransmat(p::AbstractPetri) = begin
    A = zeros(Int64,(nt(p),ns(p),2))
    for i in 1:ni(p)
        A[subpart(p,i,:it),subpart(p,i,:is),1] += 1
    end
    for i in 1:no(p)
        A[subpart(p,i,:ot),subpart(p,o,:os),2] += 1
    end
    A
end

fromtransmat(A::Array{Int64,3}) = begin
    (m,n,_) = size(A)
    p = Petri()
    add_species!(p,n)
    add_transitions!(p,m)
    for i in 1:m
        for j in 1:n
            add_inputs!(p,A[i,j,1],i,j)
            add_outputs!(p,A[i,j,2],i,j)
        end
    end
    p
end

rate_eq(p::AbstractPetri) = begin
    n,m = nt(p),ns(p)
    rates = zeros(m)
    A = totransmat(p)
    f(du,u,p,t) = begin
        for i in 1:n
            rates[i] = p[i] * prod(u[j] ^ T[i,j,1] for j in 1:m)
        end
        for j in 1:m
            du[j] = sum(rates[i] * (T[i,j,2] - T[i,j,1]) for i in 1:n)
        end
    end
    f
end

# Reaction Nets
###############

@present TheoryReactionNet <: TheoryPetri begin
    number::Ob

    rate::Hom(T, number)
    u0::Hom(S, number)
end

const AbstractReactionNet = AbstractCSetType(TheoryReactionNet, data=[:number])
const ReactionNet = CSetType(TheoryReactionNet, index=[:is,:it,:os,:ot], data=[:number])

end
