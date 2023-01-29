""" Types to compose Petri nets along transitions.
"""

module OpenTransitions

using AlgebraicPetri
using Catlab
using Catlab.CategoricalAlgebra

export OpenLabelledPetriNetObT, OpenLabelledPetriNetT, OpenT,
    OpenPetriNetObT, OpenPetriNetT

const OpenPetriNetObT, OpenPetriNetT = OpenCSetTypes(PetriNet,:T)

OpenT(p::AbstractPetriNet) = OpenPetriNetT(p, map(x->FinFunction([x], nt(p)), 1:nt(p))...)
OpenT(p::AbstractPetriNet, legs...) = OpenPetriNetT(p, map(l->FinFunction(l, nt(p)), legs)...)
OpenT(n, p::AbstractPetriNet, m) = OpenT(p, n, m)

# const OpenPetriNetOb, OpenPetriNet = OpenCSetTypes(PetriNet,:S)

# """ Open(p::AbstractPetriNet)

# Converts a PetriNet to an OpenPetriNet where each state is exposed as a leg of
# the cospan. The OpenPetriNet can be composed over an undirected wiring diagram
# (see this
# [blog post](https://www.algebraicjulia.org/blog/post/2020/11/structured-cospans-2/)
# for a description of this compositional tooling)
# """
# Open(p::AbstractPetriNet) = OpenPetriNet(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)

# """ Open(p::AbstractPetriNet, legs...)

# Generates on OpenPetriNet with legs bundled as described by `legs`
# """
# Open(p::AbstractPetriNet, legs...) = OpenPetriNet(p, map(l->FinFunction(l, ns(p)), legs)...)

# """ Open(n, p::AbstractPetriNet, m)

# Generates on OpenPetriNet with two legs, `n` and `m`
# """
# Open(n, p::AbstractPetriNet, m) = Open(p, n, m)


const OpenLabelledPetriNetObUntypedT, OpenLabelledPetriNetUntypedT = OpenACSetTypes(LabelledPetriNetUntyped,:T)
const OpenLabelledPetriNetObT, OpenLabelledPetriNetT = OpenLabelledPetriNetObUntypedT{Symbol}, OpenLabelledPetriNetUntypedT{Symbol}

OpenT(p::AbstractLabelledPetriNet) = OpenLabelledPetriNetT(p, map(x->FinFunction([x], nt(p)), 1:nt(p))...)
OpenT(p::AbstractLabelledPetriNet, legs...) = begin
    t_idx = Dict(tname(p, t)=>t for t in 1:nt(p))
    OpenLabelledPetriNetT(p, map(l->FinFunction(map(i->t_idx[i], l), nt(p)), legs)...)
end
OpenT(n, p::AbstractLabelledPetriNet, m) = OpenT(p, n, m)

# const OpenLabelledPetriNetObUntyped, OpenLabelledPetriNetUntyped = OpenACSetTypes(LabelledPetriNetUntyped,:S)
# const OpenLabelledPetriNetOb, OpenLabelledPetriNet = OpenLabelledPetriNetObUntyped{Symbol}, OpenLabelledPetriNetUntyped{Symbol}


# Open(p::AbstractLabelledPetriNet) = OpenLabelledPetriNet(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
# Open(p::AbstractLabelledPetriNet, legs...) = begin
#   s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
#   OpenLabelledPetriNet(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)),legs)...)
# end
# Open(n, p::AbstractLabelledPetriNet, m) = Open(p, n, m)


end