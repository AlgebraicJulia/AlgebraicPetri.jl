""" Types to compose Petri nets along transitions.
"""

module OpenTransitions

using AlgebraicPetri: 
    LabelledPetriNetUntyped, LabelledPetriNet, PetriNet, AbstractPetriNet, 
    ns, nt, tname
using Catlab
using Catlab.CategoricalAlgebra

export OpenLabelledPetriNetObT, OpenLabelledPetriNetT, OpenT,
    OpenPetriNetObT, OpenPetriNetT

const OpenPetriNetObT, OpenPetriNetT = OpenCSetTypes(PetriNet,:T)

""" OpenT(p::AbstractPetriNet)

Converts a PetriNet to an OpenPetriNetT where each transition is exposed as a leg of
the cospan. The OpenPetriNetT can be composed over an undirected wiring diagram.
"""
OpenT(p::AbstractPetriNet) = OpenPetriNetT(p, map(x->FinFunction([x], nt(p)), 1:nt(p))...)

""" OpenT(p::AbstractPetriNet, legs...)

Generates on OpenPetriNetT with legs bundled as described by `legs`
"""
OpenT(p::AbstractPetriNet, legs...) = OpenPetriNetT(p, map(l->FinFunction(l, nt(p)), legs)...)

""" OpenT(n, p::AbstractPetriNet, m)

Generates on OpenPetriNetT with two legs, `n` and `m`
"""
OpenT(n, p::AbstractPetriNet, m) = OpenT(p, n, m)

const OpenLabelledPetriNetObUntypedT, OpenLabelledPetriNetUntypedT = OpenACSetTypes(LabelledPetriNetUntyped,:T)
const OpenLabelledPetriNetObT, OpenLabelledPetriNetT = OpenLabelledPetriNetObUntypedT{Symbol}, OpenLabelledPetriNetUntypedT{Symbol}

OpenT(p::LabelledPetriNet) = OpenLabelledPetriNetT(p, map(x->FinFunction([x], nt(p)), 1:nt(p))...)
OpenT(p::LabelledPetriNet, legs...) = begin
    t_idx = Dict(tname(p, t)=>t for t in 1:nt(p))
    OpenLabelledPetriNetT(p, map(l->FinFunction(map(i->t_idx[i], l), nt(p)), legs)...)
end
OpenT(n, p::LabelledPetriNet, m) = OpenT(p, n, m)

end
