""" Useful tools for analysing and comparing between multiple PEtriNet models
"""

module ModelComparison

using AlgebraicPetri
using Catlab
import Catlab.CategoricalAlgebra: Subobject

export petri_homomorphisms, compare, PetriSubobject

""" petri_homomorphisms
Calculates all homomorphisms between two PetriNet models and returns an array
of ACSetTransformations that describe these homomorphisms.

NOTE:
This function does restrict to homomorphisms which preserve the transition
signatures (number of input/output wires).
"""
function petri_homomorphisms(p1::AbstractPetriNet, p2::AbstractPetriNet; kw...)
  cat = ACSetCategory(ACSetCat(p1))
  results = ACSetTransformation[]
  sigsTS = Set{Vector{Int}}()
  homs = homomorphisms(PetriNet(p1), PetriNet(p2); cat, kw...)
  for transform in homs
    if all([length(inputs(p1, t)) == length(inputs(p2, transform.components[:T](t))) &&
            length(outputs(p1, t)) == length(outputs(p2, transform.components[:T](t)))
            for t in 1:nt(p1)]) &&
        vcat(collect(transform.components[:T]), collect(transform.components[:S])) ∉ sigsTS

      push!(results, ACSetTransformation(deepcopy(transform.components), p1, p2))
      push!(sigsTS, vcat(collect(transform.components[:T]),
                         collect(transform.components[:S])))
    end
  end
  results
end

""" compare
Calculates all homomorphisms from a single Petri net (`apex`) to other Petri
nets. This is returned as a Multispan with `apex` as the apex and the
homomorphisms as the legs.
"""
function compare(apex::T, pn::T...; monic=[:S, :T, :I, :O], kw...) where {T <: AbstractPetriNet}
  legs = map(enumerate(pn)) do (i, p)
    hp = petri_homomorphisms(apex, p; monic=monic, kw...)
    isempty(hp) && error("No homomorphism from apex to leg $i which is monic in $monic")
    hp
  end
  Multispan(apex, vcat(legs...))
end

end
