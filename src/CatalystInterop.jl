""" Supports conversion from PetriNets to Catalyst ReactionSystems

This provides access to the parameter estimation, optimization, and sensitivity
tooling provided in the Catalyst library
"""

module CatalystInterop
  using AlgebraicPetri
  using Catlab.CategoricalAlgebra
  using ...Catalyst
  import ...Catalyst: ReactionSystem

  counter(a) = [count(==(i),a) for i in unique(a)]

  """ Convert a general PetriNet to a ReactionSystem

  This conversion forgets any labels or rates provided, and converts all
  parameters and variables into symbols. It does preserve the ordering of
  transitions and states though (Transition 1 has a rate of k[1], state 1 has a
  concentration of S[1])
  """
  function ReactionSystem(pn::AbstractPetriNet)
    @parameters t k[1:nt(pn)]
    @variables S[collect(1:ns(pn))](t)

    rxs = map(1:nt(pn)) do t
      inpts = pn[incident(pn, t, :it),:is]
      otpts = pn[incident(pn, t, :ot),:os]
      in_count = collect(counter(inpts))
      ot_count = collect(counter(otpts))
      Reaction(k[t], S[unique(inpts)], S[unique(otpts)], in_count, ot_count)
    end

    ReactionSystem(rxs, t, S, k)
  end
end
