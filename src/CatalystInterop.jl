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
  function ReactionSystem(pn::AbstractPetriNet;
                          fixed_rates=Dict{Symbol, Union{Nothing, Number}}(),
                          var_states=Array{Symbol,1}())

    # Define the variables (stored flat as [var_rates..., var_concs...])
    var_rates = filter(t -> !(t ∈ keys(fixed_rates)), tnames(pn))
    r_params = length(var_rates)
    c_params = length(var_states)
    num_params = r_params+c_params
    @parameters t k[1:num_params]
    @variables S[collect(1:ns(pn))](t)

    cur_param = 0
    rxs = map(1:nt(pn)) do i
      cur_t = tnames(pn)[i]
      inpts = pn[incident(pn, i, :it),:is]
      otpts = pn[incident(pn, i, :ot),:os]
      in_count = collect(counter(inpts))
      ot_count = collect(counter(otpts))
      Reaction(cur_t ∈ var_rates ? (cur_param += 1; k[cur_param]) : fixed_rates[cur_t],
               S[unique(inpts)], S[unique(otpts)], in_count, ot_count)
    end

    ReactionSystem(rxs, t, S, k)
  end
end
