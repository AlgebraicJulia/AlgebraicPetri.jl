""" Supports conversion from PetriNets to ModelingToolkit ODESystems

This provides access to the ModelingToolkit library.
"""

module ModelingToolkitInterop
  using AlgebraicPetri.BilayerNetworks
  using Catlab.CategoricalAlgebra: parts, incident
  using ...ModelingToolkit
  import ...ModelingToolkit: ODESystem

  """ Convert a general PetriNet to an ODESystem

  This conversion forgets any labels or rates provided, and converts all
  parameters and variables into symbols. It does preserve the ordering of
  transitions and states though (Transition 1 has a rate of k[1], state 1 has a
  concentration of S[1])
  """
  function ModelingToolkit.ODESystem(bn::Union{AbstractLabelledBilayerNetwork,AbstractBilayerNetwork}; name = :PetriNet)
    t = (@variables t)[1]
    D = Differential(t)
    symbolic_vars = map(bn[:variable]) do v
        (@variables $v(t))[1]
    end
    symbolic_params = map(bn[:parameter]) do p
        (@parameters $p)[1]
    end

    ϕs = map(parts(bn, :Box)) do b
      p = symbolic_params[b]
      vars = mapreduce(*, incident(bn, b, :call), init = p) do i
        j = bn[i, :arg]
        return symbolic_vars[j]
      end
    end

    infs = map(parts(bn, :Qout)) do tv
      flux = mapreduce(+, incident(bn, tv, :infusion), init = 0) do wa
        j = bn[wa, :influx]
        return ϕs[j]
      end
      flux -= mapreduce(+, incident(bn, tv, :effusion), init = 0) do wa
        j = bn[wa, :efflux]
        return ϕs[j]
      end
    end

    # We assume bn[:tanvar] ⊆ bn[:variable] here
    tanvar_idxs = indexin(bn[:tanvar], bn[:variable])
    zparts = zip(tanvar_idxs, infs)

    eqs = Equation[D(symbolic_vars[j::Int]) ~ rhs for (j, rhs) in zparts]
    ODESystem(eqs, t, symbolic_vars, symbolic_params, name=name)
  end
end
