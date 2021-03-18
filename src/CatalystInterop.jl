module CatalystInterop
  using AlgebraicPetri
  using Catlab.CategoricalAlgebra
  using ...Catalyst
  import ...Catalyst: ReactionSystem

  counter(a) = [count(==(i),a) for i in unique(a)]

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
