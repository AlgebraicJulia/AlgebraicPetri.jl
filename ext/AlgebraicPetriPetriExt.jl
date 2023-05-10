module AlgebraicPetriPetriExt

using AlgebraicPetri
using Catlab.CategoricalAlgebra

# TODO: Remove after dropping support for <Julia 1.9
isdefined(Base, :get_extension) ? (import Petri) : (import ..Petri)

Petri.Model(p::AbstractPetriNet) = begin
  ts = TransitionMatrices(p)

  if has_subpart(p, :sname) && has_subpart(p, :tname)
    snames = [sname(p, s) for s in 1:ns(p)]
    tnames = [tname(p, t) for t in 1:nt(p)]
    t_in = map(i->Dict(snames[k]=>v for (k,v) in enumerate(ts.input[i,:]) if v != 0), 1:nt(p))
    t_out = map(i->Dict(snames[k]=>v for (k,v) in enumerate(ts.output[i,:]) if v != 0), 1:nt(p))
    Δ = Dict(tnames[i]=>t for (i,t) in enumerate(zip(t_in, t_out)))
    S = collect(values(snames)) 
  else
    t_in = map(i->Dict(k=>v for (k,v) in enumerate(ts.input[i,:]) if v != 0), 1:nt(p))
    t_out = map(i->Dict(k=>v for (k,v) in enumerate(ts.output[i,:]) if v != 0), 1:nt(p))
    Δ = Dict(i=>t for (i,t) in enumerate(zip(t_in, t_out)))
    S = ns(p)
  end

  return Petri.Model(ns(p), Δ)
end

end
