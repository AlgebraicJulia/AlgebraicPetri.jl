using Requires

function __init__()
  @require Petri="4259d249-1051-49fa-8328-3f8ab9391c33" begin
    import .Petri

    # Interoperability with Petri.jl
    Petri.Model(p::AbstractPetriNet) = begin
      ts = TransitionMatrices(p)
      t_in = map(i->Dict(k=>v for (k,v) in enumerate(ts.input[i,:]) if v != 0), 1:nt(p))
      t_out = map(i->Dict(k=>v for (k,v) in enumerate(ts.output[i,:]) if v != 0), 1:nt(p))

      Δ = Dict(i=>t for (i,t) in enumerate(zip(t_in, t_out)))
      return Petri.Model(ns(p), Δ)
    end

    Petri.Model(p::Union{AbstractLabelledPetriNet, AbstractLabelledReactionNet}) = begin
      snames = [sname(p, s) for s in 1:ns(p)]
      tnames = [tname(p, t) for t in 1:nt(p)]
      ts = TransitionMatrices(p)
      t_in = map(i->LVector(;[(snames[k]=>v) for (k,v) in enumerate(ts.input[i,:]) if v != 0]...), 1:nt(p))
      t_out = map(i->LVector(;[(snames[k]=>v) for (k,v) in enumerate(ts.output[i,:]) if v != 0]...), 1:nt(p))

      Δ = LVector(;[(tnames[i]=>t) for (i,t) in enumerate(zip(t_in, t_out))]...)
      return Petri.Model(collect(values(snames)), Δ)
    end
  end
end
