using AlgebraicPetri
using Graphviz_jll
using LabelledArrays

sir = LabelledPetriNet([:S,:I,:R],
    :inf => ((:S,:I) => (:I,:I)),
    :rec => (:I => :R)
)

u0 = LVector(S=10.0, I=1.0, R=0.0);
p = LVector(inf=0.4, rec=0.4);

Graph(sir)

macro test_macro(pn)
    num_t = :(nt($(esc(pn))))
    num_s = :(ns($(esc(pn))))
    # println("our pn has $(esc(num_t)) transitions and $(esc(num_s)) places")
    quote
        ($(num_t), $(num_s))
        # println("our pn has $(num_t) transitions and $(num_s) places")
    end 
end

@test_macro(sir)

@macroexpand @test_macro(sir)

using Catlab
using Catlab.ACSetInterface
using Catlab.DenseACSets

acset_schema(sir) isa StructACSet
typeof(acset_schema(sir))

objects(acset_schema(sir))


expr = Expr(:block)
push!(expr.args, :(1 + 2))
eval(expr)


# vf
using Catlab.CategoricalAlgebra

tm = TransitionMatrices(sir)


struct TransitionMatrices
    input::Matrix{Int}
    output::Matrix{Int}
    TransitionMatrices(p::AbstractPetriNet) = begin
      input, output = zeros(Int, (nt(p), ns(p))), zeros(Int, (nt(p), ns(p)))
      for i in 1:ni(p)
        input[subpart(p, i, :it), subpart(p, i, :is)] += 1
      end
      for o in 1:no(p)
        output[subpart(p, o, :ot), subpart(p, o, :os)] += 1
      end
      new(input, output)
    end
  end

vectorfield(pn::AbstractPetriNet) = begin
    tm = TransitionMatrices(pn)
    dt = tm.output - tm.input
    f(du, u, p, t) = begin
      rates = zeros(eltype(du), nt(pn))
      u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
      p_m = [p[tname(pn, i)] for i in 1:nt(pn)]
      for i in 1:nt(pn)
        rates[i] = valueat(p_m[i], u, t) * prod(u_m[j]^tm.input[i, j] for j in 1:ns(pn))
      end
      for j in 1:ns(pn)
        du[sname(pn, j)] = sum(rates[i] * dt[i, j] for i in 1:nt(pn); init=0.0)
      end
      return du
    end
    return f
end

for i in 1:nt(sir)
    # inputs to this transition
    i_ix = incident(sir, i, :it)    
    rate = valueat(p_m[i], u, t) * prod(u[sname(sir, j)] for j in subpart(sir, i_ix, :is))
    # outputs from this transition
    o_ix = incident(sir, i, :ot)
    du[sname(sir, o_ix)] .+= rate
end

vectorfield2(pn::AbstractPetriNet) = begin
    f(du, u, p, t) = begin
        p_m = [p[tname(pn, i)] for i in 1:nt(pn)]
        for i in 1:nt(pn)
            # inputs to this transition
            it_ix = incident(pn, i, :it)    
            is_ix = subpart(pn, it_ix, :is) # input places

            rate = AlgebraicPetri.valueat(p_m[i], u, t) * prod(u[sname(pn, j)] for j in is_ix)

            # outputs from this transition
            ot_ix = incident(pn, i, :ot)
            os_ix = subpart(pn, ot_ix, :os) # output place

            # # net change
            # @inbounds du[sname(pn, os_ix)] .+= rate
            # @inbounds du[sname(pn, is_ix)] .-= rate

            for j in os_ix
                du[sname(pn, j)] += rate
            end
            for j in is_ix
                du[sname(pn, j)] -= rate
            end

            # actually for net change, must look at each input place and see if there is a corresponding
            # output place, and if so ignore it (cancels out)
                        
        end
        return du
    end
    return f
end

ff = vectorfield(sir)
du = LVector(S=0., I=0., R=0.)
ff(du, u0, p, 0.0)

ff2 = vectorfield2(sir)
du2 = LVector(S=0., I=0., R=0.)
ff2(du2, u0, p, 0.0)

using Debugger
