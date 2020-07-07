module DynamicRates
using Petri
import Petri: funcindex!

# struct Parameters{T}
#     x::T
# end

valueat(x::Number, t) = x
valueat(f::Function, t) = f(t)

function toODE(m::Model, p)
    S = m.S
    T = m.Δ
    ϕ = Dict()
    f(du, u, p,t) = begin
        for k in keys(T)
            ins = first(getindex(T, k))
            ϕ[k] = reduce((x,y)->x*getindex(u,y)/getindex(ins,y), keys(ins); init=valueat(getindex(p, k),t))
        end
        for s in S
            du[s] = 0
        end
        for k in keys(T)
            ins = first(getindex(T, k))
            outs = last(getindex(T, k))
            for s in keys(ins)
                funcindex!(du, s, -, getindex(ϕ, k) * getindex(ins, s))
            end
            for s in keys(outs)
                funcindex!(du, s, +, getindex(ϕ, k) * getindex(outs, s))
            end
        end
        return du
    end
    return f
end

end
