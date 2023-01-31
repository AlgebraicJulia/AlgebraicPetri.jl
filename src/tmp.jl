using AlgebraicPetri
using Graphviz_jll

sir = LabelledPetriNet([:S,:I,:R],
    :inf => ((:S,:I) => (:I,:I)),
    :rec => (:I => :R)
)

Graph(sir)

macro test_macro(pn)
    num_t = :(nt($(pn)))
    num_s = :(ns($(pn)))
    # println("our pn has $(esc(num_t)) transitions and $(esc(num_s)) places")
    quote
        ($(num_t), $(num_s))
        println("our pn has $(num_t) transitions and $(num_s) places")
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