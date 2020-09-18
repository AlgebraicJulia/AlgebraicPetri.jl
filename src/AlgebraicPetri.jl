""" Computing in the category of finite sets and Petri cospans
"""
module AlgebraicPetri
export PetriCospanOb, PetriFunctor, PetriCospan

using Catlab.GAT
using Catlab.Theories: BiproductCategory
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.Theories
using Catlab.CategoricalAlgebra.Limits
using Catlab.CategoricalAlgebra.FinSets
using AutoHashEquals

using Catlab.CategoricalAlgebra.FinSets: FinSet
import Catlab.Theories: dom, codom, id, compose, ⋅, ∘, otimes, ⊗, munit,
                        braid, σ, mcopy, Δ, mmerge, ∇, create, □, delete, ◊,
                        pair, copair, proj1, proj2, coproj1, coproj2

include("Types.jl")

""" Finite ordinal (natural number)

An object in the category of Open Petri Nets.
"""
@auto_hash_equals struct PetriCospanOb
    n::Int
end
PetriCospanOb(X::FinSet) = PetriCospanOb(length(X))
FinSet(X::PetriCospanOb) = FinSet(length(X))
# TODO: Extend to a type parameter for arbitrary petri net types (petri net, reaction net, labelled, etc)
EmptyPetriNet(X::PetriCospanOb) = PetriNet(length(X))

Base.iterate(X::PetriCospanOb, args...) = iterate(iterable(X), args...)
Base.length(X::PetriCospanOb) = length(iterable(X))

iterable(X::PetriCospanOb) = 1:X.n

struct PetriDecorator <: AbstractFunctor end
struct PetriLaxator <: AbstractLaxator end

""" Petri Functor

A functor from FinSet to Petri defined as a PetriDecorator and a PetriLaxator
"""
const PetriFunctor = LaxMonoidalFunctor{PetriDecorator, PetriLaxator}

@noinline id(::Type{PetriFunctor}) = PetriFunctor(PetriDecorator(), PetriLaxator())

""" AlgebraicPetri.PetriDecorator(n::FinSet)

A functor from FinSet to Set has an objects part, which given an object n in
FinSet (a natural number) should return a representation of F(n)::Set, sets can
be represented as a predicate that takes an element and returns true if the
element is in the set. Here we take any julia value and test whether it is a
Petri net on n states.
"""
function (pd::PetriDecorator)(n::FinSet)
    return p -> typeof(p) <: AbstractPetriNet && ns(p) == length(n)
end

""" AlgebraicPetri.PetriDecorator(f::FinFunction)

A functor from FinSet to Set has a hom part, which given a hom f in FinSet
(a function n::Int->m::Int) should return a representation of F(f)::F(n)->F(m),
here we implement this as a function that takes a Petri net of size n to a Petri
net of size m, such that the transitions are mapped appropriately.
"""
function (pd::PetriDecorator)(f::FinFunction) # TODO: MOVE TO AbstractPetriNet
    return (p::PetriNet) -> PetriNet(length(codom(f)), [(map(f, inputs(p, t)), map(f, outputs(p, t))) for t in 1:nt(p)]...)
end

""" AlgebraicPetri.PetriLaxator(p::Petri.Model, q::Petri.Model)

The laxitor takes a pair of decorations and returns the coproduct decoration
For Petri nets, this encodes the idea that you shift the states of q up by the
number of states in p.
"""
function (l::PetriLaxator)(p::PetriNet, q::PetriNet) # TODO: MOVE TO AbstractPetriNet
    return PetriNet(ns(p)+ns(q), [(inputs(p, t), outputs(p, t)) for t in 1:nt(p)]...,
        [(map(x->x+ns(p), inputs(q, t)), (map(x->x+ns(p), outputs(q, t)))) for t in 1:nt(q)]...)
end

""" Petri Cospan

A morphism in the category of Open Petri Nets defined as a decorated cospan with
a [`PetriFunctor`](@ref) as the decorator which maps the category of finite
ordinals to the category Petri and an AbstractPetriNet as the decoration
"""
const PetriCospan = DecoratedCospan{PetriFunctor, AbstractPetriNet}

""" AlgebraicPetri.PetriCospan(l::Vector{Int}, m::AbstractPetriNet, r::Vector{Int})

A constructor for Petri Cospans where `l` is a vector of the input states from
AbstractPetriNet `m`, and `r` is a vector of the output states from AbstractPetriNet `m`

Constructs the cospan: l → m ← r
"""
function (::Type{PetriCospan})(l::AbstractVector, p::AbstractPetriNet, r::AbstractVector)
    return PetriCospan(Cospan(FinFunction(l, ns(p)),
                              FinFunction(r, ns(p))),
                       id(PetriFunctor), p)
end

@instance BiproductCategory(PetriCospanOb, PetriCospan) begin
    dom(f::PetriCospan) = PetriCospanOb(dom(left(f)))
    codom(f::PetriCospan) = PetriCospanOb(dom(right(f)))

    compose(p::PetriCospan, q::PetriCospan) = begin
        # reimplementation of pushout of Span{FinSetFunc, FinSetFun}
        # to save the value of coeq
        f, g = right(p), left(q)
        ι1, ι2 = coproduct(codom(f), codom(g))
        π = proj(coequalizer(f⋅ι1, g⋅ι2))
        f′, g′ = ι1⋅π, ι2⋅π
        composite = Cospan(left(p)⋅f′, right(q)⋅g′)
        dpuq = decorator(p).L(decoration(p), decoration(q))
        return PetriCospan(composite, decorator(p), decorator(p).F(π)(dpuq))
    end

    id(X::PetriCospanOb) = PetriCospan(
        Cospan(id(FinSet(X)), id(FinSet(X))),
        id(PetriFunctor),
        EmptyPetriNet(X))

    otimes(X::PetriCospanOb, Y::PetriCospanOb) = PetriCospanOb(length(X) + length(Y))

    otimes(f::PetriCospan, g::PetriCospan) = begin
        fl, fr = left(f), right(f)
        gl, gr = left(g), right(g)
        # TODO: Replace with universal properties once implemented
        PetriCospan(
            Cospan(
                FinFunction(x->x > length(dom(fl)) ? gl(x-length(dom(fl)))+length(codom(fl)) : fl(x),
                               FinSet(length(dom(fl)) + length(dom(gl))),
                               FinSet(length(codom(fl)) + length(codom(gl)))),
                FinFunction(x->x > length(dom(fr)) ? gr(x-length(dom(fr)))+length(codom(fr)) : fr(x),
                               FinSet(length(dom(fr)) + length(dom(gr))),
                               FinSet(length(codom(fr)) + length(codom(gr))))
            ), decorator(f), decorator(f).L(decoration(f), decoration(g)))
    end

    munit(::Type{PetriCospanOb}) = PetriCospanOb(0)

    braid(X::PetriCospanOb, Y::PetriCospanOb) = begin
        Z = otimes(X, Y)
        PetriCospan(
            Cospan(
                id(FinSet(Z)),
                FinFunction(vcat(length(X)+1:length(Z), iterable(X)), length(Z), length(Z))
            ), id(PetriFunctor), EmptyPetriNet(Z))
    end

    mcopy(X::PetriCospanOb) = PetriCospan(
        Cospan(
            id(FinSet(X)),
            FinFunction(vcat(iterable(X),iterable(X)), 2*length(X), length(X))
        ), id(PetriFunctor), EmptyPetriNet(X))

    mmerge(X::PetriCospanOb) = PetriCospan(
        Cospan(
            FinFunction(vcat(iterable(X),iterable(X)), 2*length(X), length(X)),
            id(FinSet(X))
        ), id(PetriFunctor), EmptyPetriNet(X))

    create(X::PetriCospanOb) = PetriCospan(
        Cospan(FinFunction(Int[], 0, length(X)), id(FinSet(X))),
        id(PetriFunctor), EmptyPetriNet(X))

    delete(X::PetriCospanOb) = PetriCospan(
        Cospan(id(FinSet(X)), FinFunction(Int[], 0, length(X))),
        id(PetriFunctor), EmptyPetriNet(X))

    pair(f::PetriCospan, g::PetriCospan) = compose(mcopy(dom(f)), otimes(f, g))
    copair(f::PetriCospan, g::PetriCospan) = compose(otimes(f, g), mmerge(codom(f)))

    proj1(A::PetriCospanOb, B::PetriCospanOb) = otimes(id(A), delete(B))
    proj2(A::PetriCospanOb, B::PetriCospanOb) = otimes(delete(A), id(B))

    coproj1(A::PetriCospanOb, B::PetriCospanOb) = otimes(id(A), create(B))
    coproj2(A::PetriCospanOb, B::PetriCospanOb) = otimes(create(A), id(B))
end

# reenable after refactor
include("Epidemiology.jl")

end
