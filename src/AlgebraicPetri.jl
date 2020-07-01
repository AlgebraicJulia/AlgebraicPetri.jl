""" Computing in the category of finite sets and Petri cospans
"""
module AlgebraicPetri
export PetriCospanOb, PetriFunctor, PetriCospan

using Catlab.GAT
using Catlab.Theories: BiproductCategory
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets
using Petri
using AutoHashEquals

import Catlab.Theories: dom, codom, id, compose, ⋅, ∘, otimes, ⊗, munit,
                        braid, σ, mcopy, Δ, mmerge, ∇, create, □, delete, ◊,
                        pair, copair, proj1, proj2, coproj1, coproj2

@auto_hash_equals struct PetriCospanOb
    n::Int
end
Base.eachindex(X::PetriCospanOb) = 1:X.n

struct PetriDecorator <: AbstractFunctor end
struct PetriLaxitor <: AbstractLaxitor end

const PetriFunctor = LaxMonoidalFunctor{PetriDecorator, PetriLaxitor}

id(::Type{PetriFunctor}) = PetriFunctor(PetriDecorator(), PetriLaxitor())

""" we follow FinOrd's lead in taking the skeleton of FinSet so we work with
Petri nets where the state space is a single Int n rather than the range 1:n.
A functor from FinOrd to Set has an objects part, which given an object n in
FinOrd (a natural number) should return a representation of F(n)::Set, sets can
be represented as a predicate that takes an element and returns true if the
element is in the set. Here we take any julia value and test whether it is a
Petri net on n states.
"""
function (pd::PetriDecorator)(n::FinOrd)
    return p -> typeof(p) <: Petri.Model && length(p.S) == n.n
end

""" mapping a FinOrdFunc over the transitions of a Petri net means renumbering
the states incident to a transition where the transitions of a petri net are
stored as a list of tuples of lists. This could be easily adapted for storing
transitions in a list of pairs of multisets.
"""
map(f::Function, ts::Vector{Tuple{Vector{Int}, Vector{Int}}}) = [(f.(t[1]), f.(t[2])) for t in ts]

""" A functor from FinOrd to Set has a hom part, which given a hom f in FinOrd
(a function n::Int->m::Int) should return a representation of F(f)::F(n)->F(m), 
here we implement this as a function that takes a Petri net of size n to a Petri
net of size m, such that the transitions are mapped appropriately.
"""
function (pd::PetriDecorator)(f::FinOrdFunction)
    return (p::Petri.Model) -> Petri.Model(collect(1:codom(f).n), map(x->f(x), p.Δ))
end

""" the laxitor takes a pair of decorations and returns the coproduct decoration
For Petri nets, it encodes the idea that you shift the states of q up by the
number of states in p.
"""
function (l::PetriLaxitor)(p::Petri.Model, q::Petri.Model)
    return Petri.Model(collect(1:(length(p.S)+length(q.S))),
                       vcat(p.Δ, map(x->x+length(p.S), q.Δ))
                      )
end

const PetriCospan = DecoratedCospan{PetriFunctor, Petri.Model}

NullPetri(n::Int) = Petri.Model(collect(1:n), Vector{Tuple{Vector{Int}, Vector{Int}}}())

@instance BiproductCategory(PetriCospanOb, PetriCospan) begin
    dom(f::PetriCospan) = PetriCospanOb(dom(left(f)).n)
    codom(f::PetriCospan) = PetriCospanOb(dom(right(f)).n)

    compose(p::PetriCospan, q::PetriCospan) = begin
        # reimplementation of pushout of Span{FinOrdFunc, FinOrdFun}
        # to save the value of coeq
        f, g = right(p), left(q)
        coprod = coproduct(codom(f), codom(g))
        ι1, ι2 = left(coprod), right(coprod)
        coeq = coequalizer(f⋅ι1, g⋅ι2)
        f′, g′ = ι1⋅coeq, ι2⋅coeq
        composite = Cospan(left(p)⋅f′, right(q)⋅g′)
        dpuq = decorator(p).L(decoration(p), decoration(q))
        return PetriCospan(composite, decorator(p), decorator(p).F(coeq)(dpuq))
    end

    id(X::PetriCospanOb) = PetriCospan(
        Cospan(id(FinOrd(X.n)), id(FinOrd(X.n))),
        id(PetriFunctor),
        NullPetri(X.n)
    )

    otimes(X::PetriCospanOb, Y::PetriCospanOb) = PetriCospanOb(X.n + Y.n)

    otimes(f::PetriCospan, g::PetriCospan) = begin
        fl, fr = left(f), right(f)
        gl, gr = left(g), right(g)
        PetriCospan(
            Cospan(
                FinOrdFunction(x->x > dom(fl).n ? gl(x-dom(fl).n)+codom(fl).n : fl(x),
                               FinOrd(dom(fl).n + dom(gl).n),
                               FinOrd(codom(fl).n + codom(gl).n)),
                FinOrdFunction(x->x > dom(fr).n ? gr(x-dom(fr).n)+codom(fr).n : fr(x),
                               FinOrd(dom(fr).n + dom(gr).n),
                               FinOrd(codom(fr).n + codom(gr).n))
            ),
            decorator(f), decorator(f).L(decoration(f), decoration(g))
        )
    end

    munit(::Type{PetriCospanOb}) = PetriCospanOb(0)

    braid(X::PetriCospanOb, Y::PetriCospanOb) = begin
        Z = otimes(X, Y)
        PetriCospan(
            Cospan(
                id(FinOrd(Z.n)),
                FinOrdFunction(x->vcat(X.n+1:Z.n, 1:X.n)[x], Z.n, Z.n)
            ), id(PetriFunctor), NullPetri(Z.n)
        )
    end

    mcopy(X::PetriCospanOb) = PetriCospan(
        Cospan(
            id(FinOrd(X.n)),
            FinOrdFunction(x->vcat(1:X.n,1:X.n)[x], 2*X.n, X.n)
        ), id(PetriFunctor), NullPetri(X.n)
    )

    mmerge(X::PetriCospanOb) = PetriCospan(
        Cospan(
            FinOrdFunction(x->vcat(1:X.n,1:X.n)[x], 2*X.n, X.n),
            id(FinOrd(X.n))
        ), id(PetriFunctor), NullPetri(X.n)
    )

    create(X::PetriCospanOb) = PetriCospan(
        Cospan(FinOrdFunction(x->[][x], 0, X.n), id(FinOrd(X.n))),
        id(PetriFunctor), NullPetri(X.n)
    )

    delete(X::PetriCospanOb) = PetriCospan(
        Cospan(id(FinOrd(X.n)), FinOrdFunction(x->[][x], 0, X.n)),
        id(PetriFunctor), NullPetri(X.n)
    )

    pair(f::PetriCospan, g::PetriCospan) = compose(mcopy(dom(f)), otimes(f, g))
    copair(f::PetriCospan, g::PetriCospan) = compose(otimes(f, g), mmerge(codom(f)))

    proj1(A::PetriCospanOb, B::PetriCospanOb) = otimes(id(A), delete(B))
    proj2(A::PetriCospanOb, B::PetriCospanOb) = otimes(delete(A), id(B))

    coproj1(A::PetriCospanOb, B::PetriCospanOb) = otimes(id(A), create(B))
    coproj2(A::PetriCospanOb, B::PetriCospanOb) = otimes(create(A), id(B))
end

end
