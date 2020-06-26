module AlgebraicPetri

using Catlab.GAT
using Catlab.Theories: BiproductCategory
using Catlab.CategoricalAlgebra.ShapeDiagrams
using Catlab.CategoricalAlgebra.FinSets
using Petri

import Catlab.Theories: dom, codom, id, compose, ⋅, ∘

struct PetriCospan
    f::Cospan
end

@instance BiproductCategory(FinOrd, PetriCospan) begin
    dom(f::PetriCospan) = FinOrd(dom(left(f.f)))
    codom(f::PetriCospan) = FinOrd(dom(right(f.f)))

    compose(f::PetriCospan, g::PetriCospan) = begin
        Nothing
    end

    id(X::FinOrd) = Nothing

    otimes(X::FinOrd, Y::FinOrd) = FinOrd(X.n + Y.n)

    otimes(f::PetriCospan, g::PetriCospan) = begin
        Nothing
    end

    munit(::Type{FinOrd}) = FinOrd(0)

    braid(X::FinOrd, Y::FinOrd) = begin
        Nothing
    end

    mcopy(X::FinOrd) = begin
        Nothing
    end

    mmerge(X::FinOrd) = begin
        Nothing
    end

    create(X::FinOrd) = begin
        Nothing
    end

    delete(X::FinOrd) = begin
        Nothing
    end

    pair(f::PetriCospan, g::PetriCospan) = compose(mcopy(dom(f)), otimes(f, g))
    copair(f::PetriCospan, g::PetriCospan) = compose(otimes(f, g), mmerge(codom(f)))

    proj1(A::FinOrd, B::FinOrd) = otimes(id(A), delete(B))
    proj2(A::FinOrd, B::FinOrd) = otimes(delete(A), id(B))

    coproj1(A::FinOrd, B::FinOrd) = otimes(id(A), create(B))
    coproj2(A::FinOrd, B::FinOrd) = otimes(create(A), id(B))
end

end
