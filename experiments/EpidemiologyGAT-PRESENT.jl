using AlgebraicPetri
using Petri
using Catlab
using Catlab.GAT
using Catlab.Theories
using Catlab.Programs
using Catlab.WiringDiagrams

using Catlab.CategoricalAlgebra.ShapeDiagrams # remove later
using Catlab.Graphics # remove later

import Catlab.Theories: dom, codom, id, compose, ⋅, ∘, otimes, ⊗, munit,
                        braid, σ, mcopy, Δ, mmerge, ∇, create, □, delete, ◊,
                        pair, copair, proj1, proj2, coproj1, coproj2
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

@theory BiproductCategory(Ob,Hom) => EpidemiologyC(Ob,Hom) begin
  spontaneous(A::Ob, B::Ob)::(A → B)
  transmission(A::Ob, B::Ob)::(A ⊗ B → B)
  exposure(A::Ob, B::Ob, C::Ob)::(A ⊗ B → C ⊗ B)
end

@instance EpidemiologyC(PetriCospanOb, PetriCospan) begin
  @import dom, codom, compose, id, otimes, munit, braid, mcopy, mmerge, create, delete, pair, copair, proj1, proj2, coproj1, coproj2
  spontaneous(A::PetriCospanOb, B::PetriCospanOb) = PetriCospan([1], Petri.Model(1:2, [(Dict(1=>1), Dict(2=>1))]), [2])
  transmission(A::PetriCospanOb, B::PetriCospanOb) = PetriCospan([1,2], Petri.Model(1:2, [(Dict(1=>1, 2=>1), Dict(2=>2))]), [2])
  exposure(A::PetriCospanOb, B::PetriCospanOb, C::PetriCospanOb) = PetriCospan([1, 2], Petri.Model(1:3, [(Dict(1=>1, 2=>1), Dict(3=>1, 2=>1))]), [3, 2])
end

@syntax FreeEpidemiologyC(ObExpr,HomExpr) EpidemiologyC begin
    otimes(A::Ob, B::Ob) = associate_unit(new(A,B), munit)
    otimes(f::Hom, g::Hom) = associate(new(f,g))
    compose(f::Hom, g::Hom) = associate(new(f,g; strict=true))

    pair(f::Hom, g::Hom) = Δ(dom(f)) ⋅ (f ⊗ g)
    copair(f::Hom, g::Hom) = (f ⊗ g) ⋅ ∇(codom(f))
    proj1(A::Ob, B::Ob) = id(A) ⊗ ◊(B)
    proj2(A::Ob, B::Ob) = ◊(A) ⊗ id(B)
    coproj1(A::Ob, B::Ob) = id(A) ⊗ □(B)
    coproj2(A::Ob, B::Ob) = □(A) ⊗ id(B)
end

spontaneous(a::Ports{EpidemiologyC}, b::Ports{EpidemiologyC}) = begin
    return singleton_diagram(Box(:spontaneous, Ports([a.ports...]), Ports([b.ports...])))
end

transmission(a::Ports{EpidemiologyC}, b::Ports{EpidemiologyC}) = begin
    return singleton_diagram(Box(:transmission, Ports([a.ports...,b.ports...]), Ports([b.ports...])))
end

exposure(a::Ports{EpidemiologyC}, b::Ports{EpidemiologyC}, c::Ports{EpidemiologyC}) = begin
    dump(a)
    return singleton_diagram(Box(:exposure, Ports([a.ports...,b.ports...]), Ports([c.ports...,b.ports...])))
end

mcopy(A::Ports{EpidemiologyC}, n::Int) = implicit_mcopy(A, n)
mmerge(A::Ports{EpidemiologyC}, n::Int) = implicit_mmerge(A, n)

@present BasicEpi(FreeEpidemiologyC) begin
  S::Ob
  E::Ob
  I::Ob
  R::Ob
  D::Ob
end

S,E,I,R,D = generators(BasicEpi)

F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(map(x->x=>PetriCospanOb(1), generators(BasicEpi))))

sir = transmission(S,I) ⋅ spontaneous(I,R)

sir = @program BasicEpi (s::S, i::I) begin
    i_2 = transmission{S,I}(s, i)
    return spontaneous{I,R}(i)
end
F(to_hom_expr(FreeEpidemiologyC, sir))

F(sir)


s = PetriCospanOb(1)
i = PetriCospanOb(1)
r = PetriCospanOb(1)

sei = expose ⋅ (illness ⊗ id(I)) ⋅ ∇(I)

seir = sei ⋅ recover

seird = sei ⋅ Δ(I) ⋅ (death ⊗ recover)

display_wd(seird)

# new coexist stuff

@present Coexist(FreeEpidemiologyC) begin
    S::Ob
    E::Ob
    I1::Ob
    I2::Ob
    A::Ob
    R1::Ob
    R2::Ob
    D::Ob

    # exposes(a) := exposure(S,a,E)
    illness := spontaneous(E,I1)
    asymptomatic_illness := spontaneous(E,A)
    progression := spontaneous(I1,I2)
    death := spontaneous(I2,D)
    recover := spontaneous(I2,R1)
    asymptomatic_recover := spontaneous(A,R1)
    recover_progression := spontaneous(R1,R2)
end

S,E,I1,I2,A,R1,R2,D,illness,asymptomatic_illness,progression,death,recover,asymptomatic_recover,recover_progression = generators(Coexist)
F(ex) = functor((PetriCospanOb, PetriCospan), ex, generators=Dict(map(x->x=>PetriCospanOb(1), generators(Coexist))))

coexist = @program Coexist (s::S, e::E, i::I1, i2::I2, a::A, r::R1, r2::R2, d::D) begin
    e_2, i_2 = exposure{S,I1,E}(s, i)
    e_3, i2_2 = exposure{S,I2,E}(s, i2)
    e_4, a_2 = exposure{S,A,E}(s, a)
    e_5 = exposure{S,E,E}(s, e)
    e_all = [e_2, e_3, e_4, e_5]
    a_3 = asymptomatic_illness(e_all)
    a_all = [a_2, a_3]
    r_2 = asymptomatic_recover(a_all)
    i_3 = illness(e_all)
    i_all = [i_2, i_3]
    i2_3 = progression(i_all)
    i2_all = [i2_2, i2_3]
    d_2 = death(i2_all)
    r_3 = recover(i2_all)
    r_all = [r, r_2, r_3]
    r2_2 = recover_progression(r_all)
    r2_all = [r2, r2_2]
    d_all = [d, d_2]
    return s, e_all, i_all, i2_all, a_all, r_all, r2_all, d_all
end

display_wd(coexist)

to_hom_expr(FreeEpidemiologyC, coexist)
display_wd(to_hom_expr(FreeEpidemiologyC, coexist))