using Catlab
using Catlab.Theories
using Catlab.Graphics
using Catlab.WiringDiagrams

display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

"""

Vocabulary
Nouns:
======
    Mature cathepsin – cathepsin form able to hydrolyze substrate or other cathepsins or to bind inhibitor
    Procathepsin – Precursor to mature cathepsin that can be activated
    Inactive – not capable of hydrolyzing substrate
    Inhibited – mature cathepsin and inhibitor bind to become new inhibited species
    Complex(X,Y) – A complex of two molecules

Verbs:
======
    Activation – removal of propeptide from procathepsin to convert it to mature cathepsin (can happen on its own, or can be mediated by another protease)
    Binding – two molecules/species become one Complex
    Hydrolyzing (deg) – breaking peptide bonds
    Cannabilizing – one cathepsin hydrolyzes another
    Distracting – inactive cathepsin and active cathepsin bind and hydrolysis occurs

"""
@theory Cathepsins{Ob,Hom,Cat,Inhibitor,Pro,Clx,Deg,Inhib} <: MonoidalCategory{Ob,Hom} begin
    Cat()::TYPE
    Inhibitor()::TYPE
    Pro(cat::Cat)::TYPE
    Clx(enzyme::Cat,substrate::Ob)::TYPE
    Deg(cat::Cat)::TYPE
    Inhib(cat::Cat,inhibitor::Inhibitor)::TYPE

    ob(X::Cat)::Ob

    # sub(X::Ob)::Ob
    # gen(cat(X))::Ob ⊣ (X::Ob)
    # inact(cat(X))::Ob ⊣ (X::Ob)
    # drug(X::Ob)::Ob

    activate(X)      :: (Pro(X) → X) ⊣ (X::Cat)
    inhibit(X,Y)  :: (X ⊗ Y  → Inh(X,Y)) ⊣ (X::Cat, Y::Inhibitor)
    bind(X,Y)        :: (X ⊗ Y  → Clx(X,Y)) ⊣ (X::Ob, Y::Ob)
    # ... add more possibilities ex. Cat, Cat
    unbind(X,Y)      :: (Clx(X,Y) → X ⊗ Y) ⊣ (X::Ob, Y::Ob)
    # ... add more possibilities ex. Cat, Cat
    hydrolyze(X, Y)   :: (Clx(X,Y)  → X ⊗ Deg(Y)) ⊣ (X::Ob, Y::Ob)
    # ... add more possibilities ex. Cat, Cat
    cannabilize(X,Y) :: (X ⊗ Y  → X ⊗ Deg(Y)) ⊣ (X::Cat, Y::Cat)
    # ... add more possibilities ex. Cat, Cat
    distract(X,Y)    :: (X ⊗ Y → X ⊗ Deg(cat(Y))) ⊣ (X::Cat, Y::Pro)
    distract(X,Y)    :: (X ⊗ Y → X ⊗ Deg(cat(Y))) ⊣ (X::Cat, Y::Inhib)
    # # # TODO: What are the products of distraction?
    # cat(sub(x)) == sub(cat(x))
    # cat(sub(x)) == munit()
    # cat(drug(x)) == drug(cat(x))
    # cat(drug(x)) == munit()
    # cannabilize(X,Y) == hydrolyze(X,Y) # we should probably distinguish proteins and cathepsins as different kinds of objects
    distract(X,Y) == (bind(X,Y) ⋅ hydrolyze(X,Y)) ⊣ (X::Cat, Y::Pro)
    distract(X,Y) == (bind(X,Y) ⋅ hydrolyze(X,Y)) ⊣ (X::Cat, Y::Inhib)
end

@syntax FreeCathepsins{ObExpr,HomExpr,GATExpr,GATExpr,GATExpr,GATExpr,GATExpr,GATExpr,GATExpr,GATExpr} Cathepsins begin
    otimes(A::Ob, B::Ob) = associate_unit(new(A,B), munit)
    otimes(f::Hom, g::Hom) = associate(new(f,g))
    compose(f::Hom, g::Hom) = associate(new(f,g; strict=true))
end

clx(x::Symbol, y::Symbol) = Symbol(x,y)
bind(p::Ports{FreeCathepsins.Hom}, q::Ports{FreeCathepsins.Hom}) = begin
    a = length(p.ports) > 1 ? otimes(p.ports) : p.ports[1]
    b = length(q.ports) > 1 ? otimes(q.ports) : q.ports[1]
    c = clx(a,b)
    return singleton_diagram(Box(:bind, Ports([p.ports..., q.ports...]), Ports([c])))
end
bind(p::Ports{Cathepsins}, q::Ports{Cathepsins}) = begin
    a = length(p.ports) > 1 ? otimes(p.ports) : p.ports[1]
    b = length(q.ports) > 1 ? otimes(q.ports) : q.ports[1]
    c = clx(a,b)
    return singleton_diagram(Box(:bind, Ports([p.ports..., q.ports...]), Ports([c])))
end
clx(p::Ports{FreeCathepsins.Hom}, q::Ports{FreeCathepsins.Hom}) = begin
    z = collect(zip(p.ports, q.ports))
    Ports{FreeCathepsins.Hom}([clx(a,b) for (a,b) in z])
end



K = Ob(FreeCathepsins.Ob, :K)
E = Ob(FreeCathepsins.Ob, :E)

display_wd(bind(K,E))