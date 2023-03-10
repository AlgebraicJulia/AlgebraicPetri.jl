using AlgebraicPetri
using AlgebraicRewriting

sir = LabelledReactionNet{Float64, Float64}(
    [:S=>100, :I=>1, :R=>0], 
    (:inf,.03)=>((:S,:I)=>(:I,:I)),
    (:rec,.25)=>(:I=>:R)
)

seir = LabelledReactionNet{Float64, Float64}(
    [:S=>100,:I=>1,:E=>1,:R=>0],
    (:inf,.03)=>((:S,:I)=>(:I,:I)),
    (:rec,.25)=>(:I=>:R),
    (:inc,.1)=>(:E=>:I),
    (:exp,.1)=>((:S,:I)=>(:E,:I))
)

Rg = LabelledReactionNet{Float64, Float64}(
    [:S=>100,:I=>1,:E=>1],
    (:exp,.1)=>((:S, :I)=>(:E, :I)),
    (:inc,.1)=>(:E => :I)
)
Lg = LabelledReactionNet{Float64, Float64}(
    [:S=>100,:I=>1],
)

L = homomorphism(Lg, Lg)
R = homomorphism(Lg, Rg)
m = homomorphism(Lg, sir)

is_isomorphic(seir, rewrite_match(Rule{:DPO}(L,R), m))