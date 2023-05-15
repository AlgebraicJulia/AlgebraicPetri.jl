module TestRewriting

using AlgebraicPetri
using AlgebraicRewriting
using Catlab, Catlab.CategoricalAlgebra
using Test

const homomorphism = CategoricalAlgebra.homomorphism
const is_isomorphic = CategoricalAlgebra.is_isomorphic

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

# seir test
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

@test is_isomorphic(seir, rewrite_match(Rule{:DPO}(L,R), m))

# seirs test
Lg = LabelledReactionNet{Float64, Float64}(
  [:S=>100,:I=>1],
)
Rg = LabelledReactionNet{Float64, Float64}(
  [:S=>100,:I=>1], 
  (:sus,.1)=>(:I=>:S)
)

L = homomorphism(Lg, Lg)
R = homomorphism(Lg, Rg)
m = homomorphism(Lg, seir)

seirs = rewrite_match(Rule{:DPO}(L,R), m)

@test nt(seirs) == nt(seir) + 1
@test ns(seirs) == ns(seir)

sus_ix = only(incident(seirs, :sus, :tname))
sus_input = only(seirs[incident(seirs, sus_ix, :it), :is])
sus_output = only(seirs[incident(seirs, sus_ix, :ot), :os])

@test seirs[sus_input,:sname] == :I
@test seirs[sus_output,:sname] == :S

# seird test

Lg = LabelledReactionNet{Float64, Float64}(
  [:I=>1],
)
Rg = LabelledReactionNet{Float64, Float64}(
  [:I=>1,:D=>0],
  (:die,.1)=>(:I=>:D)
)

L = homomorphism(Lg, Lg)
R = homomorphism(Lg, Rg)
m = homomorphism(Lg, seir)

seird = rewrite_match(Rule{:DPO}(L,R), m)

@test nt(seird) == nt(seir) + 1
@test ns(seird) == ns(seir) + 1

die_ix = only(incident(seird, :die, :tname))
die_input = only(seird[incident(seird, die_ix, :it), :is])
die_output = only(seird[incident(seird, die_ix, :ot), :os])

@test seird[die_input,:sname] == :I
@test seird[die_output,:sname] == :D

end
