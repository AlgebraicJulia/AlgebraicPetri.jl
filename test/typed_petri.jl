using Test

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.CategoricalAlgebra, Catlab.Programs

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

sird_uwd = @relation () where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S,I,I,I) # inf
  disease(I,R) # recover
  disease(I,D) # die
end

typed_sird = add_params(
  oapply_typed(infectious_ontology, sird_uwd, [:inf, :recover, :die]),
  Dict(:S => 1.0, :I => 1.0, :R => 0.0, :D => 0.0),
  Dict(:inf => 0.5, :recover => 1.0, :die => 0.2)
)

@test ns(typed_sird.dom) == 4
@test nt(typed_sird.dom) == 3

# Quarantine model.

quarantine_uwd = @relation () where (Q::Pop, NQ::Pop) begin
  strata(Q,NQ) # enter quarantine
  strata(NQ,Q) # exit quarantine
end

typed_quarantine = add_params(
  oapply_typed(infectious_ontology, quarantine_uwd, [:exit_Q, :enter_Q]),
  Dict(:Q => 1.0, :NQ => 1.0),
  Dict(:exit_Q => 0.5, :enter_Q => 0.5)
)

typed_quarantine_aug = add_reflexives(
  typed_quarantine,
  [[:disease], [:disease, :infect]],
  infectious_ontology
)

typed_sird_aug = add_reflexives(
  typed_sird,
  [[:strata], [:strata], [:strata], []],
  infectious_ontology
)

stratified = typed_product(typed_quarantine_aug, typed_sird_aug)
@test ns(dom(stratified)) == 8
@test nt(dom(stratified)) == 6 + 4 + 1
@test dom(stratified)[1,:sname] isa Tuple{Symbol, Symbol}
@test dom(stratified)[1,:concentration] isa Tuple{Float64, Float64}

# Age-stratified model.

typed_age = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect,
                                    [:Young, :Mid, :Old])
@test ns(dom(typed_age)) == 3
@test nt(dom(typed_age)) == 9

typed_age_aug = add_reflexives(
  typed_age,
  repeat([[:disease]], 3),
  infectious_ontology
)

typed_sird = oapply_typed(infectious_ontology, sird_uwd, [:inf, :recover, :die])

stratified = typed_product(typed_age_aug, typed_sird)
@test ns(dom(stratified)) == 3*4
@test nt(dom(stratified)) == 3*2 + 9*1
