module TestTypedPetri

using Test

using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.CategoricalAlgebra, Catlab.Programs

# SIRD model.

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

sird_uwd = @relation (S,I,R,D) where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S,I,I,I) # inf
  disease(I,R) # recover
  disease(I,D) # die
end

typed_sird_no_params = oapply_typed(infectious_ontology, sird_uwd,
                                    [:inf, :recover, :die])
@test ns(dom(typed_sird_no_params)) == 4
@test nt(dom(typed_sird_no_params)) == 3

# Backwards compatibility: allow UWDs with no outer ports.
sird_uwd_no_outer = @relation () where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S,I,I,I) # inf
  disease(I,R) # recover
  disease(I,D) # die
end
@test oapply_typed(infectious_ontology, sird_uwd_no_outer,
                   [:inf, :recover, :die]) == typed_sird_no_params

typed_sird = add_params(
  oapply_typed(infectious_ontology, sird_uwd, [:inf, :recover, :die]),
  Dict(:S => 1.0, :I => 1.0, :R => 0.0, :D => 0.0),
  Dict(:inf => 0.5, :recover => 1.0, :die => 0.2)
)
@test ns(dom(typed_sird)) == 4
@test nt(dom(typed_sird)) == 3

# SIRD-with-quarantine model.

quarantine_uwd = @relation (Q,NQ) where (Q::Pop, NQ::Pop) begin
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
fl_strat = dom(flatten_labels(stratified))
tp = typed_product([typed_quarantine_aug, typed_sird_aug])
fl_tp = dom(flatten_labels(tp))
@test fl_strat == fl_tp
@test ns(dom(stratified)) == 8
@test nt(dom(stratified)) == 6 + 4 + 1
@test dom(stratified)[1,:sname] isa Tuple{Symbol, Symbol}
@test dom(stratified)[1,:concentration] isa Tuple{Float64, Float64}

# Age-stratified model.

typed_age = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect,
                                    [:Yng, :Mid, :Old])
@test ns(dom(typed_age)) == 3
@test nt(dom(typed_age)) == 9

typed_age_aug = add_reflexives(
  typed_age,
  repeat([[:disease]], 3),
  infectious_ontology
)

stratified = typed_product(typed_age_aug, typed_sird_no_params)
@test ns(dom(stratified)) == 3*4
@test nt(dom(stratified)) == 3*2 + 9*1

contact_mat = rand(Float64, (3, 3))
typed_age = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect,
                                    [:Yng, :Mid, :Old], [0.3, 0.4, 0.3],
                                    contact_mat, codom_net=codom(typed_sird))
@test dom(typed_age)[2, :concentration] == 0.4
@test dom(typed_age)[1, :rate] == contact_mat[1,1]

typed_age_aug = add_reflexives(
  typed_age,
  repeat([[:disease]], 3),
  infectious_ontology
)

stratified = typed_product(typed_age_aug, typed_sird)

end
