using Test

using Catlab.Programs, AlgebraicPetri, AlgebraicPetri.TypedPetri

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

sird_uwd = @relation () where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  infect(S,I,I,I)
  disease(I,R)
  disease(R,D)
end

typed_sird = oapply_typed(infectious_ontology, sird_uwd)

@test ns(typed_sird) == 4
@test nt(typed_sird) == 3
