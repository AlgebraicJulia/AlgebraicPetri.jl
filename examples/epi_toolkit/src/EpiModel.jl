module EpiModel

# The available functions, all have necessary docstrings
export EpiNet, EpiNetSchema
using AlgebraicPetri
using Semagrams

using JSON

using Catlab
using Catlab.CategoricalAlgebra
import Catlab.CategoricalAlgebra: migrate!
import AlgebraicPetri: LabelledReactionNet

@present TheoryEpiNet <: TheoryPetriNet begin
  Rate::AttrType
  Population::AttrType
  Name::AttrType

  tname::Attr(T, Name)
  sname::Attr(S, Name)

  rate::Attr(T, Rate)
  population::Attr(S, Population)
end

@abstract_acset_type AbstractEpiNet
@acset_type EpiNetUntyped(TheoryEpiNet, index=[:it, :is, :ot, :os]) <: AbstractEpiNet
const EpiNet{R,C} = EpiNetUntyped{R,C,Symbol}

@semagramschema EpiNetSchema(TheoryEpiNet) begin
  @box S Circle label=:sname
  @box T Square label=:tname
  @wire I(is, it)
  @wire O(ot, os)
  @data Name Stringlike
  @data Rate Numeric
  @data Population Numeric
end

function migrate!(pn::LabelledReactionNet, epn::EpiNet)
  migrate!(pn, epn,
           Dict(:S=>:S, :T=>:T, :I=>:I, :O=>:O, :Rate => :Rate,
                :Concentration => :Population, :Name => :Name),
           Dict(:is=>:is, :it=>:it, :os=>:os, :ot=>:ot,
                :rate=>:rate, :concentration=>:population,
                :tname=>:tname, :sname=>:sname))
end

LabelledReactionNet(epn::EpiNet{R,C}) where {R,C} = begin
  pn = LabelledReactionNet{R,C}()
  migrate!(pn, epn)
  pn
end
end

