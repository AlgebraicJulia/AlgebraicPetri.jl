module SubACSets
export mca

using DataStructures
using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphics
using Catlab.CategoricalAlgebra.Diagrams
using Catlab.CategoricalAlgebra.FreeDiagrams
using Catlab.ACSetInterface
import Catlab.CategoricalAlgebra.FinCats: FreeCatGraph, FinDomFunctor, collect_ob, collect_hom
using Catlab.Programs
using AlgebraicPetri

concatmap(f,xs) = mapreduce(f, vcat, xs; init =[])

function strip_names(p::AbstractPetriNet)
  map(p, Name = name -> nothing)
end

function strip_names(p::ACSetTransformation)
  init = NamedTuple([k=>collect(v) for (k,v) in pairs(components(p))])
  homomorphism(strip_names(dom(p)), strip_names(codom(p)), initial=init)
end

# the smaller acset is the "pattern"
function normalize_order(X::ACSet, Y::ACSet)
  size(X) ≤ size(Y) ? (X,Y) : (Y,X)
end

"""
    rm_cascade_subobj(X::ACSet, rm_subs)

Deletes parts from an ACSet in cascading fashion, e.g. deleting a vertex deletes its edges
rm_subs is a NamedTuple or Dict of parts to be removed.
"""
function rm_cascade_subobj(X::ACSet, rm_subs)
  S    = acset_schema(X)
  subs = Dict([k=>Set(parts(X,k)) for k in objects(S)])
  rm_subs = Dict([k=>Set(v) for (k,v) in pairs(rm_subs)])
  # println(keys(subs))
  while !isempty(rm_subs)
    # println(keys(rm_subs))
    curr_c = first(rm_subs)[1]
    if isempty(rm_subs[curr_c])
        delete!(rm_subs,curr_c)
    else
        # println(rm_subs[curr_c])
        curr_part = pop!(rm_subs[curr_c])
        if curr_part in subs[curr_c]
            delete!(subs[curr_c],curr_part)
            #
            # if isempty(subs[curr_c])
            #     delete!(subs,curr_c)
            # end
            # 
            for (f,c,d) in homs(S)
                if d==curr_c && c in keys(subs)
                    # println("Made it")
                    for test_part in subs[c]
                        if X[test_part,f] == curr_part
                            if c in keys(rm_subs)
                                push!(rm_subs[c], test_part)
                            else
                                rm_subs[c] = Set([test_part])
                            end
                        end
                    end
                end
            end
        end
        if isempty(rm_subs[curr_c])
            delete!(rm_subs,curr_c)
        end
    end
  end
  # return subs # dom(hom(Subobject(X,NamedTuple(subs))))
  # return Dict([k => collect(v) for (k,v) in pairs(subs)])
  return dom(hom(Subobject(X,NamedTuple(Dict([k => collect(v) for (k,v) in pairs(subs)])))))
end

"""
Defintion: let 𝐺: C → 𝐒et be a C-set, we define the _size_ of 𝐺 to be ∑_{c ∈ C}
|𝐺c|.  For example, under this definition, the size of:
  * a graph G is |GE| + |GV| (num edges + num vertices)
  * a Petri net P is |PT| + |PS| + |PI| + |PO| (num transitions + num species +
    num input arcs + num output arcs).
"""
function size(X::ACSet)
  foldl(+, [length(parts(X, oₛ)) for oₛ ∈ objects(acset_schema(X))])
end

function strip_attributes(p::ACSet; attributes::Vector{Symbol}=Symbol[])
  attributes = isempty(attributes) ? attrtypes(acset_schema(p)) : attributes
  isempty(attributes) ? p : map(p; Dict(attr=>(x->nothing) for attr in attributes)...)
end

function strip_attributes(p::ACSetTransformation; kw...)
  init = NamedTuple([k=>collect(v) for (k,v) in pairs(components(p))])
  homomorphism(strip_attributes(dom(p); kw...), strip_attributes(codom(p); kw...), initial=init)
end


"""Get all monomorphisms from an acset X to an acset Y
"""
monos(X::ACSet, Y::ACSet) =
  homomorphism(X, strip_attributes(Y); monic = true, type_components=(Name=x->nothing,),)

"""Ask: "does there exists a mono X ↪ Y ?"
"""
exists_mono(X::ACSet,Y::ACSet)::Bool =
  is_homomorphic(X, strip_attributes(Y); monic = true, type_components=(Name=x->nothing,),)


"""
    mca(XX::ACSet, YY::ACSet)

Computes the maximimum common subacsets between XX and YY, i.e., find all a with with |a| maximum possible such that there is a monic span of Acset a₁ ← a → a₂.
"""
function mca(XX::ACSet, YY::ACSet)
  (X,Y) = normalize_order(XX,YY)
  hX = BinaryHeap(Base.By(size,Base.Order.Reverse),[X])
  return mca_help([],hX,Y)
end

"""
    mca_help(mca_list, X_subs, Y)

Helper function to compute the maximimum common subacsets between X and Y. 
Assumes Y is larger than X.
Recursively forms the subacsets of X and checks for monos into Y.

Input: X_subs is a max heap containing the subacsets of X as they are formed, i.e., initially X_subs is the heap with X itself.
Output: mca_list is a vector containing the mca's as they are found, i.e., initially mca_list is an empty vector.
"""
function mca_help(mca_list, X_subs, Y)
  if !isempty(X_subs) && (isempty(mca_list) || size(mca_list[1]) <= size(first(X_subs)))
    curr_X_sub = pop!(X_subs)
    C = acset_schema(curr_X_sub) #X: C → Set
    if exists_mono(curr_X_sub,Y) 
      if length(filter(x -> curr_X_sub==x, mca_list))==0
          push!(mca_list,curr_X_sub)
      end
    else
      indiv_parts = []
      for c in objects(C)
          for p in parts(curr_X_sub,c)
              push!(indiv_parts, NamedTuple([c => [p]]))
          end
      end
      new_X_subs = concatmap(γ -> rm_cascade_subobj(curr_X_sub,γ), indiv_parts)
      for new_sub in new_X_subs
          push!(X_subs,new_sub)
      end
    end
    return mca_help(mca_list, X_subs, Y)
  else # isempty(X_subs) || AlgebraicPetri.SubACSets.size(mca_list[1]) > AlgebraicPetri.SubACSets.size(first(X_subs))
    return mca_list
  end
end

end
