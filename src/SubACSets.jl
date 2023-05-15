module SubACSets
export mca

using Catlab.CategoricalAlgebra
using DataStructures

"""
    rm_cascade_subobj(X::ACSet, rm_subs)

Deletes parts from an ACSet in cascading fashion, e.g. deleting a vertex deletes its edges
rm_subs is a NamedTuple or Dict of parts to be removed.
"""
function rm_cascade_subobj(X::ACSet, rm_subs)
  # HACK: Remove this once Catlab makes cascading delete default
  # https://github.com/AlgebraicJulia/Catlab.jl/pull/605 (old PR)
  S = acset_schema(X)
  subs = Dict([k => Set(parts(X, k)) for k ∈ objects(S)])
  rm_subs = Dict([k => Set(v) for (k, v) ∈ pairs(rm_subs)])
  while !isempty(rm_subs)
    curr_c = first(rm_subs)[1]
    if isempty(rm_subs[curr_c])
      delete!(rm_subs, curr_c)
    else
      curr_part = pop!(rm_subs[curr_c])
      if curr_part ∈ subs[curr_c]
        delete!(subs[curr_c], curr_part)
        for (f, c, d) ∈ homs(S)
          if d == curr_c && c ∈ keys(subs)
            for test_part ∈ subs[c]
              if X[test_part, f] == curr_part
                if c ∈ keys(rm_subs)
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
        delete!(rm_subs, curr_c)
      end
    end
  end
  dom(hom(Subobject(X, NamedTuple(k => collect(v) for (k, v) ∈ subs))))
end

"""
Defintion: let 𝐺: C → 𝐒et be a C-set, we define the _size_ of 𝐺 to be ∑_{c ∈ C}
|𝐺c|.  For example, under this definition, the size of:
  * a graph G is |GE| + |GV| (num edges + num vertices)
  * a Petri net P is |PT| + |PS| + |PI| + |PO| (num transitions + num species +
    num input arcs + num output arcs).
"""
size(X::ACSet) = foldl(+, [length(parts(X, oₛ)) for oₛ ∈ objects(acset_schema(X))])

function strip_attributes(p::ACSet)
  attributes = attrtypes(acset_schema(p))
  isempty(attributes) ? p : map(p; Dict(attr => (x -> nothing) for attr ∈ attributes)...)
end

# Ask: "does there exists a mono X ↪ Y ?"
exists_mono(X::ACSet, Y::ACSet)::Bool =
  is_homomorphic(X, strip_attributes(Y); monic=true, type_components=(Name=x -> nothing,))

"""
    mca(XX::ACSet, YY::ACSet)

Computes the maximimum common subacsets between XX and YY, i.e., find all a with with |a| maximum possible such that there is a monic span of Acset a₁ ← a → a₂.
"""
function mca(XX::ACSet, YY::ACSet)
  # (X, Y) = size(XX) ≤ size(YY) ? (XX, YY) : (YY, XX) # normalize order
  mca([XX, YY])
end

function mca_help(X::T, Y::Union{T,Vector{T}}; f_reverse = false) where T <: ACSet
  X_subs = BinaryHeap(Base.By(size, Base.Order.Reverse), [X])
  mca_list = Set{T}()
  if typeof(Y)==Vector{T} f_reverse = false end
  while_cond = f_reverse ? size(Y) <= size(first(X_subs)) :
                (isempty(mca_list) || size(first(mca_list)) <= size(first(X_subs)))
  while !isempty(X_subs) && while_cond
    curr_X_sub = pop!(X_subs)
    C = acset_schema(curr_X_sub) #X: C → Set
    match_cond = f_reverse ? is_isomorphic(strip_attributes(curr_X_sub),strip_attributes(Y)) : 
                  all([exists_mono(curr_X_sub, x) for x in Y])
    if match_cond
      # if all([!is_isomorphic(curr_X_sub,tmp) for tmp in mca_list])
        push!(mca_list, curr_X_sub)
      # end
    else
      indiv_parts = []
      for c ∈ objects(C)
        for p ∈ parts(curr_X_sub, c)
          push!(indiv_parts, NamedTuple([c => [p]]))
        end
      end
      new_X_subs = mapreduce(γ -> rm_cascade_subobj(curr_X_sub, γ), vcat, indiv_parts; init=[])
      for new_sub ∈ new_X_subs
        push!(X_subs, new_sub)
      end
    end
    while_cond = f_reverse ? size(Y) <= size(first(X_subs)) :
                  (isempty(mca_list) || size(first(mca_list)) <= size(first(X_subs)))
  end
  return collect(mca_list), X_subs
end

function mca(X::Vector{T}) where T <: ACSet
  acset_order = sortperm(size.(X))
  mca_match1, _ = mca_help(X[acset_order[1]],X[acset_order[2:end]])
  X_mca = mca_match1[1]

  @assert all([is_isomorphic(strip_attributes(X_mca),strip_attributes(match)) for match in mca_match1])
  mca_morphs = [Vector{ACSetTransformation}() for _ in 1:length(X)]
  C = acset_schema(X[acset_order[1]])
  mca_morphs[acset_order[1]] = [ACSetTransformation(strip_attributes(X_mca),X[acset_order[1]]; 
                                  Dict([k => parts(match, k) for k ∈ objects(C)])...) for match in mca_match1]
  for (jj, curr_X) in enumerate(X[acset_order[2:end]])
    curr_X_matches, _ = mca_help(curr_X,X_mca;f_reverse=true)
    mca_morphs[acset_order[jj+1]] = [ACSetTransformation(strip_attributes(X_mca),curr_X; 
                                      Dict([k => parts(match, k) for k ∈ objects(C)])...) for match in curr_X_matches]
  end
  
  return strip_attributes(X_mca), mca_morphs
end

function mca_spans(morphs::Vector{Vector{ACSetTransformation}})
  return [Span(v...) for v in Iterators.product(morphs...)]
end

#=
function mca_pushouts(morphs::Vector{Vector{ACSetTransformation}})
  return [pushout(v...) for v in Iterators.product(morphs...)]
end
=#

end
