# using Test

using AlgebraicPetri, AlgebraicPetri.SubACSets

m1 = LabelledPetriNet(
  [:X,:Y,:W,:Z,:V],
  :f=>((:X, :W)=>(:Y, :Z))
)

m2 = LabelledPetriNet(
  [:X,:Y,:Z],
  :f=>((:X, :Y)=>:Z),
  :g=>(:Z=>:X)
)

test = mca(m1,m2)

m3 = LabelledPetriNet(
  [:X3,:Y3,:W3,:Z3],
  :f3=>(:X3=>(:W3, :Y3, :Z3))
)

m4 = LabelledPetriNet(
  [:X4,:W4,:Y4,:Z4],
  :f4=>((:W4, :X4, :Y4)=>:Z4)
)

m5 = LabelledPetriNet(
  [:X,:Y],
  :f=>(:X=>:Y)
)
m6 = LabelledPetriNet(
  [:W,:Z],
  :g=>(:W=>:Z)
)

test = mca(m3,m4)

using Catlab.CategoricalAlgebra
import AlgebraicPetri.SubACSets: one_removed_subobs, exists_mono, concatmap, size
import AlgebraicPetri.SubACSets: strip_names

#***
# Copy of existing mca_help function
#***
function mca_help(X::ACSet, Y::ACSet)
    C = acset_schema(X) #X: C → Set
    getsmaller(Z,c) = map(dom, one_removed_subobs(Z, c))
    # enumerate all sub-acsets X' ↪ X of the acset X: C → Set obtained by removing one point from Xc for some c ∈ C
    oneRemovedFromX = concatmap(γ -> getsmaller(X,γ), filter(c -> !isempty(parts(X,c)), objects(C)))
    # keep only those X' ↪ X such that there exists mono X' ↪ Y
    Y_subs          = filter(χ -> exists_mono(χ, Y), oneRemovedFromX)
    map(println, AlgebraicPetri.SubACSets.size.(oneRemovedFromX))
    # either terminate or recurse
    if isempty(Y_subs)
        println("Yay, we are recursing.")
        concatmap(χ -> mca_help(χ, Y), oneRemovedFromX)
    else
      ω = maximum(map(AlgebraicPetri.SubACSets.size,Y_subs))
      filter(y -> AlgebraicPetri.SubACSets.size(y) == ω, Y_subs)
    end
end
  
#***
# Attempts at various syntaxes for making a heap ordered by acset size
using DataStructures
#***
#= 
struct ACSetOrder <: Base.Ordering
end

Lt((x,y) -> isless(AlgebraicPetri.SubACSets.size(a),AlgebraicPetri.SubACSets.size(b)))
lt(::ACSetOrder,a::LabelledPetriNet,b::LabelledPetriNet) = isless(AlgebraicPetri.SubACSets.size(a),AlgebraicPetri.SubACSets.size(b))
eq(::ACSetOrder,a::LabelledPetriNet,b::LabelledPetriNet) = isequal(AlgebraicPetri.SubACSets.size(a),AlgebraicPetri.SubACSets.size(b))
h = BinaryHeap(lt,[m3,m4])
h = BinaryHeap{LabelledPetriNet,ACSetOrder}()

h2 = BinaryHeap{LabelledPetriNet,ACSetOrder}()
h3 = BinaryHeap{LabelledPetriNet,LT}()
h5 = BinaryHeap{ACSet,Base.By(AlgebraicPetri.SubACSets.size,Base.Order.Reverse)}()


struct MyACSetOrder <: Base.Order.Ordering
end

Base.Order.lt(::MyACSetOrder,a::LabelledPetriNet,b::LabelledPetriNet) = isless(AlgebraicPetri.SubACSets.size(a),AlgebraicPetri.SubACSets.size(b))
=#

# This heap seems to work
h4 = BinaryHeap(Base.By(AlgebraicPetri.SubACSets.size,Base.Order.Reverse),[m3])

#***
# A version of mca_help that can use a heap to keep track of the sub-acsets for recursion.
# This does recur, but one_removed_subobs is not shrinking some of the subacsets, so the recursion isn't terminating
#***
function mca_help2(mca_list, X_subs, Y)
    println(length(mca_list)," ",length(X_subs))
    if !isempty(X_subs) && (isempty(mca_list) || AlgebraicPetri.SubACSets.size(mca_list[1]) <= AlgebraicPetri.SubACSets.size(first(X_subs)))
      curr_X_sub = pop!(X_subs)
      C = acset_schema(curr_X_sub) #X: C → Set
      getsmaller(Z,c) = map(dom, one_removed_subobs(Z, c))
      # enumerate all sub-acsets X' ↪ X of the acset X: C → Set obtained by removing one point from Xc for some c ∈ C
      oneRemovedFromX = concatmap(γ -> getsmaller(curr_X_sub,γ), filter(c -> !isempty(parts(curr_X_sub,c)), objects(C)))
      if exists_mono(curr_X_sub,Y) 
        push!(mca_list,curr_X_sub)
      else
        new_X_subs = oneRemovedFromX
        for new_sub in new_X_subs
            push!(X_subs,new_sub)
        end
      end
      return mca_help2(mca_list, X_subs, Y)
    else # isempty(X_subs) || AlgebraicPetri.SubACSets.size(mca_list[1]) > AlgebraicPetri.SubACSets.size(first(X_subs))
      return mca_list
    end
end

#***
# Once the reduction of subacsets is fixed there is still a major inefficiency:
# Not only are many isomorphic subacsets added to the heap, 
# but many duplicates may also be added.
#
# It may be better to check if a subacset is a duplicate before adding
# It should be skipped if it is a duplicate
# If subacsets are isomorphic to an existing subacset, we may want to keep track of their occurance/count, but can skip the search for a mono 
#
# These checks will likely require another structure besides the heap
# Perhaps some hash or binary indicator for overall parts included
# The only different isos need to be checked for monos and only subacsets of the same size need to be checked for isomorphisms
#***

#***
# A copy of the one_removed_subobs function for debugging
# The problem is that when trying to remove an arc/arrow, the double negation puts it back
# Instead of double negation, when it removes a part it should remove any part that refers to it in other tables
#***
C = acset_schema(m3) #X: C → Set
objs_C = objects(C)
prts_I = parts(m3,:I)
prts_O = parts(m3,:O)
  function one_removed_subobs(X::ACSet, c)
    C    = acset_schema(X)
    subs = [filter(j -> j != i, parts(X,c)) for i ∈ parts(X,c)]
    build_sub(sub) = hom(
      # need double negation to remove dangling edges
      ¬¬Subobject(
        X,
        NamedTuple(
          zip(
            objects(C),
            [o == c ? sub : parts(X, o) for o ∈ objects(C)]
          )
        )
      )
    )
    map(build_sub, subs)
  end

  X = m3
  Graph(dom(hom(Subobject(X,NamedTuple(zip(
                   objects(C),
                   [o == :O ? [2,3] : parts(X, o) for o ∈ objects(C)]
                 ))))))



#***
# Kris's original casade -- trim's the subobject of dangling edges
#***
"""Recursively delete anything, e.g. deleting a vertex deletes its edge"""
function cascade_subobj(X::ACSet, sub)
  sub = Dict([k=>Set(v) for (k,v) in pairs(sub)])
  S = acset_schema(X)
  change = true
  while change
    change = false
    for (f,c,d) in homs(S)
      for c_part in sub[c]
        if X[c_part,f] ∉ sub[d]
          change = true
          delete!(sub[c], c_part)
        end
      end
    end
  end
  return Dict([k => collect(v) for (k,v) in pairs(sub)])
end

#***
# Recursively removes subacsets and any dangling edges
#***
"""Recursively delete anything, e.g. deleting a vertex deletes its edge"""
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


# test = rm_cascade_subobj(X,Subobject(X,O=[2]))
test = rm_cascade_subobj(X,Dict(:O=>[2]))
test = rm_cascade_subobj(X,NamedTuple([:O=>[2]]))
test = rm_cascade_subobj(X,NamedTuple([:S=>[2]]))
test = rm_cascade_subobj(X,NamedTuple([:T=>[1]]))


#***
# A version of mca_help2 that uses rm_cascade_subobj instead of one_removed_subobs
#***
function mca_help3(mca_list, X_subs, Y)
    println(length(mca_list)," ",length(X_subs))
    if !isempty(X_subs) && (isempty(mca_list) || AlgebraicPetri.SubACSets.size(mca_list[1]) <= AlgebraicPetri.SubACSets.size(first(X_subs)))
      curr_X_sub = pop!(X_subs)
      C = acset_schema(curr_X_sub) #X: C → Set
      # enumerate all sub-acsets X' ↪ X of the acset X: C → Set obtained by removing one point from Xc for some c ∈ C
      if exists_mono(curr_X_sub,Y) 
        push!(mca_list,curr_X_sub)
      else
        indiv_parts = []
        for c in objects(C)
            for p in parts(curr_X_sub,c)
                push!(indiv_parts, NamedTuple([c => [p]]))
            end
        end
        new_X_subs = concatmap(γ -> rm_cascade_subobj(curr_X_sub,γ), indiv_parts)
        for new_sub in new_X_subs
            # if length(filter(x -> is_isomorphic(new_sub,x),X_subs))>0
                push!(X_subs,new_sub)
            # end
        end
      end
      return mca_help3(mca_list, X_subs, Y)
    else # isempty(X_subs) || AlgebraicPetri.SubACSets.size(mca_list[1]) > AlgebraicPetri.SubACSets.size(first(X_subs))
      return mca_list
    end
end

h4 = BinaryHeap(Base.By(AlgebraicPetri.SubACSets.size,Base.Order.Reverse),[m3])
dollar = mca_help3([],h4,m4)

#***
# mca_help3 appears to work
# But it is currently returning replicates (not really in correspondance to number of possible morphisms either)
# Could consider checking for replicates prior to adding to subacset data structure
# Likewise, could consider grouping as isomorphism classes and keep track of size to reduce comparisons
# Could also consider finding and returning the morphisms as well, in which case the isomorphism classes may be helpful
#***