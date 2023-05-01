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
  [:X,:Y,:W,:Z],
  :f=>(:X=>(:W, :Y, :Z))
)

m4 = LabelledPetriNet(
  [:X,:W,:Y,:Z],
  :f=>((:W, :X, :Y)=>:Z)
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
            if length(filter(x -> is_isomorphic(new_sub,x),X_subs))>0
                push!(X_subs,new_sub)
            end
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

