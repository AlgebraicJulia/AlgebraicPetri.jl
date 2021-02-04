# # [Cathepsin Enzyme Reactions](@id enzyme_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/enzymes/enzyme_reactions.ipynb)

using AlgebraicPetri
using Catlab.Programs
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Distributions

using DifferentialEquations
using Plots

display_uwd(ex, prog="neato") = to_graphviz(ex, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:overlap => "false"), prog=prog)
ode(x::Union{AbstractReactionNet{Distribution, Number},AbstractLabelledReactionNet{Distribution, Number}}, t) = begin
  β = mean.(rates(x))
  ODEProblem(vectorfield(x), concentrations(x), t, β)
end
ode(x, t) = ODEProblem(vectorfield(x), concentrations(x), t, rates(x));
meanRates(rxn, pred) = Dict(tname(rxn,t)=>mean(pred).nt.mean[t] for t in 1:nt(rxn))

# ## Define objects and initial conditions

ob(x) = codom(Open([first(x)], LabelledReactionNet{Distribution,Number}(x), [first(x)])).ob;

# ## Helper functions for generating primitive Petri net operations

inact(in,on::Distribution) = begin
  inact = Symbol(first(in), :inact)
  Open(LabelledReactionNet{Distribution,Number}(unique((in, inact=>0)), ((Symbol(:inact_,first(in)),on),first(in)=>inact)))
end;

bind(in1, in2, on::Distribution, off::Distribution) = begin
  out = Symbol(first(in1),first(in2))
  Open(LabelledReactionNet{Distribution,Number}(unique((in1, in2,out=>0)), ((Symbol(:bind_,first(in1),first(in2)),on),(first(in1),first(in2))=>out),
                                                           ((Symbol(:unbind_,out),off),out=>(first(in1),first(in2)))))
end;

deg(prod1,prod2,on::Distribution) = begin
  in = Symbol(first(prod1),first(prod2))
  prod2str = String(first(prod2))
  degprod2 = Symbol(endswith(prod2str, "inact") ? first(prod2str) : prod2str, :deg)
  Open(LabelledReactionNet{Distribution,Number}(unique((in=>0, prod1,degprod2=>0)), ((Symbol(:deg_,in),on),in=>(first(prod1),degprod2))));
end;

# ## Cathepsin *X* reacting with itself

catX = @relation (X, Xinact, Xdeg) where (X, Xinact, Xdeg, XX, XXinact) begin
  inactX(X, Xinact)
  bindXX(X, XX)
  degXX(XX, X, Xdeg)
  bindXXinact(X, Xinact, XXinact)
  degXXinact(XXinact, X, Xdeg)
end
display_uwd(catX)

# ## Cathepsin *X* reacting with Substrate *Y*

catXsubY = @relation (X, Xinact, Xdeg, Y, Ydeg) where (X, Xinact, Xdeg, Y, XY, Ydeg) begin
  bindXY(X, Y, XY)
  degXY(XY, X, Ydeg)
end
display_uwd(catXsubY)

# ## Cathepsin *X* reacting with Cathepsin *Y*

catXY = @relation (X, Xinact, Xdeg, Y, Yinact, Ydeg) where (X, Xinact, Xdeg, Y, Yinact, Ydeg, XY, XYinact) begin
  bindXY(X, Y, XY)
  degXY(XY, X, Ydeg)
  bindXYinact(X, Yinact, XYinact)
  degXYinact(XYinact, X, Ydeg)
end
display_uwd(catXY)

#######################
# EDIT CONSTANTS HERE #
#######################

K = :K=>33000;
S = :S=>33000;
L = :L=>33000;
Kinact = :Kinact=>0;
Sinact = :Sinact=>0;
Linact = :Linact=>0;
E = :E=>700000;
G = :G=>1300000;

kon = Uniform(6e-5, 6e-2)
koff = Uniform(6e-4, 6e3)
kcat = Uniform(0, 6e4)

rxns = Dict(
  :K => [inact(K, kcat)
         bind(K, K, kon, koff)
         deg(K, K, kcat)
         bind(K, Kinact, kon, koff)
         deg(K, Kinact, kcat)],
  :S => [inact(S, kcat)
         bind(S, S, kon, koff)
         deg(S, S, kcat)
         bind(S, Sinact, kon, koff)
         deg(S, Sinact, kcat)],
  :L => [inact(L, kcat)
         bind(L, L, kon, koff)
         deg(L, L, kcat)
         bind(L, Linact, kon, koff)
         deg(L, Linact, kcat)],
  :KE => [bind(K, E, kon, koff)
          deg(K, E, kcat)],
  :KG => [bind(K, G, kon, koff)
          deg(K, G, kcat)],
  :SE => [bind(S, E, kon, koff)
          deg(S, E, kcat)],
  :SG => [bind(S, G, kon, koff)
          deg(S, G, kcat)],
  :LE => [bind(L, E, kon, koff)
          deg(L, E, kcat)],
  :LG => [bind(L, G, kon, koff)
          deg(L, G, kcat)],
  :KS => [bind(K, S, kon, koff)
          deg(K, S, kcat)
          bind(K, Sinact, kon, koff)
          deg(K, Sinact, kcat)],
  :KL => [bind(K, L, kon, koff)
          deg(K, L, kcat)
          bind(K, Linact, kon, koff)
          deg(K, Linact, kcat)],
  :SK => [bind(S, K, kon, koff)
          deg(S, K, kcat)
          bind(S, Kinact, kon, koff)
          deg(S, Kinact, kcat)],
  :SL => [bind(S, L, kon, koff)
          deg(S, L, kcat)
          bind(S, Linact, kon, koff)
          deg(S, Linact, kcat)],
  :LK => [bind(L, K, kon, koff)
          deg(L, K, kcat)
          bind(L, Kinact, kon, koff)
          deg(L, Kinact, kcat)],
  :LS => [bind(L, S, kon, koff)
          deg(L, S, kcat)
          bind(L, Sinact, kon, koff)
          deg(L, Sinact, kcat)]
);

# ## Helper functions to generate oapply calls

cat(cat) = begin
  catsym = first(cat)
  out = oapply(catX, Dict([:inactX, :bindXX, :degXX, :bindXXinact, :degXXinact] .=> rxns[catsym]), Dict(
    :X=>ob(cat),
    :Xinact=>ob(Symbol(catsym,:inact)=>0),
    :Xdeg=>ob(Symbol(catsym,:deg)=>0),
    :XX=>ob(Symbol(catsym,catsym)=>0),
    :XXinact=>ob(Symbol(catsym,catsym,:inact)=>0)))
  bundle_legs(out, [[1,2,3]])
end

cat_sub(cat1, sub) = begin
  catsym = first(cat1)
  subsym = first(sub)
  catsub = Symbol(catsym, subsym)
  out = oapply(catXsubY, Dict([:bindXY, :degXY] .=> rxns[catsub]), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(catsym,:inact)=>0),
    :Xdeg=>ob(Symbol(catsym,:deg)=>0),
    :Y=>ob(sub),
    :XY=>ob(Symbol(catsym,subsym)=>0),
    :Ydeg=>ob(Symbol(subsym,:deg)=>0)))
  bundle_legs(out, [[1,2,3], [4,5]])
end


cat_cat(cat1, cat2) = begin
  cat1sym = first(cat1)
  cat2sym = first(cat2)
  catcat = Symbol(cat1sym, cat2sym)
  out = oapply(catXY, Dict([:bindXY, :degXY, :bindXYinact, :degXYinact] .=> rxns[catcat]), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(cat1sym,:inact)=>0),
    :Xdeg=>ob(Symbol(cat1sym,:deg)=>0),
    :Y=>ob(cat2),
    :Yinact=>ob(Symbol(cat2sym,:inact)=>0),
    :Ydeg=>ob(Symbol(cat2sym,:deg)=>0),
    :XY=>ob(catcat=>0),
    :XYinact=>ob(Symbol(catcat,:inact)=>0)))
  bundle_legs(out, [[1,2,3], [4,5,6]])
end

functor(x) = oapply(x, Dict(
  :catK=>cat(K),
  :catS=>cat(S),
  :catL=>cat(L),
  :catKcatS=>cat_cat(K,S),
  :catKcatL=>cat_cat(K,L),
  :catScatK=>cat_cat(S,K),
  :catScatL=>cat_cat(S,L),
  :catLcatK=>cat_cat(L,K),
  :catLcatS=>cat_cat(L,S),
  :catKsubE=>cat_sub(K,E),
  :catSsubE=>cat_sub(S,E),
  :catLsubE=>cat_sub(L,E),
  :catKsubG=>cat_sub(K,G),
  :catSsubG=>cat_sub(S,G),
  :catLsubG=>cat_sub(L,G)));

function enzyme_uwd(enzymes::Array{Symbol}, substrates::Array{Symbol})
  rel = RelationDiagram{Symbol}(0)

  chemicals = vcat(substrates, enzymes)

  subs = add_junctions!(rel, length(substrates), variable=substrates)
  enzs = add_junctions!(rel, length(enzymes), variable=enzymes)
  nsubs = length(subs)
  nenzs = length(enzs)

  catx = add_parts!(rel, :Box, nenzs, name=[Symbol("cat$i") for i in enzymes])
  add_parts!(rel, :Port, nenzs, junction=enzs, box=catx)

  for x in 1:nenzs
    for y in 1:nenzs
      if y != x
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])cat$(enzymes[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], enzs[y]], box=catxy)
      end
    end
  end

  for x in 1:nenzs
    for y in 1:nsubs
      catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])sub$(substrates[y])"))
      add_parts!(rel, :Port, 2, junction=[enzs[x], subs[y]], box=catxy)
    end
  end
  add_parts!(rel, :OuterPort, length(chemicals), outer_junction = vcat(subs, enzs))
  rel
end

enzyme_reaction(args...) = enzyme_uwd(args...) |> functor |> apex

# # Manually Defining Models

KSE = @relation (K, S, E) begin
  catK(K)
  catS(S)
  catKcatS(K, S)
  catScatK(S, K)
  catKsubE(K, E)
  catSsubE(S, E)
end
display_uwd(KSE)

# # Defining Models using API

# ## CatK, CatS, Elastin

KSE = enzyme_uwd([:K, :S], [:E])
display_uwd(KSE)
#-

KSE_petri = apex(functor(KSE))
ode_prob = ode(KSE_petri, (0.0, 120.0))
sol = solve(ode_prob)

plot(sol)

#-

plot(sol, lw = 2, ylims = (0, 15000), xlims = (0, 30))

# ## CatK, CatS, CatL, Elastin

KSLE = enzyme_reaction([:K, :S, :L], [:E])
ode_prob = ode(KSLE, (0.0,120.0))
sol = solve(ode_prob)
plot(sol, lw = 1, size = (1066, 600))

#-

plot(sol, ylims = (0, 20000), xlims = (0, 30), lw = 2, size = (1066, 600))

# ## CatK, CatS, CatL, Elastin, Gelatin

KSLEG = enzyme_reaction([:K, :S, :L], [:E, :G])
ode_prob = ode(KSLEG, (0.0,120.0))
sol = solve(ode_prob)
plot(sol, lw = 1, size = (1066, 600))

#-

plot(sol, ylims = (0, 30000), xlims = (0, 50), lw = 2, size = (1066, 600))
