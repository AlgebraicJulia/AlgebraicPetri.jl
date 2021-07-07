module EnzymeReactions

using AlgebraicPetri
using Catlab.Programs
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Distributions

using DifferentialEquations
using Plots

export ob, ode,
       inactivate, bindunbind, degrade,
       enzX, enzXY, enzXsubY,
       enz, enz_enz, enz_sub,
       enzyme_uwd, enzyme_generators

ob(type, x) = codom(Open([first(x)], LabelledReactionNet{type,Number}(x), [first(x)])).ob;
ob(x) = codom(Open([x], LabelledPetriNet(x), [x])).ob;

ode(x::Union{AbstractReactionNet{Distribution, Number},AbstractLabelledReactionNet{Distribution, Number}}, t) = begin
  β = mean.(rates(x))
  ODEProblem(vectorfield(x), concentrations(x), t, β)
end
ode(x, t) = ODEProblem(vectorfield(x), concentrations(x), t, rates(x));

function inactivate(in,on::T) where T
  inact = Symbol(first(in), :_inact)
  Open(LabelledReactionNet{T,Number}(unique((in, inact=>0)), ((Symbol(:inact_,first(in)),on),first(in)=>inact)))
end;

function bindunbind(in1, in2, on::T, off::T) where T
  out = Symbol(first(in1),first(in2))
  Open(LabelledReactionNet{T,Number}(unique((in1, in2,out=>0)), ((Symbol(:bind_,first(in1),first(in2)),on),(first(in1),first(in2))=>out),
                                                                ((Symbol(:unbind_,out),off),out=>(first(in1),first(in2)))))
end;

function degrade(prod1,prod2,on::T) where T
  in = Symbol(first(prod1),first(prod2))
  prod2str = String(first(prod2))
  degprod2 = Symbol(endswith(prod2str, "inact") ? first(prod2str) : prod2str, :_deg)
  Open(LabelledReactionNet{T,Number}(unique((in=>0, prod1,degprod2=>0)), ((Symbol(:deg_,in),on),in=>(first(prod1),degprod2))));
end;

function inactivate(in)
  inact = Symbol(in, :_inact)
  Open(LabelledPetriNet(unique((in, inact)), (Symbol(:inact_,in),in=>inact)))
end;

function bindunbind(in1, in2)
  out = Symbol(in1,in2)
  Open(LabelledPetriNet(unique((in1, in2,out)), (Symbol(:bind_,in1,in2),(in1,in2)=>out),
                                                (Symbol(:unbind_,out),out=>(in1,in2))))
end;

function degrade(prod1,prod2)
  in = Symbol(prod1,prod2)
  prod2str = String(prod2)
  degprod2 = Symbol(endswith(prod2str, "inact") ? first(prod2str) : prod2str, :_deg)
  Open(LabelledPetriNet(unique((in, prod1,degprod2)), (Symbol(:deg_,in),in=>(prod1,degprod2))));
end;

# ## Cathepsin *X* reacting with itself

enzX = @relation (X, Xinact, Xdeg) where (X, Xinact, Xdeg, XX, XXinact) begin
  inactX(X, Xinact)
  bindXX(X, XX)
  degXX(XX, X, Xdeg)
  bindXXinact(X, Xinact, XXinact)
  degXXinact(XXinact, X, Xdeg)
end

# ## Cathepsin *X* reacting with Substrate *Y*

enzXsubY = @relation (X, Xinact, Xdeg, Y, Ydeg) where (X, Xinact, Xdeg, Y, XY, Ydeg) begin
  bindXY(X, Y, XY)
  degXY(XY, X, Ydeg)
end

# ## Cathepsin *X* reacting with Cathepsin *Y*

enzXY = @relation (X, Xinact, Xdeg, Y, Yinact, Ydeg) where (X, Xinact, Xdeg, Y, Yinact, Ydeg, XY, XYinact) begin
  bindXY(X, Y, XY)
  degXY(XY, X, Ydeg)
  bindXYinact(X, Yinact, XYinact)
  degXYinact(XYinact, X, Ydeg)
end

function enz(rxns, cat)
  catsym = first(cat)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzX, Dict([:inactX, :bindXX, :degXX, :bindXXinact, :degXXinact] .=> rxns[catsym]), Dict(
    :X=>ob(obtype, cat),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :XX=>ob(obtype, Symbol(catsym,catsym)=>0),
    :XXinact=>ob(obtype, Symbol(catsym,catsym,:_inact)=>0)))
  bundle_legs(out, [[1,2,3]])
end

function enz_sub(rxns, cat1, sub)
  catsym = first(cat1)
  subsym = first(sub)
  catsub = Symbol(catsym, subsym)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXsubY, Dict([:bindXY, :degXY] .=> rxns[catsub]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :Y=>ob(obtype, sub),
    :XY=>ob(obtype, Symbol(catsym,subsym)=>0),
    :Ydeg=>ob(obtype, Symbol(subsym,:_deg)=>0)))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function enz_enz(rxns, cat1, cat2)
  cat1sym = first(cat1)
  cat2sym = first(cat2)
  catcat = Symbol(cat1sym, cat2sym)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXY, Dict([:bindXY, :degXY, :bindXYinact, :degXYinact] .=> rxns[catcat]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(cat1sym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(cat1sym,:_deg)=>0),
    :Y=>ob(obtype, cat2),
    :Yinact=>ob(obtype, Symbol(cat2sym,:_inact)=>0),
    :Ydeg=>ob(obtype, Symbol(cat2sym,:_deg)=>0),
    :XY=>ob(obtype, catcat=>0),
    :XYinact=>ob(obtype, Symbol(catcat,:_inact)=>0)))
  bundle_legs(out, [[1,2,3], [4,5,6]])
end

function enz(cat::Symbol)
  catsym = cat
  out = oapply(enzX, Dict(:inactX=>inactivate(cat), :bindXX=>bindunbind(cat, cat), :degXX=>degrade(cat, cat),
                          :bindXXinact=>bindunbind(cat, Symbol(cat,:_inact)),
                          :degXXinact=>degrade(cat, Symbol(cat, :_inact))), Dict(
    :X=>ob(cat),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :XX=>ob(Symbol(catsym,catsym)),
    :XXinact=>ob(Symbol(catsym,catsym,:_inact))))
  bundle_legs(out, [[1,2,3]])
end

function enz_sub(cat1::Symbol, sub::Symbol)
  catsym = cat1
  subsym = sub
  catsub = Symbol(catsym, subsym)
  out = oapply(enzXsubY, Dict(:bindXY=>bindunbind(cat1, sub), :degXY=>degrade(cat1, sub)), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :Y=>ob(sub),
    :XY=>ob(Symbol(catsym,subsym)),
    :Ydeg=>ob(Symbol(subsym,:_deg))))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function enz_enz(cat1::Symbol, cat2::Symbol)
  cat1sym = cat1
  cat2sym = cat2
  catcat = Symbol(cat1sym, cat2sym)
  out = oapply(enzXY, Dict(:bindXY=>bindunbind(cat1, cat2), :degXY=>degrade(cat1, cat2), :bindXYinact=>bindunbind(cat1, Symbol(cat2, :_inact)), :degXYinact=>degrade(cat1, Symbol(cat2, :_inact))), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(cat1sym,:_inact)),
    :Xdeg=>ob(Symbol(cat1sym,:_deg)),
    :Y=>ob(cat2),
    :Yinact=>ob(Symbol(cat2sym,:_inact)),
    :Ydeg=>ob(Symbol(cat2sym,:_deg)),
    :XY=>ob(catcat),
    :XYinact=>ob(Symbol(catcat,:_inact))))
  bundle_legs(out, [[1,2,3], [4,5,6]])
end

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

function enzyme_generators(enzymes::Array{Symbol}, substrates::Array{Symbol})
  gens = Dict{Symbol, Any}()
  for e1 in enzymes
    for e2 in enzymes
      if e1 == e2
        gens[Symbol(:cat, e1)] = enz(e1)
      else
        gens[Symbol(:cat, e1, :cat, e2)] = enz_enz(e1, e2)
      end
    end
    for s in substrates
      gens[Symbol(:cat, e1, :sub, s)] = enz_sub(e1, s)
    end
  end
  gens
end
end