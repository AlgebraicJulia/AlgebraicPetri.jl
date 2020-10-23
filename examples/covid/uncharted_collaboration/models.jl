using AlgebraicPetri

using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.CategoricalAlgebra.Squares
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.Graphics.Graphviz: run_graphviz

function expand(p::Presentation, wd::WiringDiagram)
  eqs = Dict(pair for pair in p.equations if head(first(pair)) == :generator)
  substitute(functor(wd, identity, box->begin
    val = p[box.value]
    if val in keys(eqs)
      val = eqs[val]
    end
    to_wiring_diagram(val)
  end))
end
expand(p::Presentation, ex) = expand(p, to_wiring_diagram(ex))
display_wd(ex) = to_graphviz(ex, orientation=LeftToRight, labels=true);

# Define some helper types where transition rates are
# real numbers and populations are natural numbers

# Define helper functions for defining the two types of
# reactions in an epidemiology Model. Either a state
# spontaneously changes, or one state causes another to change

ob(x::Symbol) = codom(Open([x], LabelledPetriNet(x), [x]));
spontaneous_petri(transition::Symbol, s::Symbol, t::Symbol) =
  Open([s], LabelledPetriNet(collect(Set([s,t])), (transition,s=>t)), [t])
exposure_petri(transition::Symbol, s::Symbol, e::Symbol, t::Symbol) =
  Open([s, e], LabelledPetriNet(collect(Set([s,e,t])), (transition,((s,e)=>(t,e)))), [t])

# Extend the infectious disease presentation to handle the more
# complex version of SEIRD that the COEXIST model uses with
# asymptomatic infection, multiple stages of infection, and
# multiple stages of recovery

@present Epidemiology(FreeBiproductCategory) begin
  (S, E, I, R)::Ob

  infection::Hom(S⊗I,I)
  exposure::Hom(S⊗I,E)
  illness::Hom(E,I)
  recovery::Hom(I,R)
end;
exp_ill = @program Epidemiology (s::S, i::I) begin
  e = exposure(s, i)
  i2 = illness(e)
  return [i, i2]
end;
exp_ill = to_hom_expr(FreeBiproductCategory, exp_ill);
add_definition!(Epidemiology, :exp_ill, exp_ill)

# Define a functor from the presentation to the building block
# Petri nets that define these operations

epi_generators(n) = begin
  s,e,i,r = Symbol(:S, n), Symbol(:E, n), Symbol(:I, n), Symbol(:R, n);
  Dict(
    Epidemiology[:S]=>ob(s),
    Epidemiology[:E]=>ob(e),
    Epidemiology[:I]=>ob(i),
    Epidemiology[:R]=>ob(r),

    Epidemiology[:infection]=>exposure_petri(:inf, s, i, i),
    Epidemiology[:exposure]=>exposure_petri(:exp, s, i, e),
    Epidemiology[:illness]=>spontaneous_petri(:ill, e, i),
    Epidemiology[:recovery]=>spontaneous_petri(:rec, i, r),
    Epidemiology[:exp_ill]=>(exposure_petri(:exp, s, i, e) ⋅ spontaneous_petri(:ill, e, i))
  );
end;

F_epi(ex, n) = functor((OpenLabelledPetriNetOb, OpenLabelledPetriNet), ex, generators=epi_generators(n));

sir = Epidemiology[:infection] ⋅ Epidemiology[:recovery];
display_wd(sir)
Graph(F_epi(sir, 1))

# Define SEIR by the primitives

seir = @program Epidemiology (s::S, i::I) begin
  e = exposure(s, i)
  i2 = illness(e)
  return recovery([i, i2])
end;
seir = to_hom_expr(FreeBiproductCategory, seir);
display_wd(seir)
Graph(F_epi(seir, 1))

# Define SEIR like sir, but replace infection => exp_ill

seir = Epidemiology[:exp_ill] ⋅ Epidemiology[:recovery];
display_wd(seir)
display_wd(expand(Epidemiology, seir))
Graph(F_epi(seir, 1))

@present EpidemiologyMasks <: Epidemiology begin
  (Ms, Me, Mi)::Ob

  mask_s::Hom(S,Ms)
  mask_e::Hom(E,Me)
  mask_i::Hom(I,Mi)
  unmask_ms::Hom(Ms,S)
  unmask_me::Hom(Me,E)
  unmask_mi::Hom(Mi,I)
  exposure_s_mi::Hom(S⊗Mi,E)
  exposure_ms_i::Hom(Ms⊗I,Me)
  exposure_ms_mi::Hom(Ms⊗Mi,Me)
  illness_me::Hom(Me,Mi)
  recovery_mi::Hom(Mi,R)
end;

masks_generators(n) = begin
  s,e,i,r = Symbol(:S, n), Symbol(:E, n), Symbol(:I, n), Symbol(:R, n);
  ms,me,mi = Symbol(:Ms, n), Symbol(:Me, n), Symbol(:Mi, n);
  Dict(
    EpidemiologyMasks[:Ms]=>ob(ms),
    EpidemiologyMasks[:Me]=>ob(me),
    EpidemiologyMasks[:Mi]=>ob(mi),

    EpidemiologyMasks[:mask_s]=>spontaneous_petri(:mask_s, s, ms),
    EpidemiologyMasks[:mask_e]=>spontaneous_petri(:mask_e, e, me),
    EpidemiologyMasks[:mask_i]=>spontaneous_petri(:mask_i, i, mi),
    EpidemiologyMasks[:unmask_ms]=>spontaneous_petri(:unmask_ms, ms, s),
    EpidemiologyMasks[:unmask_me]=>spontaneous_petri(:unmask_me, me, e),
    EpidemiologyMasks[:unmask_mi]=>spontaneous_petri(:unmask_mi, mi, i),
    EpidemiologyMasks[:exposure_s_mi]=>exposure_petri(:exp_s_mi, s, mi, e),
    EpidemiologyMasks[:exposure_ms_i]=>exposure_petri(:exp_ms_i, ms, i, me),
    EpidemiologyMasks[:exposure_ms_mi]=>exposure_petri(:exp_ms_mi, ms, mi, me),
    EpidemiologyMasks[:illness_me]=>spontaneous_petri(:ill_me, me, mi),
    EpidemiologyMasks[:recovery_mi]=>spontaneous_petri(:rec_mi, mi, r)
  );
end;

F_masks(ex, n) = functor((OpenLabelledPetriNetOb, OpenLabelledPetriNet), ex, generators=merge(epi_generators(n), masks_generators(n)));

seirm_no_unmask = @program EpidemiologyMasks (s::S, i::I) begin
  e = exposure(s,i)
  ms = mask_s(s)
  me = mask_e(e)
  mi = mask_i(i)
  e2 = exposure_s_mi(s, mi)
  e_all = [e, e2]
  me2 = exposure_ms_i(ms, i)
  me3 = exposure_ms_mi(ms, mi)
  me_all = [me, me2, me3]
  i2 = illness(e_all)
  i_all = [i, i2]
  mi2 = illness_me(me_all)
  mi_all = [mi, mi2]
  return [recovery(i_all), recovery_mi(mi_all)]
end
seirm_no_unmask = to_hom_expr(FreeBiproductCategory, seirm_no_unmask);
display_wd(seirm_no_unmask)
Graph(F_masks(seirm_no_unmask, 1))

mask_exp_ill = @program EpidemiologyMasks (s::S, i::I) begin
  e = exposure(s,i)
  ms = mask_s(s)
  me = mask_e(e)
  mi = mask_i(i)
  e2 = exposure_s_mi(s, mi)
  me2 = exposure_ms_i(ms, i)
  me3 = exposure_ms_mi(ms, mi)
  me_all = [me, me2, me3]
  s_all = [s, unmask_ms(ms)]
  e_all = [e, e2, unmask_me(me_all)]
  i2 = unmask_mi(mi)
  i_all = [i, i2, illness(e_all)]
  mi_all = [mi, illness_me(me_all)]

  return i_all, mi_all
end
mask_exp_ill = to_hom_expr(FreeBiproductCategory, mask_exp_ill);

mask_rec = (EpidemiologyMasks[:recovery] ⊗ EpidemiologyMasks[:recovery_mi]) ⋅ ∇(EpidemiologyMasks[:R])

@present EpidemiologyMasksBoxed <: EpidemiologyMasks begin end
add_definition!(EpidemiologyMasksBoxed, :mask_exp_ill, mask_exp_ill)
add_definition!(EpidemiologyMasksBoxed, :mask_rec, mask_rec)

masksboxed_generators(n) = Dict(
  EpidemiologyMasksBoxed[:mask_exp_ill]=>F_masks(mask_exp_ill, n),
  EpidemiologyMasksBoxed[:mask_rec]=>F_masks(mask_rec, n)
)
F_boxed(ex, n) = functor((OpenLabelledPetriNetOb, OpenLabelledPetriNet), ex, generators=merge(epi_generators(n), masks_generators(n), masksboxed_generators(n)));

seirm = EpidemiologyMasksBoxed[:mask_exp_ill] ⋅ EpidemiologyMasksBoxed[:mask_rec]
display_wd(seirm)
display_wd(expand(EpidemiologyMasksBoxed, seirm))
Graph(F_boxed(seirm, 1))

# Building transformations by hand:

t = id(EpidemiologyMasksBoxed[:S] ⊗ EpidemiologyMasksBoxed[:I])
b = id(EpidemiologyMasksBoxed[:I])
l = EpidemiologyMasksBoxed[:infection]
r = EpidemiologyMasksBoxed[:exp_ill]

add_exposure = SquareDiagram(t, b, l, r)

t2 = id(EpidemiologyMasksBoxed[:S] ⊗ EpidemiologyMasksBoxed[:I])
b2 = Hom(:add_mask_inf, EpidemiologyMasksBoxed[:I], (EpidemiologyMasksBoxed[:I] ⊗ EpidemiologyMasksBoxed[:Mi]))
l2 = r
r2 = EpidemiologyMasksBoxed[:mask_exp_ill]

add_mask_exp = SquareDiagram(t2, b2, l2, r2)

t3 = b
b3 = id(EpidemiologyMasksBoxed[:R])
l3 = EpidemiologyMasksBoxed[:recovery]
r3 = EpidemiologyMasksBoxed[:recovery]

id_rec = SquareDiagram(t3, b3, l3, r3)

t4 = b2
b4 = id(EpidemiologyMasksBoxed[:R])
l4 = r3
r4 = EpidemiologyMasksBoxed[:mask_rec]

add_mask_rec = SquareDiagram(t4, b4, l4, r4)

infection2mask_exp_ill = composeH(add_exposure, add_mask_exp)
rec2mask_rec = composeH(id_rec, add_mask_rec)

sir_seirmH = composeV(infection2mask_exp_ill, rec2mask_rec)

sir2seir = composeV(add_exposure, id_rec)
seir2seirm = composeV(add_mask_exp, add_mask_rec)

sir_seirmV = composeH(sir2seir, seir2seirm)

# Using a presentation for double category transformations

@present EpiTransforms(FreeSymmetricMonoidalDoubleCategory) begin
  (SI, I, MI, R)::Ob

  add_mask_inf::HomH(I,MI)
  infection::HomV(SI,I)
  exp_ill::HomV(SI,I)
  mask_exp_ill::HomV(SI,MI)
  rec::HomV(I,R)
  mask_rec::HomV(MI,R)

  inf2exp_ill::Hom2(idH(SI), idH(I), infection, exp_ill)
  exp_ill2mask_exp_ill::Hom2(idH(SI), add_mask_inf, exp_ill, mask_exp_ill)
  rec2mask_rec::Hom2(add_mask_inf, idH(R), rec, mask_rec)

  sir2seir := inf2exp_ill ⋅ id2H(rec)
  seir2seirm := exp_ill2mask_exp_ill ⋅ rec2mask_rec
  inf2mask_exp_ill := inf2exp_ill ⋆ exp_ill2mask_exp_ill
  recid2mask_rec := id2H(rec) ⋅ rec2mask_rec
end;

epitransforms_obs_homs = Dict(
  EpiTransforms[:SI]=>EpidemiologyMasksBoxed[:S] ⊗ EpidemiologyMasksBoxed[:I],
  EpiTransforms[:I]=>EpidemiologyMasksBoxed[:I],
  EpiTransforms[:MI]=>EpidemiologyMasksBoxed[:I] ⊗ EpidemiologyMasksBoxed[:Mi],
  EpiTransforms[:R]=>EpidemiologyMasksBoxed[:R],
  EpiTransforms[:infection]=>EpidemiologyMasksBoxed[:infection],
  EpiTransforms[:exp_ill]=>EpidemiologyMasksBoxed[:exp_ill],
  EpiTransforms[:mask_exp_ill]=>EpidemiologyMasksBoxed[:mask_exp_ill],
  EpiTransforms[:rec]=>EpidemiologyMasksBoxed[:recovery],
  EpiTransforms[:mask_rec]=>EpidemiologyMasksBoxed[:mask_rec],
  EpiTransforms[:add_mask_inf]=>Hom(:add_mask_inf, EpidemiologyMasksBoxed[:I], (EpidemiologyMasksBoxed[:I] ⊗ EpidemiologyMasksBoxed[:Mi]))
)

epitransforms_squares = begin
  bd(sym) = epitransforms_obs_homs[EpiTransforms[sym]]
  Dict(
    EpiTransforms[:inf2exp_ill]=>SquareDiagram(id(bd(:SI)), id(bd(:I)), bd(:infection), bd(:exp_ill)),
    EpiTransforms[:exp_ill2mask_exp_ill]=>SquareDiagram(id(bd(:SI)), bd(:add_mask_inf), bd(:exp_ill), bd(:mask_exp_ill)),
    EpiTransforms[:rec2mask_rec]=>SquareDiagram(bd(:add_mask_inf), id(bd(:R)), bd(:rec), bd(:mask_rec)),
  )
end

epitransforms_definitions = begin
  bd(sym) = epitransforms_obs_homs[EpiTransforms[sym]]
  sq(sym) = epitransforms_squares[EpiTransforms[sym]]
  Dict(
    EpiTransforms[:sir2seir]=>composeV(sq(:inf2exp_ill), SquareDiagram(id(bd(:I)), id(bd(:R)), bd(:rec), bd(:rec))),
    EpiTransforms[:seir2seirm]=>composeV(sq(:exp_ill2mask_exp_ill), sq(:rec2mask_rec)),
    EpiTransforms[:inf2mask_exp_ill]=>composeH(sq(:inf2exp_ill), sq(:exp_ill2mask_exp_ill)),
    EpiTransforms[:recid2mask_rec]=>composeH(SquareDiagram(id(bd(:I)), id(bd(:R)), bd(:rec), bd(:rec)), sq(:rec2mask_rec))
  )
end

F_transform(ex) = functor((FreeBiproductCategory.Ob, FreeBiproductCategory.Hom, FreeBiproductCategory.Hom, SquareDiagram), ex, generators=merge(epitransforms_obs_homs, epitransforms_squares, epitransforms_definitions));

et(n) = EpiTransforms[n]

seir2seirmH = F_transform(et(:sir2seir) ⋆ et(:seir2seirm))

seir2seirmV = F_transform(et(:inf2mask_exp_ill) ⋅ et(:recid2mask_rec))