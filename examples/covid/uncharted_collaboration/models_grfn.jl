using Catlab
using Catlab.Theories
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.Graphics.Graphviz: run_graphviz
using JSON

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
display_wd(ex) = to_graphviz(ex, orientation=TopToBottom);
save_wd(ex, fname::AbstractString, format="svg") = begin
    g = display_wd(ex)
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, g, format=format)
    end
end
save_wd_json(ex, fname::AbstractString) = begin
    g = display_wd(ex)
    open(string(fname, ".json"), "w") do io
        run_graphviz(io, g, format="json0")
    end
end
save_json(c, fname::AbstractString) = begin
  open(string(fname, ".json"), "w") do io
    println(io, json(c.tables))
  end
end


@present GrFN_SIR(FreeBiproductCategory) begin
  F::Ob

  γ::Hom(munit(), F)
  β::Hom(munit(), F)
  dt::Hom(munit(), F)

  recovery::Hom(otimes(F,F,F), F)
  infection::Hom(otimes(F,F,F,F,F), F)
  update_S::Hom(otimes(F,F), F)
  update_I::Hom(otimes(F,F,F), F)
  update_R::Hom(otimes(F,F), F)
end;

sir = @program GrFN_SIR (s::F, i::F, r::F) begin
  dt = dt()
  β = β()
  γ = γ()
  infected = infection(s, i, r, β, dt)
  recovered = recovery(i, γ, dt)
  return update_S(s, infected), update_I(i, infected, recovered), update_R(r, recovered)
end;
display_wd(sir)
save_wd(sir, "sir_grfn", "png")
save_wd(sir, "sir_grfn")
save_wd_json(sir, "sir_grfn")

@present GrFN_SEIR(FreeBiproductCategory) begin
    F::Ob

    S::Hom(munit(), F)
    E::Hom(munit(), F)
    I::Hom(munit(), F)
    R::Hom(munit(), F)
    t_rec::Hom(munit(), F)
    rt::Hom(munit(), F)
    t_inf::Hom(munit(), F)
    n::Hom(otimes(F,F,F,F),F)
    mu_e::Hom(otimes(F,F),F)
    mu_s::Hom(otimes(F,F),F)
    mu_i::Hom(otimes(F,F),F)
    drdt::Hom(otimes(F,F),F)
    dedt::Hom(otimes(F,F,F,F,F,F,F,F),F)
    dsdt::Hom(otimes(F,F,F,F,F,F),F)
    didt::Hom(otimes(F,F,F,F,F,F),F)
end

seir = @program GrFN_SEIR (s::F,e::F,i::F,r::F) begin
    t_rec = t_rec()
    rt = rt()
    t_inf = t_inf()
    n = n(s,e,i,r)
    mu_e = mu_e(rt, t_inf)
    mu_s = mu_s(rt, t_inf)
    mu_i = mu_i(rt, t_inf)
    drdt = drdt(t_rec, i)
    dedt = dedt(rt,t_rec,t_inf,n,mu_e,s,e,i)
    dsdt = dsdt(rt,t_inf,mu_s,n,s,i)
    didt = didt(t_rec,t_inf,mu_i,n,s,e)
    return dsdt,dedt,didt,drdt
end
display_wd(seir)
save_wd(seir, "seir_grfn", "png")
save_wd(seir, "seir_grfn")
save_wd_json(seir, "seir_grfn")

@present GrFN_SEIRP(FreeBiproductCategory) begin
    F::Ob

    S::Hom(munit(),F)
    E::Hom(munit(),F)
    I::Hom(munit(),F)
    R::Hom(munit(),F)

    t_a::Hom(munit(),F)
    et_a::Hom(munit(),F)
    r_a::Hom(munit(),F)
    t_b::Hom(munit(),F)
    r_b::Hom(munit(),F)

    inc_inf_a::Hom(otimes(F,F),F)
    inc_inf_b::Hom(otimes(F,F,F),F)
    inc_exp::Hom(otimes(F,F,F),F)
    inc_rec::Hom(otimes(F,F),F)

    dsdt::Hom(otimes(F,F),F)
    dedt::Hom(otimes(F,F),F)
    didt::Hom(otimes(F,F),F)
    drdt::Hom(F,F)
end

seirp = @program GrFN_SEIRP (s::F,ea::F,ia::F,ib::F) begin
  t_a = t_a()
  et_a = et_a()
  r_a = r_a()
  t_b = t_b()
  r_b = r_b()

  exp_a = inc_exp(s, ia, t_a)
  inf_a = inc_inf_a(ea, et_a)
  rec_a = inc_rec(ia, r_a)

  inf_b = inc_inf_b(s, ib, t_b)
  rec_b = inc_rec(ib, r_b)

  dsadt = dsdt(exp_a, inf_b)
  deadt = dedt(exp_a, inf_a)
  diadt = didt(inf_a, rec_a)
  dradt = drdt(rec_a)
  dibdt = didt(inf_b, rec_b)
  drbdt = drdt(rec_b)

  return dsadt, deadt, diadt, dradt, dibdt, drbdt
end
display_wd(seirp)
save_wd(seirp, "seirp_grfn", "png")
save_wd(seirp, "seirp_grfn")
save_wd_json(seirp, "seirp_grfn")

@present GrFN_CHIME(FreeBiproductCategory) begin
    F::Ob
    B::Ob

    S::Hom(munit(),F)
    I::Hom(munit(),F)
    R::Hom(munit(),F)

    γ::Hom(munit(),F)
    n::Hom(munit(),F)
    doubling_time::Hom(munit(),F)
    relative_contact_rate::Hom(munit(),F)
    growth_rate::Hom(munit(),F)

    inv_contact_rate::Hom(F,F)
    double_to_growth_rate::Hom(F,F)

    cond_0_0::Hom(F,B)

    growth_rate_cond::Hom(otimes(F,F,B),F)
    update_growth_rate::Hom(otimes(F,F),F)
    beta::Hom(otimes(F,F,F),F)
    scale::Hom(otimes(F,F,F,F),F)
    apply_scale::Hom(otimes(F,F),F)

    update_S::Hom(otimes(F,F,F), F)
    update_I::Hom(otimes(F,F,F,F), F)
    update_R::Hom(otimes(F,F,F), F)
end

get_growth_rate = @program GrFN_CHIME (doublingtime::F) begin
  growth_true = growth_rate()
  grouth_false = double_to_growth_rate(doublingtime)
  cond = cond_0_0(doublingtime)
  return growth_rate_cond(growth_true, grouth_false, cond)
end;
save_wd(get_growth_rate, "get_growth_rate_grfn", "png")
save_wd(get_growth_rate, "get_growth_rate_grfn")
save_wd_json(get_growth_rate, "get_growth_rate_grfn")
add_definition!(GrFN_CHIME, :get_growth_rate, to_hom_expr(FreeBiproductCategory, get_growth_rate))

get_beta = @program GrFN_CHIME (s::F, gamma::F, doublingtime::F, relativecontactrate::F) begin
  growth_rate = get_growth_rate(doublingtime)
  inv_contact_rate = inv_contact_rate(relativecontactrate)
  updated_growth_rate = update_growth_rate(gamma, growth_rate)
  return beta(s, updated_growth_rate, inv_contact_rate)
end
save_wd(get_beta, "get_beta_grfn", "png")
save_wd(get_beta, "get_beta_grfn")
save_wd_json(get_beta, "get_beta_grfn")
add_definition!(GrFN_CHIME, :get_beta, to_hom_expr(FreeBiproductCategory, get_beta))

chime = @program GrFN_CHIME (s::F, i::F, r::F) begin
  n = n()
  γ = γ()
  doubling_time = doubling_time()
  relative_contact_rate = relative_contact_rate()
  β = get_beta(s, γ, doubling_time, relative_contact_rate)

  new_s = update_S(s, i, β)
  new_i = update_I(s, i, β, γ)
  new_r = update_R(i, r, γ)

  scale = scale(new_s, new_i, new_r, n)

  return apply_scale(new_s, scale), apply_scale(new_i, scale), apply_scale(new_r, scale)
end
display_wd(chime)
save_wd(chime, "chime_grfn_high_abstract", "png")
save_wd(chime, "chime_grfn_high_abstract")
save_wd_json(chime, "chime_grfn_high_abstract")

chime_med_abstract = expand(GrFN_CHIME, chime)
display_wd(chime_med_abstract)
save_wd(chime_med_abstract, "chime_grfn_med_abstract", "png")
save_wd(chime_med_abstract, "chime_grfn_med_abstract")
save_wd_json(chime_med_abstract, "chime_grfn_med_abstract")

chime_low_abstract = expand(GrFN_CHIME, chime_med_abstract)
display_wd(chime_low_abstract)
save_wd(chime_low_abstract, "chime_grfn_low_abstract", "png")
save_wd(chime_low_abstract, "chime_grfn_low_abstract")
save_wd_json(chime_low_abstract, "chime_grfn_low_abstract")