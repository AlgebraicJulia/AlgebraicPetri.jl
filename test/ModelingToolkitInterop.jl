module ModelingToolkitInterop
  using Test
  using ModelingToolkit
  using ModelingToolkit: unwrap, term
  using AlgebraicPetri
  using AlgebraicPetri.Epidemiology
  using AlgebraicPetri.BilayerNetworks

  using Catlab
  using Catlab.CategoricalAlgebra
  import Catlab.CategoricalAlgebra: migrate!
  using Catlab.WiringDiagrams
  using Catlab.Programs.RelationalPrograms

  # Define SIR Model
  sir = @relation (s, i, r) begin
    infection(s, i)
    recovery(i, r)
  end

  # Convert to Epidemiology petri net
  psir = apex(oapply_epi(sir))

  # Create empty bilayer network
  bnsir = LabelledBilayerNetwork()
  # migrate petri model to bilayer network
  migrate!(bnsir, psir)

  @variables t S(t) I(t) R(t)
  @parameters inf rec
  D = Differential(t)
  ϕ1 = inf * S * I
  ϕ2 = rec * I
  plus(x, y) = term(+, unwrap(x), unwrap(y); type = Real)
  minus(x, y) = term(-, unwrap(x), unwrap(y); type = Real)
  minus(x) = term(-, unwrap(x); type = Real)
  eqs = [D(S) ~ -ϕ1, D(I) ~ ϕ1 + ϕ1 - (ϕ1 + ϕ2), D(R) ~ ϕ2]
  eqs2 = [D(S) ~ minus(ϕ1), D(I) ~ minus(plus(ϕ1, ϕ1), plus(ϕ1, ϕ2)), D(R) ~ ϕ2]
  bilayer_example = ODESystem(eqs2, t, name=:BilayerNetwork)
  simp_bilayer_example = ODESystem(eqs, t, name=:BilayerNetwork)
  petri_example = ODESystem(eqs, t, name=:PetriNet)

  @test ODESystem(psir) == petri_example
  @test ODESystem(bnsir) == bilayer_example
  @test ODESystem(bnsir, simplify = true) == simp_bilayer_example
  @test typeof(ODESystem(PetriNet(psir))) == ODESystem # Sanity check that non-labelled Petri Nets resolve correctly
end
