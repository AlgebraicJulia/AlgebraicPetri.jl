module ModelingToolkitInterop
  using Test
  using ModelingToolkit
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

  make_depvar(p, t) = :($p($t))

  @assert make_depvar(:a, :t) == :(a(t))
  @assert make_depvar(:y, :x) == :(y(x))

  #= /Users/fairbanksj/github/AlgebraicJulia/ASKEM-demos/MTK/petri.jl:39 =#
  @variables t S(t) I(t) R(t)
  @parameters inf rec
  D = Differential(t)
  ϕ1 = inf * S * I
  ϕ2 = rec * I
  eqs = [D(S) ~ +(-ϕ1), D(I) ~ ϕ1 + ϕ1 + -ϕ1 + -ϕ2, D(R) ~ +ϕ2]
  example = ODESystem(eqs, t, name=:PetriNet)
  @test ODESystem(bnsir) == example
end
