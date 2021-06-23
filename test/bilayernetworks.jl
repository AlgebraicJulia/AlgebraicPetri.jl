using Test

using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks
using LabelledArrays

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs.RelationalPrograms
display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));


# #### SIR Model:

# define model
sir = @relation (s,i,r) begin
    infection(s,i)
    recovery(i,r)
end
display_uwd(sir)
#-
#

# Extract the Petri network form of the model.

psir = apex(oapply_epi(sir))
psir
Graph(psir)

# Convert the Petri network into a bilayer network and draw it.
# This model uses a computation graph to express the computation of the vector field of the Petri net.

bnsir = BilayerNetwork()
lab_bnsir = LabelledBilayerNetwork()
migrate!(bnsir, psir)
migrate!(lab_bnsir, psir)
bnsir
@test typeof(to_graphviz(bnsir)) == Graph
@test typeof(to_graphviz(lab_bnsir)) == Graph

bnsir_test = @acset BilayerNetwork begin
    Qin = 3
    Qout = 3
    Win = 3
    Wa = 3
    Wn = 3
    Box = 2
    arg = [1,2,2]
    call = [1,1,2]
    efflux = [1,1,2]
    effusion = [1,2,2]
    influx = [1,1,2]
    infusion = [2,2,3]
end

@test bnsir == bnsir_test


######################
# Round Trip testing #
######################
function roundtrip(pn::AbstractPetriNet, bn::AbstractBilayerNetwork)
    roundtrippetri = PetriNet()
    migrate!(roundtrippetri, bn)
    pn_structure = PetriNet()
    copy_parts!(pn_structure, pn)
    return roundtrippetri, pn_structure
end

function roundtrip(pn::AbstractLabelledPetriNet, bn::AbstractLabelledBilayerNetwork)
    roundtrippetri = LabelledPetriNet()
    migrate!(roundtrippetri, bn)
    return roundtrippetri, pn
end

function test_roundtrip(pn::AbstractPetriNet, bn::AbstractBilayerNetwork)
    roundtrippetri, pn_structure = roundtrip(pn, bn)
    @test roundtrippetri == pn_structure
end

bnsir = BilayerNetwork()
migrate!(bnsir, psir)
test_roundtrip(psir, bnsir)

seir = @relation (s,e,i,r) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
end

pseir = apex(oapply_epi(seir))
bnseir = BilayerNetwork()
migrate!(bnseir, pseir)
test_roundtrip(pseir, bnseir)

lbnseir = LabelledBilayerNetwork()
migrate!(lbnseir, pseir)
test_roundtrip(pseir, lbnseir)

########################
# Compile and Evaluate #
########################

du = LVector(S=0.0,I=0.0,R=0.0)
u = LVector(S=10.0,I=1.0,R=0.0)
params = LVector(inf=0.1,rec=0.05)

bn_du =  AlgebraicPetri.BilayerNetworks.evaluate(lab_bnsir, u; inf=0.1,rec=0.05)
vectorfield(psir)(du, u, params, 0)

@test all(abs.(bn_du .- du[[:S,:I,:R]]) .< 1e-9)

eval(AlgebraicPetri.BilayerNetworks.compile(lab_bnsir, :du, :ϕ, :u; inf=0.1, rec=0.05))
comp_du = [0.0,0.0,0.0]
comp_ϕ = [0.0,0.0,0.0]

f!(comp_du, comp_ϕ, u, 0)

@test all(abs.(comp_du .- du[[:S,:I,:R]]) .< 1e-9)

eval(AlgebraicPetri.BilayerNetworks.compile(lab_bnsir, :du, :ϕ, :u, :params))
comp_du = [0.0,0.0,0.0]
comp_ϕ = [0.0,0.0,0.0]

f!(comp_du, comp_ϕ, u, params, 0)

@test all(abs.(comp_du .- du[[:S,:I,:R]]) .< 1e-9)

##############
# Edge cases #
##############

function test_sir_equivalence(bn, pn)
  du = LVector(S=0.0,I=0.0,R=0.0)
  u = LVector(S=10.0,I=1.0,R=0.0)
  params = LVector(inf=0.1,rec=0.05)

  bn_du =  AlgebraicPetri.BilayerNetworks.evaluate(bn, u; inf=0.1,rec=0.05)
  vectorfield(pn)(du, u, params, 0)

  @test all(abs.(bn_du .- du[[:S,:I,:R]]) .< 1e-9)
end

# Ensure that the migration function properly balances the BLN
bnsir_edge = @acset LabelledBilayerNetwork begin
    Qin = 3
    Qout = 3
    Win = 3
    Wa = 2
    Wn = 2
    Box = 2
    arg = [1,2,2]
    call = [1,1,2]
    efflux = [1,2]
    effusion = [1,2]
    influx = [1,2]
    infusion = [2,3]
    parameter=[:inf, :rec]
    variable=[:S, :I, :R]
    tanvar=[:S, :I, :R]
end

test_sir_equivalence(bnsir_edge, psir)

edge_sir = LabelledPetriNet()
migrate!(edge_sir, bnsir_edge)
test_sir_equivalence(bnsir_edge, edge_sir)

bnsir_edge_add = @acset LabelledBilayerNetwork begin
    Qin = 3
    Qout = 3
    Win = 3
    Wa = 4
    Wn = 4
    Box = 2
    arg = [1,2,2]
    call = [1,1,2]
    efflux = [1,1,1,2]
    effusion = [1,2,2,2]
    influx = [1,1,1,2]
    infusion = [2,2,2,3]
    parameter=[:inf, :rec]
    variable=[:S, :I, :R]
    tanvar=[:S, :I, :R]
end

test_sir_equivalence(bnsir_edge_add, psir)

edge_sir = LabelledPetriNet()
migrate!(edge_sir, bnsir_edge_add)
test_sir_equivalence(bnsir_edge, edge_sir)



bnsir_migrate_break = @acset LabelledBilayerNetwork begin
    Qin = 3
    Qout = 3
    Win = 2
    Wa = 2
    Wn = 2
    Box = 2
    arg = [1,2]
    call = [1,1]
    efflux = [1,2]
    effusion = [1,2]
    influx = [1,2]
    infusion = [2,3]
    parameter=[:inf, :rec]
    variable=[:S, :I, :R]
    tanvar=[:S, :I, :R]
end
bn_error = "Mass action does not allow species $(bnsir_migrate_break[2, :variable]) to be "*
           "removed without contributing to the rate."
@test_throws ErrorException(bn_error) migrate!(edge_sir, bnsir_migrate_break)
