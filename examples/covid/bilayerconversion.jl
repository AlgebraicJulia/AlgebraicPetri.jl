# # [Converting Models to Computation Graphs](@id bilayernetwork_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/bilayerconversion.ipynb)


# using Pkg
# Pkg.status()
# Pkg.activate(".")

# Pkg.instantiate()

# Pkg.develop(path="../../")
# Pkg.instantiate()

using AlgebraicPetri
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs.RelationalPrograms
display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));

using PrettyTables
function printsoln(bn::AbstractLabelledBilayerNetwork, soln::Vector)
    pretty_table(soln)
end

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
#-

# Convert the Petri network into a bilayer network and draw it.
# This model uses a computation graph to express the computation of the vector field of the Petri net.
bnsir = BilayerNetwork()
migrate!(bnsir, psir)
bnsir
to_graphviz(bnsir)
#-

# We can hand code a Bilayer netowork using the `@acset` macro provided by Catlab. As you can see from the code,
# there is a lot of typing to specify the incidence of all these wires. The Petri Net is more compact.
# This notion of Bilayer network comes from the definition of mass action kinetics for reaction networks.

# hand coded Bilayer network
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

@assert bnsir == bnsir_test

function roundtrip(pn::AbstractPetriNet, bn::AbstractBilayerNetwork)
    roundtrippetri = PetriNet()
    migrate!(roundtrippetri, bn)
    pn_structure = PetriNet()
    copy_parts!(pn_structure, pn)
    return roundtrippetri, pn_structure
end
#-

# We can specify the SEIR model with the relation macro and functorial semantics as usual.

seir = @relation (s,e,i,r) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
end

pseir = apex(oapply_epi(seir))

# By converting this to a Bilayer network, we are able to visualize differences in the computational pattern
# of data flow between different reaction network models.

#
bnseir = BilayerNetwork()
migrate!(bnseir, pseir)

bnrt,pnstr = roundtrip(pseir, bnseir)

display_uwd(sir)

Graph(psir)
#-

to_graphviz(bnsir)
#-

to_graphviz(bnseir)
#-

Graph(bnrt)
#-

qm = @relation (s,q) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
    quarexp(e,q)
    quarinf(i,q)
end
display_uwd(qm)
#-

import AlgebraicPetri.Epidemiology: exposure_petri, spontaneous_petri
semantics = Dict(
    :infection => exposure_petri(:S, :I, :I, :inf),
    :exposure  => exposure_petri(:S, :I, :E, :exp),
    :illness   => spontaneous_petri(:E,:I,:ill),
    :recovery  => spontaneous_petri(:I,:R,:rec),
    :death     => spontaneous_petri(:I,:D,:death),
    :quarexp   => spontaneous_petri(:E, :Q, :qe),
    :quarinf   => spontaneous_petri(:I, :Q, :qi),
    :quarrec   => spontaneous_petri(:Q, :R, :qr)
)
pn_quar = oapply(qm, semantics)  |> apex
Graph(pn_quar)
#-

bnquar = LabelledBilayerNetwork()
migrate!(bnquar, pn_quar)
to_graphviz(bnquar)

qm = @relation (s,i,q) begin
    exposure(s,i,e)
    illness(e,i)
    recovery(i,r)
    quarexp(e,q)
    quarinf(i,q)
    quarrec(q,r)
end
display_uwd(qm)
#-

pn_quar = oapply(qm, semantics)  |> apex
Graph(pn_quar)
#-

bnquar = LabelledBilayerNetwork()
migrate!(bnquar, pn_quar)
to_graphviz(bnquar)
#-

quarrt = LabelledPetriNet()
migrate!(quarrt, bnquar) |> Graph

bnquar
#-

import AlgebraicPetri.BilayerNetworks: evaluate!, evaluate

bnsir = LabelledBilayerNetwork()
migrate!(bnsir, psir)

evaluate(bnsir, [10.0, 1, 0.0], inf=0.1, rec=0.3)

evaluate(bnquar, [10.0, 1, 0, 0,0,0], exp=0.1, rec=0.3, qi=0.2, ill=0.7, qe=0.23, qr=0.3)

function euler(bn::AbstractLabelledBilayerNetwork, state, nsteps::Integer, stepsize::Real; params...)
    #preallocate storage so that each step is nonallocating
    #create a storage space for steps of euler's method store intermediate
    #states as named tuples so that you can integrate with julia Tables.jl data  analysis ecosystem
    du = zeros(length(state))
    ϕ = ones(nparts(bn, :Box))
    u = tuple(state...)
    results = Vector{NamedTuple{tuple(bn[:,:variable]...)}}()
    for i in 1:nsteps
        u = u .+ stepsize.*evaluate!(du, ϕ, bn, u; params...)
        push!(results, NamedTuple{tuple(bn[:,:variable]...)}(u))
    end
    return results
end

soln = euler(bnquar, (S=10.0, I=1, E=0, R=0, Q=0), 30, 0.15, exp=0.1, rec=0.03, qi=0.37, ill=0.7, qe=0.23, qr=0.03)
printsoln(bnquar, soln)
#-

import AlgebraicPetri.BilayerNetworks: compile
compile(bnquar, :du, :ϕ, :u, :p)
#-
#
compile(bnquar, :du, :ϕ, :u, exp=0.1, rec=0.03, qi=0.37, ill=0.7, qe=0.23, qr=0.03)
#-

function eulers(bn::AbstractLabelledBilayerNetwork, funcname::Symbol; params...)
    f = compile(bn, :du, :ϕ, :u; params...)
    varnames = tuple(bn[:,:variable]...)
    nϕ = nparts(bn, :Box)
    quote
    function $funcname(state, nsteps::Integer, stepsize::Real)
        $f
        #preallocate storage so that each step is nonallocating
        du = zeros(length(state))
        ϕ = ones($nϕ)
        u = tuple(state...)
        #create a storage space for steps of euler's method
        #store intermediate states as named tuples so that you can integrate
        #with julia Tables.jl data analysis ecosystem
        results = Vector{NamedTuple{$varnames}}()
        for i in 1:nsteps
            Δ = f!(du, ϕ, u, 0)
            u = u .+ stepsize.*Δ
            push!(results, NamedTuple{$varnames}(u))
        end
        return results
    end 
    end
end

eulers(bnsir, :eulsir, inf=0.3, rec=0.2)
#-

eulseirqexp = eulers(bnquar, :eulseirq, exp=0.1, rec=0.03, qi=0.37, ill=0.7, qe=0.23, qr=0.03)
#-

eulsirexp = eulers(bnsir, :eulsir, inf=0.3, rec=0.2)
#-

eval(eulsirexp)
soln_codegen = eulsir((S=10.0, I=1, R=0), 30, 0.15)
pretty_table(soln_codegen)
#-

eval(eulseirqexp)
soln_codegen = eulseirq((S=10.0, I=1, E=0, R=0, Q=0), 30, 0.15)
pretty_table(soln_codegen)
#-
