# # [Converting Models to Computation Graphs](@id bilayernetwork_example)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/covid/bilayernetworks.ipynb)
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
migrate!(bnsir, psir)
bnsir
to_graphviz(bnsir)
