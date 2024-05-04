# # [Maximum Common Sub C-Set](@id max_common_subobject)
#
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/other/max_common_subobject.ipynb)

using AlgebraicPetri
using Catlab.CategoricalAlgebra, Catlab.Graphics

# ### Define epidemiology models

# #### SIR

sir   = LabelledPetriNet([:S, :I, :R],
                         :inf => ((:S, :I) => (:I, :I)),
                         :rec => (:I => :R)
                        )
to_graphviz(sir)

# #### SIRD

sird  = LabelledPetriNet([:S, :I, :R, :D],
                         :inf => ((:S, :I) => (:I, :I)),
                         :rec => (:I => :R),
                         :death => (:R => :D)
                        )
to_graphviz(sird)

# #### SEIR

seir  = LabelledPetriNet([:S, :E, :I, :R],
                         :exp => ((:S, :I) => (:E, :I)),
                         :ill => (:E => :I),
                         :rec => (:I => :R)
                        )
to_graphviz(seir)

# #### SEIRD

seird = LabelledPetriNet([:S, :E, :I, :R, :D],
                         :exp => ((:S, :I) => (:E, :I)),
                         :ill => (:E => :I),
                         :rec => (:I => :R),
                         :death => (:R => :D)
                        )
to_graphviz(seird)

# ### Calculate the Maximum Common C-Set

sub, overlaps = maximum_common_subobject(sir, sird, seir, seird) |> first
morphisms = legs(overlaps[1])

to_graphviz(sub)

# ### Visualize Sub C-Sets

# The morphism from each model to the maximum common C-Set is depicted in bold. 

# #### SIR

morphisms[1](sub) |> to_graphviz

# #### SIRD

morphisms[2](sub) |> to_graphviz

# #### SEIR

morphisms[3](sub) |> to_graphviz

# #### SEIRD

morphisms[4](sub) |> to_graphviz
