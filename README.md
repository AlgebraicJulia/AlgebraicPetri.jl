# AlgebraicPetri.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicPetri.jl/stable)
[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://AlgebraicJulia.github.io/AlgebraicPetri.jl/dev)
[![CI/CD](https://github.com/AlgebraicJulia/AlgebraicPetri.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/AlgebraicPetri.jl/actions/workflows/julia_ci.yml)
[![DOI](https://zenodo.org/badge/275202510.svg)](https://zenodo.org/badge/latestdoi/275202510)

`AlgebraicPetri.jl` is a Julia package for building Petri Net models
compositionally. This package depends on [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl)
which provides a framework for applied category theory (ACT) in Julia, and is
part of the [AlgebraicJulia organization](https://www.algebraicjulia.org/),
which develops ACT-based software to improve scientific and technical
computing.

`AlgebraicPetri.jl` defines the category of open Petri nets as described in [[Baez 2018](https://arxiv.org/abs/1808.05415)],
and implements composition and stratification methods for such Petri nets
from [[Libkind 2022]](https://doi.org/10.1098/rsta.2021.0309).

## Getting started

Please visit the [documentation](https://algebraicjulia.github.io/AlgebraicPetri.jl/dev/)
to learn more about how to use the package. For more details on the underlying theory,
please consult [[Baez 2018](https://arxiv.org/abs/1808.05415)] or [[Libkind 2022]](https://doi.org/10.1098/rsta.2021.0309),
with the latter being more oriented towards applied practitioners.

Tutorials are organized by theme. Currently examples from epidemiology are the most
developed. The basic tutorial introduces the main themes of modeling with open Petri nets,
followed by a detailed tutorial on model stratification, and a detailed example
of developing a complex model with multiple strains of a directly-transmitted disease
with vaccines.