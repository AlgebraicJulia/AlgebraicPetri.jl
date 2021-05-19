# AlgebraicPetri Examples

This directory contains several projects which rely on the tooling present in
AlgebraicPetri, and which can be taken as examples of what can be done with the
tooling in AlgebraicPetri. These projects are not large enough to warrant their
own library, yet they do not contribute to the core functionality of AlgebraicPetri.

## Execution Instructions

There are some basic instructions which should work for running the demo
file present in each example directory. 

1. Clone repository
1. Change to the desired example's folder (e.g. `AlgebraicPetri/examples/covid`)
1. From that directory run: `julia --project setup.jl`
1. Run the example script: `julia --project demo.jl`

## Standards

Each examples directory should include the following files:
- Project.toml: This file stores any external dependencies of the example
- setup.jl: This file sets up the user environment, adding any necessary
  packages and generating potential test data
- demo.jl: This file provides a demonstration of the tooling developed in this package

Extra files and directories may be included as it is relevant to the particular example
