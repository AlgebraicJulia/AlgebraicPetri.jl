#!/bin/bash

# julia --project -e 'using Pkg; Pkg.instantiate(;verbose=true)'
# julia --project -e 'using Pkg; Pkg.precompile'
# julia --project -e 'using Pkg; Pkg.status'


julia --project test/runtests.jl > log_test.md
