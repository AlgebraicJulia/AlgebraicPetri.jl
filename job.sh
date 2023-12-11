#!/bin/bash

julia --project -e 'using Pkg; Pkg.status; Pkg.test()' > log_test.md
