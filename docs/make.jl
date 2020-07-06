using Documenter

@info "Loading AlgebraicPetri"
using AlgebraicPetri

@info "Building Documenter.jl docs"
makedocs(
  modules   = [AlgebraicPetri],
  format    = Documenter.HTML(),
  sitename  = "AlgebraicPetri.jl",
  doctest   = false,
  checkdocs = :none,
  pages     = Any[
    "AlgebraicPetri.jl" => "index.md",
    "Basic Usage" => "usage.md",
    "Library Reference" => "api.md",
  ]
)

@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/AlgebraicPetri.jl.git",
  branch = "gh-pages"
)
