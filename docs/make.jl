using Documenter
using Literate

@info "Loading AlgebraicPetri"
using AlgebraicPetri

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/AlgebraicPetri.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/AlgebraicPetri.jl/blob/master/docs"
end

const literate_dir = joinpath(@__DIR__, "..", "examples")
const generated_dir = joinpath(@__DIR__, "src", "examples")

for (root, dirs, files) in walkdir(literate_dir)
  out_dir = joinpath(generated_dir, relpath(root, literate_dir))
  for file in files
    f,l = splitext(file)
    if l == ".jl" && !startswith(f, "_")
      Literate.markdown(joinpath(root, file), out_dir;
        config=config, documenter=true, credit=false)
      Literate.notebook(joinpath(root, file), out_dir;
        execute=true, documenter=true, credit=false)
    end
  end
end

@info "Building Documenter.jl docs"
makedocs(
  modules   = [AlgebraicPetri],
  format    = Documenter.HTML(
    assets = ["assets/analytics.js"],
  ),
  sitename  = "AlgebraicPetri.jl",
  doctest   = false,
  checkdocs = :none,
  pages     = Any[
    "AlgebraicPetri.jl" => "index.md",
    "Examples" => Any[
      "examples/predation/lotka-volterra.md",
      "examples/covid/epidemiology.md",
      "examples/covid/coexist/coexist.md",
      "examples/enzymes/enzyme_reactions.md",
    ],
    "Library Reference" => "api.md",
  ]
)

@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/AlgebraicPetri.jl.git",
  branch = "gh-pages"
)
