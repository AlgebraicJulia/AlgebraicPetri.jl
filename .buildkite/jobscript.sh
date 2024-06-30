pwd; hostname; date

module load julia

echo "Running Tests..."
julia --project -t 32 -e 'using Pkg; Pkg.status(); Pkg.test()'

echo "Building Documentation..."
julia --project=docs -t 32 -e'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.status(); Pkg.instantiate(); include("docs/make.jl")'
