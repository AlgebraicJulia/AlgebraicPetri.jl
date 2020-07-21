using Remark, FileWatching

while true
    Remark.slideshow(@__DIR__; options = Dict("ratio" => "16:9"), title = "AlgebraicPetri.jl: COEXIST")
    @info "Rebuilt"
    FileWatching.watch_folder(joinpath(@__DIR__, "src"))
end
