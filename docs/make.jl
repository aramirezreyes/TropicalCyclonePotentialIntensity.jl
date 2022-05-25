using TropicalCyclonePotentialIntensity
using Documenter

DocMeta.setdocmeta!(TropicalCyclonePotentialIntensity, :DocTestSetup, :(using TropicalCyclonePotentialIntensity); recursive=true)

makedocs(;
    modules=[TropicalCyclonePotentialIntensity],
    authors="Argel Ramirez Reyes <aramirezreyes@ucdavis.edu> and contributors",
    repo="https://github.com/aramirezreyes/TropicalCyclonePotentialIntensity.jl/blob/{commit}{path}#{line}",
    sitename="TropicalCyclonePotentialIntensity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aramirezreyes.github.io/TropicalCyclonePotentialIntensity.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aramirezreyes/TropicalCyclonePotentialIntensity.jl",
    devbranch="main",
)
