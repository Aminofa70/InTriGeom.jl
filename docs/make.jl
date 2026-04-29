using InTriGeom
using Documenter

DocMeta.setdocmeta!(InTriGeom, :DocTestSetup, :(using InTriGeom); recursive=true)

makedocs(;
    modules=[InTriGeom],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="InTriGeom.jl",
    format=Documenter.HTML(;
        canonical="https://Aminofa70.github.io/InTriGeom.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Aminofa70/InTriGeom.jl",
    devbranch="main",
)
