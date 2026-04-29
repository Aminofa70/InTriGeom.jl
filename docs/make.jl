using InTriGeom
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(InTriGeom, :DocTestSetup, :(using InTriGeom); recursive=true)

makedocs(;
    modules=[InTriGeom],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="InTriGeom.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="https://github.com/Aminofa70/InTriGeom.jl",
        devbranch="main",
       devurl="dev",
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Getting started" => "getting-started.md",
            
             "Demos" => [
                "Demo 0001" => "demo_0001.md",
            ],
        ], "API Reference" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo="https://github.com/Aminofa70/InTriGeom.jl.git",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true
)