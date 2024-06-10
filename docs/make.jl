using CodonMolecularEvolution
using Documenter

DocMeta.setdocmeta!(CodonMolecularEvolution, :DocTestSetup, :(using CodonMolecularEvolution); recursive=true)

makedocs(;
    modules=[CodonMolecularEvolution],
    authors="Ben Murrell <benjamin.murrell@ki.se> and contributors",
    repo="https://github.com/MurrellGroup/CodonMolecularEvolution.jl/blob/{commit}{path}#{line}",
    sitename="CodonMolecularEvolution.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MurrellGroup.github.io/CodonMolecularEvolution.jl",
        edit_link="main",
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
        "difFUBAR.md",
    ],
)

deploydocs(;
    repo="github.com/MurrellGroup/CodonMolecularEvolution.jl",
    devbranch="main",
)
