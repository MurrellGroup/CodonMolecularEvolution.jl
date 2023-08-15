using CodonMolecularEvolution
using Documenter

DocMeta.setdocmeta!(CodonMolecularEvolution, :DocTestSetup, :(using CodonMolecularEvolution); recursive=true)

makedocs(;
    modules=[CodonMolecularEvolution],
    authors="Ben Murrell <murrellb@gmail.com> and contributors",
    repo="https://github.com/murrellb/CodonMolecularEvolution.jl/blob/{commit}{path}#{line}",
    sitename="CodonMolecularEvolution.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://murrellb.github.io/CodonMolecularEvolution.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/murrellb/CodonMolecularEvolution.jl",
    devbranch="main",
)
