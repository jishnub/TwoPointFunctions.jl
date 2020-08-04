using Documenter
using TwoPointFunctions

DocMeta.setdocmeta!(TwoPointFunctions, :DocTestSetup, :(using TwoPointFunctions); recursive=true)

makedocs(;
    modules=[TwoPointFunctions],
    authors="Jishnu Bhattacharya",
    repo="https://github.com/jishnub/TwoPointFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="TwoPointFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jishnub.github.io/TwoPointFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Reference" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jishnub/TwoPointFunctions.jl",
)
