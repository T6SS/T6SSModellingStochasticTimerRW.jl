using ModellingStochasticTimerRW
using Documenter

DocMeta.setdocmeta!(ModellingStochasticTimerRW, :DocTestSetup, :(using ModellingStochasticTimerRW); recursive=true)

makedocs(;
    modules=[ModellingStochasticTimerRW],
    authors="Jonathan Miller",
    repo="https://github.com/fieldofnodes/ModellingStochasticTimerRW.jl/blob/{commit}{path}#{line}",
    sitename="ModellingStochasticTimerRW.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fieldofnodes.github.io/ModellingStochasticTimerRW.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fieldofnodes/ModellingStochasticTimerRW.jl",
    devbranch="main",
)
