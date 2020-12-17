using HierarchicalGeometry
using Documenter

makedocs(;
    modules=[HierarchicalGeometry],
    authors="kylebrown <kylejbrown17@gmail.com> and contributors",
    repo="https://github.com/kylejbrown17/HierarchicalGeometry.jl/blob/{commit}{path}#L{line}",
    sitename="HierarchicalGeometry.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kylejbrown17.github.io/HierarchicalGeometry.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kylejbrown17/HierarchicalGeometry.jl",
)
