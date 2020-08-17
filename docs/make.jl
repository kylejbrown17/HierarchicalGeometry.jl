using PolyhedralOverapproximation
using Documenter

makedocs(;
    modules=[PolyhedralOverapproximation],
    authors="kylebrown <kylejbrown17@gmail.com> and contributors",
    repo="https://github.com/kylejbrown17/PolyhedralOverapproximation.jl/blob/{commit}{path}#L{line}",
    sitename="PolyhedralOverapproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kylejbrown17.github.io/PolyhedralOverapproximation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kylejbrown17/PolyhedralOverapproximation.jl",
)
