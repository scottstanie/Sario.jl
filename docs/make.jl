using Documenter, Sario

makedocs(;
    modules=[Sario],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/scottstanie/Sario.jl/blob/{commit}{path}#L{line}",
    sitename="Sario.jl",
    authors="scott <scott.stanie@gmail.com>",
    assets=String[],
)

deploydocs(
    repo = "github.com/scottstanie/Sario.jl.git",
)
