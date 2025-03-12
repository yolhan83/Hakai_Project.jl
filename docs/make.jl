using Hakai_Project
using Documenter

DocMeta.setdocmeta!(Hakai_Project, :DocTestSetup, :(using Hakai_Project); recursive=true)

makedocs(;
    modules=[Hakai_Project],
    authors="yolhan83 <yolhan@laposte.net>",
    sitename="Hakai_Project.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
