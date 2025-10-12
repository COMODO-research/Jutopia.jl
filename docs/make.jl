using Documenter
using DocumenterVitepress
using Jutopia

makedocs(
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/COMODO-research/Jutopia.jl.git"
    ),
    sitename = "Jutopia.jl",
    modules = [Jutopia],
    pages = [
        "Home" =>[
            "index.md"
            ],
        "Tutorials"  =>[
        "tutorials/getting-started.md"
        ],
        "Examples"  =>[

        ]
    ]
)

DocumenterVitepress.deploydocs(
    repo = "github.com/COMODO-research/Jutopia.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)

