using Documenter
makedocs(
    sitename = "SpinSim.jl",
    pages=[
      "Home" => "index.md"
    ])

deploydocs(; 
  repo="github.com/javinoram/SpinSim.jl",
)