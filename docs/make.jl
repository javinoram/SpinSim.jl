push!(LOAD_PATH, "../src/")
using Documenter
makedocs(
    sitename = "SpinSim.jl",
    pages=[
      "Home" => "index.md"
      "Hamiltonian Module" => "hamiltonian.md"
    ])

deploydocs(; 
  repo="github.com/javinoram/SpinSim.jl",
)
