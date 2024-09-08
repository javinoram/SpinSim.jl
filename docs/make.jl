using Documenter

makedocs(; format=Documenter.HTML(),
         sitename="SpinSim.jl",
         pages=["Home" => "index.md"])

deploydocs(; repo="github.com/javinoram/SpinSim.jl.git", target="build")
