__precompile__()

module SpinSim
  using LinearAlgebra

  include("./hamiltonian/hamiltonian.jl")
  include("./information/info.jl")
  include("./operators/operators.jl")
  include("./thermodinamic/thermal.jl")
end
