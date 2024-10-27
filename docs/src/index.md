# SpinSim.jl

In the different sections you can find some examples and explanations of how to use the code.
Welcome to this small package to construct the matrix representation of some spin model.
This matrices are expected to be used for magnetocaloric studies or analize the spectral decomposition.


This is the recopilation of all the work done to simulate and publish the result of the 3Ni SMM.

## Installation

`SpinSim.jl` will be soon registed in the `General` registry.

```julia-repl
pkg> add SpinSim
```

## Toy example
For our toy example we will use a Ising model of four spin sites and open boundary conditions. 
The spin value of each site will ve 0.5 and 1.0.

```julia
using SpinSim

#Hamiltonian description
const spin_list = [0.5, 0.5, 0.5]
const J = 1.0
const hamiltonian_descripction = [(-J, "ZZII"), (-J, "IZZI"), (-J, "IIZZ")]

#Matriz representation of the hamiltonian
matrix_representation = SpinSim.hamiltonian.construct_hamiltonian(hamiltonian_descripction, spin_list)
```

