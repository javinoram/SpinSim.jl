The objetive of this module is give some functions that facilite the process of construct 
spin matrices. Some functions help to construct some spin operators and commonly use spin terms that 
it can be see in the literature.

## Contruct Spin matrices
There are 5 function to construct the spin matrices. Four of them are to construct the pauli matrix.

```julia
x_matrix( size::Int )
y_matrix( size::Int )
z_matrix( size::Int )
i_matrix( size::Int )
```

Where 'size' is a integer that represent the size of the square matrix, this size is obtained from $2S +1$ where $S$
is the spin.

The next function return a dictionary of all the four pauli matrices of a particular spin, this function get as an input 
a float value that represent the spin site.


```julia
pauli_matrices( spin::Float64 )
```

## Construct terms
The next function allows us to compute the matrix asociated a particular term in a spin hamiltonian.

```julia
construct_term( operator::String, spins::Array{Float64,1} )
```

This function require two parameters, an string that represent the action of the spin matrices in each site of the 
system (for example "ZZI", that can be readed as the action of the Z operator on sites 1 and 2, and the I operator 
on the site 3). The second one is an 1D array that represent the spin of each site of the system.

The lenght of both had to be the same, because they represent the number of sites in the system.


## Construct hamiltonian
The final step is construct the complete matrix of the system. This function use all the before functions.

```julia
construct_hamiltonian( list_operators::Array{Tuple{Float64, String},1}, spins::Vector{Float64} )
```

This function require two parameters, one is a 1D array of tuples, the tuple correspond a term of the hamiltonian with their 
respective numerical value (see the toy example). The second parameter is also a 1D array with the spin value in each site.
As same as the above function, the lenght of the string and the lenght of the array of spin must match.

