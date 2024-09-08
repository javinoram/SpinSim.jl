using LinearAlgebra

module hamiltonian

function delta_kronecker( i::Int, j::Int )
  if (i==j)
    return 1
  end
  return 0
end


function x_matrix( size::Int )
  spin = (size-1)/2.0
  base = zeros( (size, size) )

  for a in 1:size
    for b in 1:size
      base[a][b] = 0.5*(delta_kronecker(a, b+1) + delta_kronecker(a+1, b))*sqrt( (spin+1)*(a+b-1) -a*b )
    end
  end

  return base
end

function y_matrix( size::Int )
  spin = (size-1)/2.0
  base = zeros( (size, size) )
  for a in 1:size
    for b in 1:size
      base[a][b] = 0.5im*(delta_kronecker(a, b+1) - delta_kronecker(a+1, b))*sqrt( (spin+1)*(a+b-1) -a*b )
    end
  end

end

function z_matrix( size::Int )
  spin = (size-1)/2.0
  base = zeros( (size, size) )
  for a in 1:size
    for b in 1:size
      base[a][b] = delta_kronecker(a, b)*(spin+1-a) 
    end
  end

end

function i_matrix( size::Int )
  base = Matrix{Float64}(I, size, size)
  return base
end


function pauli_matrices( spin::Float64 )
  matrix_size = Integer( 2*spin +1 )
  dict = Dict( "I"=>i_matrix( matrix_size ), "X"=>x_matrix( matrix_size ),
      "Y"=>y_matrix( matrix_size ), "Z"=>z_matrix( matrix_size ) )
  return dict
end


function construct_term( operator::String, spins::Array{Float64,1} )
  unique_spins = unique( spins )
  pauli_dict = Dict()
  #Basis element for the construction of the operator
  for s in unique_spins
    pauli_dict[s] = pauli_matrices(s) 
  end
  
  #Sequential Kronecker product 
  result = pauli_dict[ spins[1] ][ operator[1] ]
  for i in 2:length(operator)
    result = kron( result, pauli_dict[ spins[i] ][ operator[i] ] )
  end
  
  return result 
end


function construct_hamiltonian( list_operators::Array{Tuple{Float64, String},1}, spins::Array{Float64,1} )

  size = Integer( prod( 2*spins .+ 1 ) )
  base = zeros( (size, size) )

  for (j,op) in list_operators
    base += j*construct_term( op, spins )
  end

  return base
end


end

