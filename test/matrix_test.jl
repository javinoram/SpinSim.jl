using LinearAlgebra

@testset "Pauli matrices" begin
    ϵ = 1e-7 

    @testset "Spin 0.5" begin
      @test [1 0; 0 1] == SpinSim.hamiltonian.i_matrix(2)
      @test [0 0.5; 0.5 0] == SpinSim.hamiltonian.x_matrix(2)
      @test [0 -0.5im; 0.5im 0] == SpinSim.hamiltonian.y_matrix(2)
      @test [0.5 0; 0 -0.5] == SpinSim.hamiltonian.z_matrix(2)
    end
    
    @testset "Spin 1" begin
      @test real( sum( [1 0 0; 0 1 0; 0 0 1] - SpinSim.hamiltonian.i_matrix(3) ) ) ≤ ϵ
      @test real( sum( (1.0/LinearAlgebra.sqrt(2))*[0 1 0; 1 0 1; 0 1 0] - SpinSim.hamiltonian.x_matrix(3) ) ) ≤ ϵ
      @test real( sum( (1.0im/LinearAlgebra.sqrt(2))*[0 -1 0; 1 0 -1; 0 1 0] - SpinSim.hamiltonian.y_matrix(3) ) ) ≤ ϵ
      @test real( sum( [1 0 0; 0 0 0; 0 0 -1] - SpinSim.hamiltonian.z_matrix(3) ) ) ≤ ϵ
    end
end

@testset "Spin operators" begin
  ϵ = 1e-7
  
  #tests with three 0.5 spin system XII, ZII and YII 
  @testset "Single operators" begin
    matrix_z = 0.5*[1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 -1 0 0 0; 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 -1 0; 0 0 0 0 0 0 0 -1]
    matrix_x = 0.5*[0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1; 1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0]
    matrix_y = 0.5im*[0 0 0 0 -1 0 0 0; 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 -1 0; 0 0 0 0 0 0 0 -1; 1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0]

    @test matrix_z == SpinSim.hamiltonian.construct_term( "ZII", [0.5, 0.5, 0.5] )
    @test matrix_x == SpinSim.hamiltonian.construct_term( "XII", [0.5, 0.5, 0.5] )
    @test matrix_y == SpinSim.hamiltonian.construct_term( "YII", [0.5, 0.5, 0.5] )
  end
  
  #tests with three 0.5 spin system XII, ZII and YII  
  @testset "Double operators" begin
    matrix_zz = 0.25*[1 0 0 0 0 0 0 0; 0 -1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 -1 0 0 0 0; 0 0 0 0 -1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 -1 0; 0 0 0 0 0 0 0 1]
    matrix_xx = 0.25*[0 0 0 1 0 0 0 0; 0 0 1 0 0 0 0 0; 0 1 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 1 0; 0 0 0 0 0 1 0 0; 0 0 0 0 1 0 0 0]
    matrix_yy = -0.25*[0 0 0 1 0 0 0 0; 0 0 -1 0 0 0 0 0; 0 -1 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1; 0 0 0 0 0 0 -1 0; 0 0 0 0 0 -1 0 0; 0 0 0 0 1 0 0 0]

    @test sum( matrix_zz - real( SpinSim.hamiltonian.construct_term( "ZZI", [0.5, 0.5, 0.5] ) ) ) ≤ ϵ
    @test sum( matrix_xx - real( SpinSim.hamiltonian.construct_term( "XXI", [0.5, 0.5, 0.5] ) ) ) ≤ ϵ
    @test sum( matrix_yy - real( SpinSim.hamiltonian.construct_term( "YYI", [0.5, 0.5, 0.5] ) ) ) ≤ ϵ
  end
  
end

@testset "Hamiltonians" begin
  spin_list = [0.5, 0.5, 0.5]
  ϵ = 1e-7

  @testset "Ising model" begin
    hamil_base_terms = [(-1.0, "ZZI"), (-1.0, "IZZ")]
    hamil_func = SpinSim.hamiltonian.construct_hamiltonian(hamil_base_terms, spin_list)

    hamil_exact = -1.0*SpinSim.hamiltonian.construct_term( "ZZI", spin_list ) -1.0*SpinSim.hamiltonian.construct_term( "IZZ", spin_list )

    @test real( sum( hamil_exact - hamil_func ) ) ≤ ϵ
  end

  @testset "Heisenberg Model" begin
    hamil_base_terms = [(-1.0, "ZZI"), (-1.0, "IZZ"), (-1.0, "XXI"), (-1.0, "IXX"), (-1.0, "YYI"), (-1.0, "IYY")]
    hamil_func = SpinSim.hamiltonian.construct_hamiltonian(hamil_base_terms, spin_list)
    
    hamil_exact = -1.0*SpinSim.hamiltonian.construct_term( "ZZI", spin_list ) -1.0*SpinSim.hamiltonian.construct_term( "IZZ", spin_list )
    hamil_exact = hamil_exact -1.0*SpinSim.hamiltonian.construct_term( "XXI", spin_list ) -1.0*SpinSim.hamiltonian.construct_term( "IXX", spin_list )
    hamil_exact = hamil_exact -1.0*SpinSim.hamiltonian.construct_term( "YYI", spin_list ) -1.0*SpinSim.hamiltonian.construct_term( "IYY", spin_list )

    @test real( sum( hamil_exact - hamil_func ) ) ≤ ϵ
  end

end
