using LinearAlgebra

@testset "Pauli matrices" begin
    ϵ = 1e-7 

    @testset "0.5" begin
      @test [1 0; 0 1] == SpinSim.hamiltonian.i_matrix(2)
      @test [0 0.5; 0.5 0] == SpinSim.hamiltonian.x_matrix(2)
      @test [0 -0.5im; 0.5im 0] == SpinSim.hamiltonian.y_matrix(2)
      @test [0.5 0; 0 -0.5] == SpinSim.hamiltonian.z_matrix(2)
    end
    
    @testset "1" begin
      @test real( sum( [1 0 0; 0 1 0; 0 0 1] - SpinSim.hamiltonian.i_matrix(3) ) ) ≤ ϵ
      @test real( sum( (1.0/LinearAlgebra.sqrt(2))*[0 1 0; 1 0 1; 0 1 0] - SpinSim.hamiltonian.x_matrix(3) ) ) ≤ ϵ
      @test real( sum( (1.0im/LinearAlgebra.sqrt(2))*[0 -1 0; 1 0 -1; 0 1 0] - SpinSim.hamiltonian.y_matrix(3) ) ) ≤ ϵ
      @test real( sum( [1 0 0; 0 0 0; 0 0 -1] - SpinSim.hamiltonian.z_matrix(3) ) ) ≤ ϵ
    end
end
