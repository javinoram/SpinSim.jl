using Test
using LinearAlgebra
using SpinSim


@testset begin
    filenames = ["matrix_test.jl"]
    for filename in filenames
        t = @elapsed include(filename)
        println("$(filename): $t sec")
    end
end
