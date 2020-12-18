using HierarchicalGeometry
using Test

using LazySets
using StaticArrays

@testset "HierarchicalGeometry.jl" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    @time @testset "HierarchicalGeometry.Overapproximation" begin
        include(joinpath(testdir, "test_approximations.jl"))
    end
end
