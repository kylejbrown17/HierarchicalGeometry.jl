using HierarchicalGeometry

using LazySets
using StaticArrays

using Test, Logging
# Set logging level
global_logger(SimpleLogger(stderr, Logging.Debug))

# Check equality of two arrays
@inline function array_isapprox(x::AbstractArray{F},
                  y::AbstractArray{F};
                  rtol::F=sqrt(eps(F)),
                  atol::F=zero(F)) where {F<:AbstractFloat}

    # Easy check on matching size
    if length(x) != length(y)
        return false
    end

    for (a,b) in zip(x,y)
        @test isapprox(a,b, rtol=rtol, atol=atol)
    end
end

# Check if array equals a single value
@inline function array_isapprox(x::AbstractArray{F},
                  y::F;
                  rtol::F=sqrt(eps(F)),
                  atol::F=zero(F)) where {F<:AbstractFloat}

    for a in x
        @test isapprox(a, y, rtol=rtol, atol=atol)
    end
end

@testset "HierarchicalGeometry.jl" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    @time @testset "HierarchicalGeometry.Overapproximation" begin
        include(joinpath(testdir, "test_approximations.jl"))
    end
end
