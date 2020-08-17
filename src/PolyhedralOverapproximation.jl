module PolyhedralOverapproximation

# Write your package code here.
using StaticArrays
using LinearAlgebra
using LazySets
using Polyhedra
# using QHull

"""
    cross_product_operator(x)

Computes `Ωₓ` such that `Ωₓy = x ⨯ y`
"""
cross_product_operator(x) = SMatrix{3,3}(
    [0.0    -x[3]   x[2];
     x[3]   0.0     -x[1];
     -x[2]  x[1]    0.0]
)

abstract type OverapproxModel end

"""
    RegularPolyhedronOverapprox{D,N}

Used to overrapproximate convex sets with a pre-determined set of support
vectors.
"""
struct RegularPolyhedronOverapprox{D,N} <: OverapproxModel
    support_vectors::NTuple{N,SVector{D,Float64}}
end
function RegularPolyhedronOverapprox(dim::Int,N::Int)
    vecs = Vector{SVector{dim,Float64}}()
    # v - the initial support vector for a given dimension
    v = zeros(dim)
    v[end] = 1.0
    for i in 1:dim
        d = zeros(dim)
        d[i] = 1.0
        A = cross_product_operator(d/N)
        for j in 1:N
            v = exp(A*j)*v
            push!(vecs,v)
        end
        v = d
    end
    RegularPolyhedronOverapprox(tuple(vecs))
end

function get_support_vectors(model::)

function LazySets.overapproximate(model::RegularPolyhedronOverapprox{N},pts) where {D,N}
end

end
