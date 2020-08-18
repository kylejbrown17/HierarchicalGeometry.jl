module PolyhedralOverapproximation

# Write your package code here.
using StaticArrays
using LinearAlgebra
using Polyhedra
using LazySets
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

export
    RegularPolyhedronOverapprox

"""
    RegularPolyhedronOverapprox{D,N}

Used to overrapproximate convex sets with a pre-determined set of support
vectors.
"""
struct RegularPolyhedronOverapprox{D,N} <: OverapproxModel
    support_vectors::NTuple{N,SVector{D,Float64}}
end
get_support_vecs(model::RegularPolyhedronOverapprox) = [v for v in model.support_vectors]

"""
    RegularPolyhedronOverapprox(dim::Int,N::Int,epsilon=0.1)

Construct a regular polyhedron overapproximation model by specifying the number
of dimensions and the number of support vectors to be arranged radially about
the axis formed by the unit vector along each dimension.
epsilon if the distance between support vector v and an existing support vector
is less than epsilon, the new vector will not be added.
"""
function RegularPolyhedronOverapprox(dim::Int,N::Int,epsilon=0.1)
    vecs = Vector{SVector{dim,Float64}}()
    # v - the initial support vector for a given dimension
    v = zeros(dim)
    v[end] = 1.0
    for i in 1:dim
        d = zeros(dim)
        d[i] = 1.0
        A = cross_product_operator(d)
        for j in 1:N
            v = normalize(exp(A*j*2*pi/N)*v)
            add = true
            for vp in vecs
                if norm(v-vp) < epsilon
                    add = false
                    break
                end
            end
            if add
                @show v
                push!(vecs,v)
            end
        end
        v = d
    end
    RegularPolyhedronOverapprox(tuple(vecs...))
end


function LazySets.overapproximate(lazy_set,model::RegularPolyhedronOverapprox{N},epsilon::Float64=0.1) where {D,N}
    halfspaces = map(v->LazySets.HalfSpace(v,ρ(v,lazy_set)),get_support_vecs(model))
    sort!(halfspaces; by=h->ρ(h.a, lazy_set))
    poly = HPolyhedron()
    while !isempty(halfspaces) && !isbounded(poly)
        poly = intersection(poly,halfspaces[1])
        deleteat!(halfspaces,1)
    end
    @assert isbounded(poly)
    hpoly = convert(HPolytope,poly)
    # vpoly = tovrep(hpoly)
    # continue adding new halfspace constraints until the distance "clipped off"
    # by the new halfspace constraint is less than epsilon
    # while !isempty(halfspaces)
        for v in vertices_list(hpoly)

        end
    # end
    hpoly
end

# """
#     sort_facet_vtxs
#
# Returns a counterclockwise ordering of 'vtx_ids' around the face defined by
# the halfspace associated with h_idx. It is assumed that the vertices form a
# convex polygon
# """
# function sort_facet_vtxs(poly,h_idx,vtx_idxs)
#     if !isempty(vtx_idxs)
#         vtxs = map(v->get(poly,v), vtx_idxs)
#         h = get(poly,h_idx) # halfspace
#         ctr = sum(vtxs)
#         # axes
#         a1 = normalize(vtxs[1]-ctr)
#         a2 = normalize(cross(h.a,a1))
#         # angles = map(v->atan(dot(a2,v),dot(a1,v)))
#         return sort(vtx_idxs;by=v->atan(dot(a2,get(poly,v)),dot(a1,get(poly,v))))
#     end
#     return vtx_idxs
# end

# """
#     polytope_vtxs_and_edges
#
# Identify the edges of each facet of a polytope
# """
# function polytope_vtxs_and_edges(hpolytope::AbstractPolyhedron)
#     poly = polyhedron(hpolytope)
#     vtxs = collect(points(poly))
#     edges = Vector{Tuple{Int,Int}}()
#     for idx in eachindex(halfspaces(poly))
#         vtx_idxs = incidentpointindices(poly,idx)
#         sorted_vtx_ids = sort_facet_vtxs(poly,idx,vtx_idxs)
#         for i in 1:length(sorted_vtx_ids)-1
#             v1 = sorted_vtx_ids[i]
#             v2 = sorted_vtx_ids[i+1]
#             push!(edges, (v1.value, v2.value))
#         end
#         push!(edges,(sorted_vtx_ids[end].value,sorted_vtx_ids[1].value))
#     end
#     return vtxs, edges
# end
#
# function write_vertices(filename,vtxs::Vector{Vector{Float64}})
#     open(filename, "w") do io
#         for v in vtxs
#             for c in v
#                 print(io,c," ")
#             end
#             print(io,"\n")
#         end
#     end
# end
#
# function write_edges(filename,edge_list::Vector{Tuple{Int,Int}})
#     open(filename, "w") do io
#         for e in edge_list
#             for v in e
#                 print(io,v," ")
#             end
#             print(io,"\n")
#         end
#     end
# end
#
# function write_vtxs_and_edges(filename,poly)
#     vtx_file = join([filename,".nodes"])
#     edge_file = join([filename,".edges"])
#     vtxs, edge_list = polytope_vtxs_and_edges(poly)
#     write_vertices(vtx_file,vtxs)
#     write_edges(edge_file,edge_list)
# end

end
