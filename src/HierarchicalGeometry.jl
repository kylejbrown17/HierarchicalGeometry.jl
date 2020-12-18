module HierarchicalGeometry

using StaticArrays
using LinearAlgebra
using Polyhedra
using LazySets
using GeometryTypes

"""
    cross_product_operator(x)

Computes `Ωₓ` such that `Ωₓy = x ⨯ y`
"""
cross_product_operator(x) = SMatrix{3,3}(
    [0.0    -x[3]   x[2];
     x[3]   0.0     -x[1];
     -x[2]  x[1]    0.0]
)

"""
    GeometryHierarchy

A hierarchical representation of geometry
Fields:
* graph - encodes the hierarchy of geometries
* geoms -
"""
# struct GeometryHierarchy
#
# end
export
    GridDiscretization,
    GridOccupancy,
    has_overlap

struct GridDiscretization{N,T}
    origin::SVector{N,T}
    discretization::SVector{N,T}
end
get_hyperrectangle(m::GridDiscretization,idxs) = Hyperrectangle(m.origin .+ idxs.*m.discretization, m.discretization/2)
"""
    cell_indices(m::GridDiscretization,v)

get indices of cell of `m` in which `v` falls
"""
cell_indices(m::GridDiscretization,v) = SVector(ceil.(Int,(v .- m.origin .- m.discretization/2)./m.discretization)...)
struct GridOccupancy{N,T,A<:AbstractArray{Bool,N}}
    grid::GridDiscretization{N,T}
    occupancy::A
    offset::SVector{N,Int}
end
GridOccupancy(m::GridDiscretization{N,T},o::AbstractArray) where {N,T} = GridOccupancy(m,o,SVector(zeros(Int,N)...))
Base.:(+)(o::GridOccupancy,v) = GridOccupancy(o.grid,o.occupancy,SVector(o.offset.+v...))
Base.:(-)(o::GridOccupancy,v) = o+(-v)
get_hyperrectangle(m::GridOccupancy,idxs) = get_hyperrectangle(m.grid,idxs .+ m.offset)
function Base.intersect(o1::G,o2::G) where {G<:GridOccupancy}
    offset = o2.offset - o1.offset
    starts = max.(1,offset .+ 1)
    stops = min.(SVector(size(o1.occupancy)),size(o2.occupancy) .+ offset)
    idxs = CartesianIndex(starts...):CartesianIndex(stops...)
    overlap = o1.occupancy[idxs] .* o2.occupancy[idxs .- CartesianIndex(offset...)]
    G(o1.grid,overlap,o2.offset)
end
has_overlap(o1::G,o2::G) where {G<:GridOccupancy} = any(intersect(o1,o2).occupancy)
LazySets.is_intersection_empty(o1::G,o2::G) where {G<:GridOccupancy} = !has_overlap(o1,o2)
function LazySets.overapproximate(o::GridOccupancy,::Type{Hyperrectangle})
    origin = o.grid.origin
    start = findnext(o.occupancy,CartesianIndex(ones(Int,size(origin))...))
    finish = findprev(o.occupancy,CartesianIndex(size(o.occupancy)...))
    s = get_hyperrectangle(o,start.I .- 1)
    f = get_hyperrectangle(o,finish.I .- 1)
    center = (s.center .+ f.center) / 2
    radii = (f.center .- s.center .+ o.grid.discretization) / 2
    Hyperrectangle(center,radii)
end
function LazySets.overapproximate(lazy_set,grid::GridDiscretization)
    rect = overapproximate(lazy_set,Hyperrectangle)
    starts = center(rect) .- radius_hyperrectangle(rect)
    stops = center(rect) .+ radius_hyperrectangle(rect)
    @show start_idxs = cell_indices(grid,starts)
    @show stop_idxs = cell_indices(grid,stops)
    @show offset = SVector(start_idxs...) .- 1
    occupancy = falses((stop_idxs .- offset)...)
    approx = GridOccupancy(grid,occupancy,offset)
    for idx in CartesianIndices(occupancy)
        r = get_hyperrectangle(approx,[idx.I...])
        if !is_intersection_empty(lazy_set,r)
            approx.occupancy[idx] = true
        end
    end
    GridOccupancy(grid,occupancy,offset)
end

abstract type OverapproxModel end

export PolyhedronOverapprox
"""
    PolyhedronOverapprox{D,N}

Used to overrapproximate convex sets with a pre-determined set of support
vectors.
"""
struct PolyhedronOverapprox{D,N} <: OverapproxModel
    support_vectors::NTuple{N,SVector{D,Float64}}
end
get_support_vecs(model::PolyhedronOverapprox) = [v for v in model.support_vectors]

"""
    PolyhedronOverapprox(dim::Int,N::Int,epsilon=0.1)

Construct a regular polyhedron overapproximation model by specifying the number
of dimensions and the number of support vectors to be arranged radially about
the axis formed by the unit vector along each dimension.
epsilon if the distance between support vector v and an existing support vector
is less than epsilon, the new vector will not be added.
"""
function PolyhedronOverapprox(dim::Int,N::Int,epsilon=0.1)
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
    PolyhedronOverapprox(tuple(vecs...))
end

export equatorial_overapprox_model
"""
    equatorial_overapprox_model(lat_angles=[-π/4,0.0,π/4],lon_angles=collect(0:π/4:2π),epsilon=0.1)

Returns a PolyhedronOverapprox model generated with one face for each
combination of pitch and yaw angles specified by lat_angles and lon_angles,
respectively. There are also two faces at the two poles
"""
function equatorial_overapprox_model(lat_angles=[-π/4,0.0,π/4],lon_angles=collect(0:π/4:2π),epsilon=0.1)
    vecs = Vector{SVector{3,Float64}}()
    push!(vecs,[0.0,0.0,1.0])
    push!(vecs,[0.0,0.0,-1.0])
    for phi in lat_angles
        for theta in lon_angles
            v = normalize([
                cos(phi)*cos(theta),
                cos(phi)*sin(theta),
                sin(phi)
            ])
            add = true
            for vp in vecs
                if norm(v-vp) < epsilon
                    add = false
                    break
                end
            end
            if add
                push!(vecs,v)
            end
        end
    end
    PolyhedronOverapprox(tuple(vecs...))
end

function LazySets.overapproximate(lazy_set,model::PolyhedronOverapprox{N},epsilon::Float64=0.1) where {D,N}
    halfspaces = map(v->LazySets.HalfSpace(v,ρ(v,lazy_set)),get_support_vecs(model))
    sort!(halfspaces; by=h->ρ(h.a, lazy_set))
    poly = HPolyhedron()
    while !isempty(halfspaces) #&& !isbounded(poly)
        poly = intersection(poly,halfspaces[1])
        deleteat!(halfspaces,1)
    end
    @assert isbounded(poly)
    hpoly = convert(HPolytope,poly)
    # vpoly = tovrep(hpoly)
    # continue adding new halfspace constraints until the distance "clipped off"
    # by the new halfspace constraint is less than epsilon
    # while !isempty(halfspaces)
    #     for v in vertices_list(hpoly)
    #
    #     end
    # end
    hpoly
end

function LazySets.overapproximate(lazy_set::Polyhedron,sphere::H) where {H<:HyperSphere}
    r = maximum(map(v->norm(v-origin(sphere)), points(vrep(lazy_set))))
    HyperSphere(origin(sphere),r)
end
# function LazySets.overapproximate(lazy_set::Polyhedron,rect::H) where {H<:HyperRectangle}
#     r_p = maximum(points(vrep(lazy_set)))
#     r_n = minimum(points(vrep(lazy_set)))
#     HyperSphere(origin(sphere),r)
# end

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
