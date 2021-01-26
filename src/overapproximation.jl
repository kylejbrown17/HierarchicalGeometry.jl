
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

LazySets.HPolytope(m::PolyhedronOverapprox) = HPolytope(map(v->LazySets.HalfSpace(v,1.0),get_support_vecs(m)))

# LazySets.ρ(d::AbstractVector,geom::AbstractVector{N}) where {N<:GeometryBasics.Ngon} = mapfoldl(v->ρ(d,v),max,geom)
LazySets.ρ(d::AbstractVector,geom::AbstractVector{N}) where {N<:GeometryBasics.Ngon} = maximum(map(v->ρ(d,v),geom))


function LazySets.overapproximate(lazy_set,model::H,epsilon::Float64=0.1) where {V,T,H<:HPolytope{T,V}}
    halfspaces = map(h->LazySets.HalfSpace(h.a, ρ(h.a, lazy_set)), constraints_list(model))
    sort!(halfspaces; by = h->h.b)
    halfspaces = sort(LazySets.constraints_list(model); by=h->ρ(h.a, lazy_set))
    # poly = HPolyhedron{T,V}()
    # while !isempty(halfspaces) #&& !isbounded(poly)
    #     poly = intersection(poly,halfspaces[1])
    #     deleteat!(halfspaces,1)
    # end
    # @assert isbounded(poly)
    # hpoly = convert(H,poly)
    hpoly = H()
    for h in halfspaces
        addconstraint!(hpoly,h)
    end
    hpoly
end
LazySets.overapproximate(lazy_set,m::PolyhedronOverapprox,args...) = overapproximate(lazy_set,HPolytope(m),args...)

for TYPE in (:VPolytope,:HPolytope)
    @eval begin
        function LazySets.overapproximate(p::H,m::Hyperrectangle) where {V,T,H<:$TYPE{T,V}}
            high = -Inf*ones(V)
            low = Inf*ones(V)
            for v in vertices_list(p)
                high = max.(high,v)
                low = min.(low,v)
            end
            ctr = (high .+ low) / 2
            widths = (high .- low) / 2
            Hyperrectangle(V(ctr),V(widths))
        end
    end
end

"""
    LazySets.center(p::AbstractPolytope)

A hacky way of choosing a reasonable center for a polytope.
"""
LazySets.center(p::AbstractPolytope) = LazySets.center(overapproximate(p))
# function LazySets.overapproximate(lazy_set::AbstractPolytope,sphere::H) where {H<:BallType}
function LazySets.overapproximate(lazy_set::AbstractPolytope,sphere::H) where {V,T,H<:Ball2{T,V}}
    r = maximum(map(v->norm(v-get_center(sphere)), vertices_list(lazy_set)))
    H(get_center(sphere),r)
end
# function LazySets.overapproximate(s::AbstractPolytope,sphere::Type{H}) where {H<:BallType}
function LazySets.overapproximate(s::AbstractPolytope,sphere::Type{H}) where {V,T,H<:Ball2{T,V}}
    overapproximate(s,H(LazySets.center(s),1.0))
end

export
    GridDiscretization,
    GridOccupancy

struct GridDiscretization{N,T}
    origin::SVector{N,T}
    discretization::SVector{N,T}
end
get_hyperrectangle(m::GridDiscretization,idxs) = Hyperrectangle(m.origin .+ idxs.*m.discretization, [m.discretization/2...])
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
    ctr = (s.center .+ f.center) / 2
    radii = (f.center .- s.center .+ o.grid.discretization) / 2
    Hyperrectangle(ctr ,radii)
end
function LazySets.overapproximate(lazy_set,grid::GridDiscretization)
    rect = overapproximate(lazy_set,Hyperrectangle)
    starts = LazySets.center(rect) .- radius_hyperrectangle(rect)
    stops = LazySets.center(rect) .+ radius_hyperrectangle(rect)
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

# select robot carry locations
"""
    construct_support_placement_aggregator

Construct an objective function that scores a selection of indices into `pts`.
Balances a neighbor-neighbor metric with a all-all metric
Args
* pts - a vector of points
* n - the number of indices to be selected
* [f_neighbor = v->1.0*minimum(v)] - a function mapping a vector of distances
    to a scalar value.
* [f_inner = v->1.0*minimum(v)] - a function mapping a vector of distances
    to a scalar value.
"""
function construct_support_placement_aggregator(pts, n,
        f_neighbor=v->1.0*minimum(v)+(0.5/n)*sum(v),
        f_inner=v->(0.1/(n^2))*minimum(v)
    )
    D = [norm(v-vp) for (v,vp) in Base.Iterators.product(pts,pts)]
    d_neighbor = (idxs)->f_neighbor(
        map(i->wrap_get(D,(idxs[i],wrap_get(idxs,i+1))),1:length(idxs))
        )
    d_inner = (idxs)->f_inner(
        [wrap_get(D,(i,j)) for (i,j) in Base.Iterators.product(idxs,idxs)]
        )
    d = (idxs)->d_neighbor(idxs)+d_inner(idxs)
end
"""
    spaced_neighbors(polygon,n::Int,aggregator=sum)

Return the indices of the `n` vertices of `polygon` whose neighbor distances
maximize the utility metric defined by `aggregator`. Uses local optimization,
so there is no guarantee of global optimality.
"""
function spaced_neighbors(polygon,n::Int,
        score_function=construct_support_placement_aggregator(vertices_list(polygon), n),
        ϵ=1e-8)
    pts = vertices_list(polygon)
    @assert length(pts) >= n "length(pts) = $(length(pts)), but n = $n"
    if length(pts) == n
        return collect(1:n)
    end
    best_idxs = SVector{n,Int}(collect(1:n)...)
    d_hi = score_function(best_idxs)
    idx_list = [best_idxs]
    while true
        updated = false
        for deltas in Base.Iterators.product(map(i->(-1,0,1),1:n)...)
            idxs = sort(map(i->wrap_idx(length(pts),i),best_idxs .+ deltas))
            @show idxs
            if length(unique(idxs)) == n
                if score_function(idxs) > d_hi + ϵ
                    best_idxs = map(i->wrap_idx(length(pts),i),idxs)
                    d_hi = score_function(idxs)
                    @show best_idxs, d_hi
                    push!(idx_list,best_idxs)
                    updated = true
                end
            end
        end
        if !updated
            break
        end
    end
    return best_idxs, d_hi
end
function extremal_points(pts)
    D = [norm(v-vp) for (v,vp) in Base.Iterators.product(pts,pts)]
    val, idx = findmax(D)
    return val, idx.I
end
proj_to_line(v,vec) = vec*dot(v,vec)/(norm(vec)^2)
function proj_to_line_between_points(p,p1,p2)
    v = p.-p1
    vec = normalize(p2-p1)
    p1 .+ proj_to_line(v,vec)
end
perimeter(pts) = sum(map(i->norm(wrap_get(pts,(i,i+1))),1:length(pts)))
perimeter(p::LazySets.AbstractPolygon) = perimeter(vertices_list(p))
get_pts(p::LazySets.AbstractPolytope) = vertices_list(p)
get_pts(v::AbstractVector) = v
"""
    select_support_locations(geom,transport_model)

Given some arbitrary 3D geometry, select a set of locations to support it from
beneath. Requires specification of an aggregator function.
"""
function select_support_locations(geom,transport_model,
        score_function_constructor=construct_support_placement_aggregator,
    )
    r       = transport_model.robot_radius
    a_r_max     = transport_model.max_area_per_robot
    v_r_max = transport_model.max_volume_per_robot

    zvec = SVector{3,Float64}(0,0,1)
    proj_mat = one(SMatrix{3,3,Float64})[1:2,1:3]
    # gpts = get_pts(geom)
    polygon = VPolygon(convex_hull(map(v->proj_mat*v,geom)))
    # compute geom height, projection area, and bounding volume
    H = maximum(map(v->sum(dot(v,zvec)),geom)) # height
    A = LazySets.area(polygon) # area
    P = perimeter(polygon)
    V = H*A # volume
    # compute rquired number of robots
    n = Int(ceil(max(A / a_r_max, V / v_r_max)))
    n = min(n,Int(ceil(P/(2*r))))
    if n == 1
        support_pts = [LazySets.center(polygon)]
    else
        pts = vertices_list(polygon)
        score_function=score_function_constructor(pts,n)
        best_idxs, _ = spaced_neighbors(polygon,n,score_function)
        @show best_idxs
        support_pts = pts[best_idxs]
    end
    # length
    # L, (i,j) = extremal_points(pts)
    # width
    # p1 = pts[i]
    # p2 = pts[j]
    # vecs = map(p->p.-p1,pts)
    # vec = normalize(p2-p1)
    # dists = map(v->norm(v-proj_to_line(v,vec)),vecs)
    # k = argmax(dists)
    support_pts

end
select_support_locations(p::AbstractPolytope,args...) = select_support_locations(vertices_list(p),args...)

# spread out sub assemblies
