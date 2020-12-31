module HierarchicalGeometry

using StaticArrays
using LinearAlgebra
using Polyhedra
using LazySets
using GeometryBasics
using CoordinateTransformations
using Rotations
using LightGraphs, GraphUtils
using Parameters

export
    transform,
    transform!,
    distance_lower_bound,
    has_overlap

const BaseGeometry = Union{Ball2,Hyperrectangle,AbstractPolytope}
const BallType = Ball2
get_center(s::Ball2) = LazySets.center(s)
get_radius(s::Ball2) = s.radius
GeometryBasics.Sphere(s::Ball2) = GeometryBasics.Sphere(Point(s.center...),s.radius)
const RectType = Hyperrectangle
get_center(s::Hyperrectangle) = s.center
get_radius(s::Hyperrectangle) = s.radius
GeometryBasics.HyperRectangle(s::Hyperrectangle) = GeometryBasics.HyperRectangle(Vec((s.center .- s.radius)...),2*Vec(s.radius...))

"""
    transform(geom,t)

Transform geometry `geom` according to the transformation `t`.
TODO: It is important for all base geometries to be defined with respect to
their own origin.
"""
transform(v::V,t) where {V<:AbstractVector} = V(t(v))
function transform!(v::V,t) where {V<:AbstractVector}
    v .= transform(v,t)
end
transform(h::LazySets.HalfSpace,t::CoordinateTransformations.LinearMap{R}) where {R<:Rotation} = LazySets.HalfSpace(transform(Vector(h.a),t),h.b)
function transform!(h::LazySets.HalfSpace,t)
    h.a .= transform!(h).a
    return h
end
transform(g::BaseGeometry,t::CoordinateTransformations.Translation) = LazySets.translate(g,Vector(t.translation))
transform!(g::BaseGeometry,t::CoordinateTransformations.Translation) = LazySets.translate!(g,t.translation)
transform(g::Ball2,t::CoordinateTransformations.LinearMap) = Ball2(transform(g.center,t),g.radius)
function transform!(g::Ball2,t::CoordinateTransformations.LinearMap)
    transform!(g.center,t)
    return g
end
"""
    "rotating" a Hyperrectangle `g` results in a new Hyperrectangle that bounds the
    transformed version `g`.
"""
transform(g::Hyperrectangle,t::CoordinateTransformations.LinearMap) = overapproximate(transform(convert(LazySets.VPolytope,g),t),Hyperrectangle)
function transform(g::VPolytope,t::CoordinateTransformations.LinearMap{R}) where {R<:Rotation}
    VPolytope(map(v->transform(v,t), vertices_list(g)))
end
function transform(g::HPolytope,t::CoordinateTransformations.LinearMap{R}) where {R<:Rotation}
    HPolytope(map(h->transform(h,t), constraints_list(g)))
end
function transform(g::BaseGeometry,t::CoordinateTransformations.AffineMap{R,T}) where {R<:Rotation,T}
    transform(
        transform(g,CoordinateTransformations.LinearMap(t.linear)),
        CoordinateTransformations.Translation(t.translation)
        )
end

export
    distance_lower_bound,
    has_overlap

LazySets.ρ(a::AbstractArray,p::GeometryBasics.Point) = ρ(a,Singleton([p.data...]))
LazySets.ρ(a::AbstractArray,v::AbstractArray{V,1}) where {V<:GeometryBasics.Point} = minimum(map(x->ρ(a,x),v))
LazySets.ρ(a::AbstractArray,x::GeometryBasics.Ngon) = ρ(a,coordinates(x))
# Distance between sets
LazySets.distance(a::BallType,b::BallType,p::Real=2.0) = norm(get_center(a)-get_center(b),p) - (get_radius(a)+get_radius(b))
function LazySets.distance(a::Hyperrectangle,b::Hyperrectangle,p::Real=2.0)
    m = minkowski_sum(a,b)
    d = map(h->dot(h.a,get_center(a))-h.b, constraints_list(m))
    sort!(d;rev=true)
    idx = findlast(d .> 0.0)
    if idx === nothing
        return d[1]
    else
        return norm(d[1:idx])
    end
end
distance_lower_bound(a::BallType,b::BallType) = LazySets.distance(a,b)
distance_lower_bound(a::Hyperrectangle,b::Hyperrectangle) = LazySets.distance(a,b)
has_overlap(a::BaseGeometry,b::BaseGeometry) = !is_intersection_empty(a,b) # distance_lower_bound(a,b) <= 0.0

export GeometryCollection
struct GeometryCollection{G}
    elements::Vector{G}
end
for op in [:distance_lower_bound,:(LazySets.distance)]
    @eval $op(a::GeometryCollection,b) = minimum($op(x,b) for x in a.elements)
    @eval $op(a::GeometryCollection,b::GeometryCollection) = minimum($op(a,y) for y in b.elements)
end

export
    TransformNode,
    local_transform,
    global_transform,
    set_local_transform!,
    set_global_transform!

mutable struct TransformNode
    local_transform::CoordinateTransformations.Transformation # transform from the parent frame to the child frame
    global_transform::CoordinateTransformations.Transformation # transform from the base frame to the child frame
end
function TransformNode()
    TransformNode(
    CoordinateTransformations.LinearMap(RotZ(0)) ∘ CoordinateTransformations.Translation(0,0,0),
    CoordinateTransformations.LinearMap(RotZ(0)) ∘ CoordinateTransformations.Translation(0,0,0),
    )
end
local_transform(n::TransformNode) = n.local_transform
global_transform(n::TransformNode) = n.global_transform
"""
    relative_transform(a::AffineMap,b::AffineMap)

compute the relative transform between two frames, ᵂT𝐀 and ᵂT𝐁.
"""
function relative_transform(a::CoordinateTransformations.AffineMap,b::CoordinateTransformations.AffineMap,error_map=MRPMap())
    pos_err = CoordinateTransformations.Translation(b.translation - a.translation)
    rot_err = Rotations.rotation_error(a,b,error_map)
    q = UnitQuaternion(rot_err)
    t = pos_err ∘ CoordinateTransformations.LinearMap(q) # IMPORTANT: rotation on right side
end
function Rotations.rotation_error(
    a::CoordinateTransformations.AffineMap,
    b::CoordinateTransformations.AffineMap,
    error_map=MRPMap())
    rot_err = Rotations.rotation_error(
        UnitQuaternion(a.linear),
        UnitQuaternion(b.linear),error_map)
end

function set_local_transform!(n::TransformNode,t)
    n.local_transform = t
end
function set_global_transform!(n::TransformNode,t)
    n.global_transform = t
end
const transform_node_accessor_interface = [:local_transform,:global_transform]
const transform_node_mutator_interface = [:set_local_transform!,:set_global_transform!]

@with_kw struct TransformTree{ID} <: AbstractCustomTree{TransformNode,ID}
    graph               ::DiGraph               = DiGraph()
    nodes               ::Vector{TransformNode} = Vector{TransformNode}()
    vtx_map             ::Dict{ID,Int}          = Dict{ID,Int}()
    vtx_ids             ::Vector{ID}            = Vector{ID}()
end
for op in transform_node_accessor_interface
    @eval $op(g::TransformTree,v) = $op(get_node(g,v))
end
function get_parent_transform(g::TransformTree,v)
    vp = get_parent(g,v)
    if has_vertex(g,vp)
        return global_transform(get_node(g,vp))
    end
    return global_transform(TransformNode())
end
"""
    update_transform_tree!(g::TransformTree,v)

Updates the global transforms of `v` and its successors based on their local
transforms.
"""
function update_transform_tree!(g::TransformTree,v)
    n = get_node(g,v)
    set_global_transform!(n,get_parent_transform(g,v) ∘ local_transform(n))
    for vp in outneighbors(g,v)
        update_transform_tree!(g,vp)
    end
    return g
end
function set_local_transform!(g::TransformTree,v,t,update=true)
    set_local_transform!(get_node(g,v),t)
    if update
        update_transform_tree!(g,v)
    end
    return g
end

"""
    check_collision(a,b)

Recursive collision checking
"""
function check_collision(a,b)
    if has_overlap(a,b)
        for p in components(b)
            if check_collision(p,a)
                return true
            end
        end
    end
    return false
end

export
    CachedElement,
    get_element,
    is_up_to_date,
    set_up_to_date!,
    set_element!,
    update_element!

"""
    CachedElement{G}

A mutable container for geometries
"""
mutable struct CachedElement{G}
    element::G
    is_up_to_date::Bool
end
CachedElement(element) = CachedElement(element,false)
get_element(n::CachedElement) = n.element
is_up_to_date(n::CachedElement) = n.is_up_to_date
function set_up_to_date!(n::CachedElement,val::Bool=true)
    n.is_up_to_date = val
end
function set_element!(n::CachedElement,g)
    n.element = g
    return g
end
function update_element!(n::CachedElement,g)
    set_element!(n,g)
    set_up_to_date!(n,true)
    return g
end
# accessor_interface(::Type{CachedElement}) = [:is_up_to_date]
# mutator_interface(::Type{CachedElement}) = [:set_computed_geometry!,:set_up_to_date!]
const cached_element_accessor_interface = [:is_up_to_date]
const cached_element_mutator_interface = [:update_element!,:set_element!,:set_up_to_date!]

export
    GeomNode,
    get_base_geom,
    get_cached_geom

struct GeomNode{G}
    base_geom::G
    transform_node::TransformNode
    cached_geom::CachedElement{G}
end
function GeomNode(geom)
    GeomNode(geom,TransformNode(),CachedElement(geom))
end
get_base_geom(n::GeomNode) = n.base_geom
for op in cached_element_accessor_interface
    @eval $op(g::GeomNode) = $op(g.cached_geom)
end
for op in cached_element_mutator_interface
    @eval $op(g::GeomNode,val) = $op(g.cached_geom,val)
end
for op in transform_node_accessor_interface
    @eval $op(g::GeomNode) = $op(g.transform_node)
end
function set_local_transform!(n::GeomNode,t)
    set_local_transform!(n.transform_node,t)
    set_up_to_date!(n,false)
end
function set_global_transform!(n::GeomNode,t)
    set_global_transform!(n.transform_node,t)
    set_up_to_date!(n,false)
end
function get_cached_geom(n::GeomNode)
    if !is_up_to_date(n)
        transformed_geom = transform(get_base_geom(n),global_transform(n))
        update_element!(n,transformed_geom)
    end
    return get_element(n.cached_geom)
end
distance_lower_bound(a::GeomNode{G},b::GeomNode{G}) where {G<:Union{BallType,RectType}} = distance_lower_bound(get_cached_geom(a),get_cached_geom(b))


export GeometryHierarchy
"""
    GeometryHierarchy

A hierarchical representation of geometry
Fields:
* graph - encodes the hierarchy of geometries
* nodes - geometry nodes
"""
@with_kw struct GeometryHierarchy <: AbstractCustomDiGraph{GeomNode,Symbol}
    graph               ::DiGraph               = DiGraph()
    nodes               ::Vector{GeomNode}      = Vector{GeomNode}()
    vtx_map             ::Dict{Symbol,Int}      = Dict{Symbol,Int}()
    vtx_ids             ::Vector{Symbol}        = Vector{Symbol}() # maps vertex uid to actual graph node
end

distance_lower_bound(a::GeometryHierarchy,b::GeometryHierarchy) = distance_lower_bound(get_node(a,:Hypersphere),get_node(b,:Hypersphere))
function has_overlap(a::GeometryHierarchy,b::GeometryHierarchy,leaf_id=:Hypersphere)
    v = get_vtx(a,leaf_id)
    while has_vertex(a,v)
        k = get_vtx_id(a,v)
        if has_vertex(b,k)
            if !has_overlap(get_node(a,k),get_node(b,k))
                return false
            end
        end
        v = get_parent(a,k)
    end
    return true
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

LazySets.HPolytope(m::PolyhedronOverapprox) = HPolytope(map(v->LazySets.HalfSpace(v,1.0),get_support_vecs(m)))

function LazySets.overapproximate(lazy_set,model::HPolytope,epsilon::Float64=0.1)
    halfspaces = map(h->LazySets.HalfSpace(h.a, ρ(h.a, lazy_set)), constraints_list(model))
    sort!(halfspaces; by = h->h.b)
    # halfspaces = sort(LazySets.constraints_list(model); by=h->ρ(h.a, lazy_set))
    poly = HPolyhedron()
    while !isempty(halfspaces) #&& !isbounded(poly)
        poly = intersection(poly,halfspaces[1])
        deleteat!(halfspaces,1)
    end
    @assert isbounded(poly)
    hpoly = convert(HPolytope,poly)
    hpoly
end
LazySets.overapproximate(lazy_set,m::PolyhedronOverapprox,args...) = overapproximate(lazy_set,HPolytope(m),args...)




"""
    LazySets.center(p::AbstractPolytope)

A hacky way of choosing a reasonable center for a polytope.
"""
LazySets.center(p::AbstractPolytope) = LazySets.center(overapproximate(p))
function LazySets.overapproximate(lazy_set::AbstractPolytope,sphere::H) where {H<:BallType}
    r = maximum(map(v->norm(v-get_center(sphere)), vertices_list(lazy_set)))
    Ball2(get_center(sphere),r)
end
function LazySets.overapproximate(s::AbstractPolytope,sphere::Type{H}) where {H<:BallType}
    overapproximate(s,H(LazySets.center(s),1.0))
end

export add_child_approximation!
function add_child_approximation!(g,model,parent_id,child_id)
    @assert has_vertex(g,parent_id)
    @assert !has_vertex(g,child_id)
    geom = overapproximate(get_node(g,parent_id).geom,model)
    add_node!(g,GeomNode(geom),child_id)
    add_edge!(g,parent_id,child_id)
    return g
end

export construct_geometry_tree!
function construct_geometry_tree!(g,geom)
    add_node!(g,GeomNode(geom),:BaseGeom)
    add_child_approximation!(g,equatorial_overapprox_model(),:BaseGeom,:Polyhedron)
    add_child_approximation!(g,Hyperrectangle,:Polyhedron,:Hyperrectangle)
    add_child_approximation!(g,Ball2,         :Polyhedron,:Hypersphere)
end

### Collision Table
const CollisionStack = Dict{Int,Set{ID}} where {ID}
get_active_collision_ids(c::CollisionStack,t::Int) = get(c,t,valtype(c)())

"""
    CollisionTable <: AbstractCustomGraph{Graph,CollisionStack,I}

A Data structure for efficient discrete-time collision checking between large
    numbers of objects.
An edge between two nodes indicates that those nodes are "eligible" to collide.
If there is no edge, we don't need to check collisions.
Each node of the graph is a `CollisionStack`--a kind of timer bank that keeps
track of when collision checking needs to be performed between two objects.
The timing depends on the distance between the objects at a given time step, as
well as the maximum feasible speed of each object.
"""
@with_kw struct CollisionTable{ID} <: AbstractCustomGraph{Graph,CollisionStack{ID},ID}
    graph               ::Graph                 = Graph()
    nodes               ::Vector{CollisionStack{ID}} = Vector{CollisionStack{ID}}()
    vtx_map             ::Dict{ID,Int}           = Dict{ID,Int}()
    vtx_ids             ::Vector{ID}             = Vector{ID}()
end

function get_transformed_config end
function get_max_speed end

"""
    update_collision_table!(table,env_state,env,i,t=0)

Updates the counters in a CollisionTable. Should be called each time the modeled
process "steps forward" in time.
"""
function update_collision_table!(table,env_state,env,i,t=0)
    c = get_tranformed_config(env,env_state,i,t)
    stack = get_node(table,i)
    active_ids = get_active_collision_ids(stack,t)
    while !isempty(active_ids)
        j = pop!(active_ids)
        has_edge(table,i,j) ? nothing : continue
        cj = get_transformed_config(env,env_state,j,t)
        d_min = distance_lower_bound(c,cj)
        v_max = get_max_speed(env,i) + get_max_speed(env,j)
        dt = d_min / v_max
        push!(get!(stack,t+Int(floor(dt)),valtype(stack)()),j)
    end
end

function init_collision_table!(table,env_state,env,i,t=0)
    stack = get_node(table,i)
    for j in outneighbors(table,i)
        push!(get!(stack,t,valtype(stack)()),get_vtx_id(table,j))
    end
    update_collision_table!(table,env_state,env,i,t)
end

"""
    find_collision(table,env_state,env,i,t=0)

Checks for collisions between object "i" and all objects that may be (according
    to the table) near enough to warrant collision checking at time `t`.
"""
function find_collision(table,env_state,env,i,t=0)
    c = get_tranformed_config(env,env_state,i,t)
    stack = get_node(table,i)
    active_ids = get_active_collision_ids(stack,t)
    for j in active_ids
        has_edge(table,i,j) ? nothing : continue
        cj = get_transformed_config(env,env_state,j,t)
        if has_overlap(c,cj)
            return true, j
        end
    end
    return false, -1
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

end
