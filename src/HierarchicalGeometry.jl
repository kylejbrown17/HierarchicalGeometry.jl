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
    identity_linear_map,
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

include("overapproximation.jl")

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
function transform(g::VPolytope,t::CoordinateTransformations.LinearMap{R}) where {R}#{R<:Rotation}
    VPolytope(map(v->transform(v,t), vertices_list(g)))
end
function transform(g::HPolytope,t::CoordinateTransformations.LinearMap{R}) where {R}#{R<:Rotation}
    HPolytope(map(h->transform(h,t), constraints_list(g)))
end
function transform(g::BaseGeometry,t::CoordinateTransformations.AffineMap{R,T}) where {R,T}#{R<:Rotation,T}
    transform(
        transform(g,CoordinateTransformations.LinearMap(t.linear)),
        CoordinateTransformations.Translation(t.translation)
        )
end
identity_linear_map3() = compose(CoordinateTransformations.Translation(zero(SVector{3,Float64})),CoordinateTransformations.LinearMap(one(SMatrix{3,3,Float64})))
identity_linear_map2() = compose(CoordinateTransformations.Translation(zero(SVector{3,Float64})),CoordinateTransformations.LinearMap(one(SMatrix{3,3,Float64})))
identity_linear_map() = identity_linear_map3()

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
    TransformNode(identity_linear_map(),identity_linear_map())
end
local_transform(n::TransformNode) = n.local_transform
global_transform(n::TransformNode) = n.global_transform
function set_local_transform!(n::TransformNode,t)
    n.local_transform = t
end
function set_global_transform!(n::TransformNode,t)
    n.global_transform = t
end
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
const transform_node_accessor_interface = [:local_transform,:global_transform]
const transform_node_mutator_interface = [:set_local_transform!,:set_global_transform!]

export
    # TransformTree,
    set_child!,
    get_parent_transform,
    update_transform_tree!

# @with_kw struct TransformTree{ID} <: AbstractCustomTree{TransformNode,ID}
#     graph               ::DiGraph               = DiGraph()
#     nodes               ::Vector{TransformNode} = Vector{TransformNode}()
#     vtx_map             ::Dict{ID,Int}          = Dict{ID,Int}()
#     vtx_ids             ::Vector{ID}            = Vector{ID}()
# end
for op in transform_node_accessor_interface
    @eval $op(g::AbstractCustomTree,v) = $op(get_node(g,v))
end
function get_parent_transform(g::AbstractCustomTree,v)
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
function update_transform_tree!(g::AbstractCustomTree,v)
    n = get_node(g,v)
    set_global_transform!(n,get_parent_transform(g,v) ∘ local_transform(n))
    for vp in outneighbors(g,v)
        update_transform_tree!(g,vp)
    end
    return g
end
function set_local_transform!(g::AbstractCustomTree,v,t,update=true)
    set_local_transform!(get_node(g,v),t)
    if update
        update_transform_tree!(g,v)
    end
    return g
end
"""
    set_child!(tree,parent,child,new_local_transform)

Replace existing edge `old_parent` → `child` with new edge `parent` → `child`.
Set the `child.local_transform` to `new_local_transform`.
"""
function set_child!(tree::AbstractCustomTree,parent,child,
        t=relative_transform(global_transform(tree,parent),
            global_transform(tree,child))
    )
    rem_edge!(tree,get_parent(tree,child),child)
    add_edge!(tree,parent,child)
    @assert !is_cyclic(tree) "adding edge $(parent) → $(child) made tree cyclic"
    set_local_transform!(tree,child,t)
end
const transform_tree_mutator_interface = [:set_local_transform!,:set_global_transform!]

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
"""
    Base.copy(e::CachedElement)

Shares the e.element, since it doesn't need to be replaced until `set_element!`
is called. Copies `e.is_up_to_date` to preserve the cache state.
"""
Base.copy(e::CachedElement) = CachedElement(e.element,copy(is_up_to_date(e)))

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
    # TODO implement a version of this that includes the Tree, not just the node
    if !is_up_to_date(n)
        transformed_geom = transform(get_base_geom(n),global_transform(n))
        update_element!(n,transformed_geom)
    end
    return get_element(n.cached_geom)
end
distance_lower_bound(a::GeomNode{G},b::GeomNode{G}) where {G<:Union{BallType,RectType}} = distance_lower_bound(get_cached_geom(a),get_cached_geom(b))
const geom_node_accessor_interface = [
    cached_element_accessor_interface...,
    transform_node_accessor_interface...,
    :get_base_geom, :get_cached_geom,
    ]
const geom_node_mutator_interface = [
    cached_element_mutator_interface...,
    transform_node_mutator_interface...
]
"""
    Base.copy(n::GeomNode)

Shares `n.base_geom`, deepcopies `n.transform_node`, and copies `n.cached_geom`
(see documenation for `Base.copy(::CachedElement)`).
"""
Base.copy(n::GeomNode) = GeomNode(
    n.base_geom,
    deepcopy(n.transform_node),
    copy(n.cached_geom)
)

export GeometryHierarchy
"""
    GeometryHierarchy

A hierarchical representation of geometry
Fields:
* graph - encodes the hierarchy of geometries
* nodes - geometry nodes
"""
@with_kw struct GeometryHierarchy <: AbstractCustomTree{GeomNode,Symbol}
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
"""
    TODO: The default behavior should dispatch on the geometry type. We may not
    care for all successive approximations with e.g., a large assembly.
"""
function construct_geometry_tree!(g,geom)
    add_node!(g,GeomNode(geom),:BaseGeom)
    add_child_approximation!(g,equatorial_overapprox_model(),:BaseGeom,:Polyhedron)
    add_child_approximation!(g,Hyperrectangle,:Polyhedron,:Hyperrectangle)
    add_child_approximation!(g,Ball2,         :Polyhedron,:Hypersphere)
end

export
    AssemblyID,
    TransportUnitID,
    SceneNode,
    RobotNode,
    ObjectNode,
    AssemblyNode,
        components,
        add_component!,
    TransportUnitNode,
        robot_team,
        add_robot!,
    SceneTree

const TransformDict{T} = Dict{T,CoordinateTransformations.Transformation}
"""
    abstract type SceneNode end

An Abstract type, of which all nodes in a SceneTree are concrete subtypes.
"""
@with_kw struct AssemblyID <: AbstractID
    id::Int = -1
end
@with_kw struct TransportUnitID <: AbstractID
    id::Int = -1
end
# struct BaseNode <: SceneNode
#     id::VtxID
# end
abstract type SceneNode end
for op in geom_node_accessor_interface
    @eval $op(n::SceneNode) = $op(n.geom)
end
for op in geom_node_mutator_interface
    @eval $op(n::SceneNode,val) = $op(n.geom,val)
end
"""
    Base.copy(n::N) where {N<:SceneNode}

Shares `n.base_geom`, deepcopies `n.transform_node`, and copies `n.cached_geom`
(see documenation for `Base.copy(::CachedElement)`).
"""
Base.copy(n::N) where {N<:SceneNode} = N(n,copy(n.geom))

struct RobotNode{R} <: SceneNode
    id::BotID{R}
    geom::GeomNode
end
RobotNode{R}(n::RobotNode,geom::GeomNode) where {R} = RobotNode(n.id,geom)
RobotNode(n::RobotNode,geom::GeomNode) = RobotNode(n.id,geom)
struct ObjectNode <: SceneNode
    id::ObjectID
    geom::GeomNode
end
ObjectNode(n::ObjectNode,geom::GeomNode) = ObjectNode(n.id,geom)
has_component(n::SceneNode,id) = false
struct AssemblyNode <: SceneNode
    id::AssemblyID
    geom::GeomNode
    components::TransformDict{Union{ObjectID,AssemblyID}}
end
AssemblyNode(n::AssemblyNode,geom::GeomNode) = AssemblyNode(n.id,geom,n.components)
AssemblyNode(id,geom) = AssemblyNode(id,geom,TransformDict{Union{ObjectID,AssemblyID}}())
components(n::AssemblyNode) = n.components
add_component!(n::AssemblyNode,p) = push!(n.components,p)
has_component(n::AssemblyNode,id) = haskey(components(n),id)
child_transform(n::AssemblyNode,id) = components(n)[id]
struct TransportUnitNode <: SceneNode
    id::TransportUnitID
    geom::GeomNode
    assembly::Pair{AssemblyID,CoordinateTransformations.Transformation}
    robots::TransformDict{BotID}
end
TransportUnitNode(n::TransportUnitNode,geom::GeomNode) = TransportUnitNode(n.id,geom,n.assembly,n.robots)
TransportUnitNode(id,geom,assembly) = TransportUnitNode(id,geom,assembly,TransformDict{BotID}())
TransportUnitNode(id,geom,assembly_id::AssemblyID) = TransportUnitNode(id,geom,assembly_id=>identity_linear_map())
robot_team(n::TransportUnitNode) = n.robots
assembly_id(n::TransportUnitNode) = n.assembly.first
has_component(n::TransportUnitNode,id::AssemblyID) = id == assembly_id(n)
has_component(n::TransportUnitNode,id::BotID) = haskey(robot_team(n),id)
child_transform(n::TransportUnitNode,id::AssemblyID) = n.assembly.second
child_transform(n::TransportUnitNode,id::BotID) = robot_team(n)[id]
add_robot!(n::TransportUnitNode,p) = push!(robot_team(n),p)
add_robot!(n::TransportUnitNode,r,t) = add_robot!(n,r=>t)
"""
    SceneTree

A tree data structure for describing the state of a manufacturing project.
`AssemblyNode`s define how objects/subassemblies fit together, and
`TransportUnitNode` define how robots fit together to form a transport team.
"""
@with_kw struct SceneTree <: AbstractCustomTree{SceneNode,AbstractID}
    graph               ::DiGraph               = DiGraph()
    nodes               ::Vector{SceneNode}     = Vector{SceneNode}()
    vtx_map             ::Dict{AbstractID,Int}  = Dict{AbstractID,Int}()
    vtx_ids             ::Vector{AbstractID}    = Vector{AbstractID}()
end
GraphUtils.add_node!(tree::SceneTree,node::SceneNode) = add_node!(tree,node,node.id)
GraphUtils.get_vtx(tree::SceneTree,n::SceneNode) = get_vtx(tree,n.id)
function Base.copy(tree::SceneTree)
    SceneTree(
        graph = deepcopy(tree.graph),
        nodes = map(copy, tree.nodes),
        vtx_map = deepcopy(tree.vtx_map), # TODO would a shallow copy be fine?
        vtx_ids = deepcopy(tree.vtx_ids) # TODO would a shallow copy be fine?
    )
end

set_child!(tree::SceneTree,parent::SceneNode,args...) = set_child!(tree,parent.id,args...)
set_child!(tree::SceneTree,parent::AbstractID,child::SceneNode,args...) = set_child!(tree,parent,child.id,args...)
"""
    set_child!(tree::SceneTree,parent::AbstractID,child::AbstractID)

Only eligible children can be added to a candidate parent in a `SceneTree`.
Eligible edges are those for which the child node is listed in the components of
the parent. `ObjectNode`s and `RobotNode`s may not have any children.
"""
function set_child!(tree::SceneTree,parent::AbstractID,child::AbstractID)
    node = get_node(tree,parent)
    @assert has_component(node,child) "$(get_node(tree,child)) cannot be a child of $(node)"
    t = child_transform(node,child)
    set_child!(tree,parent,get_vtx(tree,child),t)
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


end
