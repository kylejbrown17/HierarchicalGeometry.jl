module HierarchicalGeometry

using StaticArrays
using LinearAlgebra
using Polyhedra
using LazySets
using GeometryBasics
# using RigidBodyDynamics
using CoordinateTransformations
using Rotations
using LightGraphs, GraphUtils
using Parameters
using Logging

# TODO convert RigidBodyDynamics.CartesianFrame3D to CoordinateTransformations.Transformation
# TODO add RigidBodyDynamics.Mechanism(s) to SceneTree
# TODO replace TransformNode with RigidBodyDynamics equivalent. Each SceneNode should have a RigidBody included in it.

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
transform(v,t) = t(v)
# transform(v::V,t) where {V<:AbstractVector} = V(t(v))

(t::CoordinateTransformations.LinearMap)(::Nothing) = nothing
(t::CoordinateTransformations.Translation)(::Nothing) = nothing
(t::CoordinateTransformations.LinearMap)(h::LazySets.HalfSpace) = LazySets.HalfSpace(t(Vector(h.a)),h.b)
(t::CoordinateTransformations.Translation)(g::BaseGeometry) = LazySets.translate(g,Vector(t.translation))
(t::CoordinateTransformations.LinearMap)(g::Ball2) = Ball2(t(g.center),g.radius)
"""
    (t::CoordinateTransformations.LinearMap)(g::Hyperrectangle)

"rotating" a Hyperrectangle `g` results in a new Hyperrectangle that bounds the
transformed version `g`.
"""
(t::CoordinateTransformations.LinearMap)(g::Hyperrectangle) = overapproximate(t(convert(LazySets.VPolytope,g)),Hyperrectangle)
(t::CoordinateTransformations.LinearMap)(g::VPolytope) = VPolytope(map(t, vertices_list(g)))
(t::CoordinateTransformations.LinearMap)(g::HPolytope) = HPolytope(map(t, constraints_list(g)))

identity_linear_map3() = compose(CoordinateTransformations.Translation(zero(SVector{3,Float64})),CoordinateTransformations.LinearMap(one(SMatrix{3,3,Float64})))
identity_linear_map2() = compose(CoordinateTransformations.Translation(zero(SVector{3,Float64})),CoordinateTransformations.LinearMap(one(SMatrix{3,3,Float64})))
identity_linear_map() = identity_linear_map3()
Base.convert(::Type{Hyperrectangle{Float64,T,T}},rect::Hyperrectangle) where {T} = Hyperrectangle(T(rect.center),T(rect.radius))
# function transform!(v::V,t) where {V<:AbstractVector}
#     v .= transform(v,t)
# end
# function transform!(h::LazySets.HalfSpace,t)
#     h.a .= transform(h,t).a
#     return h
# end
# transform!(g::BaseGeometry,t::CoordinateTransformations.Translation) = LazySets.translate!(g,t.translation)
# function transform!(g::Ball2,t::CoordinateTransformations.LinearMap)
#     transform!(g.center,t)
#     return g
# end

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
for op in (:relative_transform,:(Rotations.rotation_error))
    @eval $op(tree::AbstractCustomTree,parent,child,args...) = $op(
        global_transform(tree,parent),
        global_transform(tree,child),
        args...
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

const cached_element_accessor_interface = [
    :(GraphUtils.is_up_to_date),
    :(GraphUtils.time_stamp)
    ]
const cached_element_mutator_interface = [
    :(GraphUtils.update_element!),
    :(GraphUtils.set_time_stamp!),
    :(GraphUtils.set_element!),
    :(GraphUtils.set_up_to_date!)
    ]


export
    TransformNode,
    local_transform,
    global_transform,
    set_local_transform!,
    set_global_transform!

@with_kw struct TransformNodeID <: AbstractID
    id::Int = -1
end

"""
    TransformNode

Has a parent field that points to another TransformNode.
"""
mutable struct TransformNode <: CachedTreeNode{TransformNodeID}
    id::TransformNodeID
    local_transform::CoordinateTransformations.Transformation # transform from the parent frame to the child frame
    global_transform::CachedElement{CoordinateTransformations.Transformation}
    parent::TransformNode # parent node
    children::Dict{AbstractID,CachedTreeNode}
    function TransformNode(
            a::CoordinateTransformations.Transformation,
            b::CachedElement)
        t = new()
        t.id = get_unique_id(TransformNodeID)
        t.local_transform = a
        t.global_transform = b
        t.parent = t
        # t.children = Dict{TransformNodeID,TransformNodeID}()
        t.children = Dict{AbstractID,CachedTreeNode}()
        return t
    end
end
function TransformNode()
    TransformNode(identity_linear_map(),CachedElement(identity_linear_map()))
end
function TransformNode(
    a::CoordinateTransformations.Transformation,
    b::CoordinateTransformations.Transformation
    )
    TransformNode(a,CachedElement(b))
end
# TODO copying as currently implemented will cause major issues with 
function Base.copy(n::TransformNode)
    node = TransformNode(
    deepcopy(n.local_transform),
    deepcopy(n.global_transform))
    # Don't want to copy parent or children. Just need to reconnect later.
    return node
end
# Cached Tree interface
GraphUtils.cached_element(n::TransformNode) = n.global_transform
tf_up_to_date(n::TransformNode) = GraphUtils.cached_node_up_to_date(n)
set_tf_up_to_date!(n::TransformNode,val=true) = GraphUtils.set_cached_node_up_to_date!(n,val)
local_transform(n::TransformNode) = n.local_transform
set_global_transform!(n::TransformNode,t,args...) = GraphUtils.update_element!(n,t,args...)
global_transform(n::TransformNode) = GraphUtils.get_cached_value!(n)
function GraphUtils.propagate_forward!(parent::TransformNode,child::TransformNode) 
    if parent === child
        set_global_transform!(child,local_transform(child))
    else
        set_global_transform!(child,global_transform(parent) ∘ local_transform(child))
    end
end
function set_local_transform!(n::TransformNode,t,update=false)
    n.local_transform = t
    set_tf_up_to_date!(n,false)
    if update
        global_transform(n) # just get the global transform
    end
    return n.local_transform
end
const transform_node_accessor_interface = [:tf_up_to_date,:local_transform,:global_transform]
const transform_node_mutator_interface = [:set_tf_up_to_date!,:set_local_transform!,:set_global_transform!]

for op in transform_node_accessor_interface
    @eval $op(tree::AbstractCustomTree,v) = $op(get_node(tree,v))
    @eval $op(tree::AbstractCustomTree,n::TransformNode) = $op(n)
end
for op in transform_node_mutator_interface
    @eval $op(tree::AbstractCustomTree,v,args...) = $op(get_node(tree,v),args...)
    @eval $op(tree::AbstractCustomTree,n::TransformNode,args...) = $op(n,args...)
end

export
    set_child!,
    get_parent_transform,
    update_transform_tree!

"""
    set_child!(tree,parent,child,new_local_transform)

Replace existing edge `old_parent` → `child` with new edge `parent` → `child`.
Set the `child.local_transform` to `new_local_transform`.
"""
function set_child!(tree::AbstractCustomTree,parent,child,
        t=relative_transform(tree,parent,child),
        edge=nothing
    )
    rem_edge!(tree,get_parent(tree,child),child)
    if add_edge!(tree,parent,child,edge)
        set_parent!(get_node(tree,child),get_node(tree,parent)) # for TransformNode
        @assert !is_cyclic(tree) "adding edge $(parent) → $(child) made tree cyclic"
        set_local_transform!(tree,child,t)
        return true
    end
    return false
end
const transform_tree_mutator_interface = [:set_local_transform!,:set_global_transform!]
function GraphUtils.add_child!(tree::AbstractCustomTree,parent,child,child_id)
    n = add_node!(tree,child,child_id)
    if set_child!(tree,get_node(tree,parent),child_id)
        return n
    else
        rem_node!(tree,child_id)
        return nothing
    end
end

export
    GeomNode,
    get_base_geom,
    get_cached_geom

@with_kw struct GeomID <: AbstractID
    id::Int = -1
end

mutable struct GeomNode{G} <: CachedTreeNode{GeomID}
    id::GeomID
    base_geom::G
    parent::TransformNode 
    cached_geom::CachedElement{G}
end
function GeomNode(geom)
    GeomNode(get_unique_id(GeomID),geom,TransformNode(),CachedElement(geom))
end
function GeomNode(geom,tf)
    GeomNode(get_unique_id(GeomID),geom,tf,CachedElement(geom))
end
GraphUtils.get_children(n::GeomNode) = Dict{AbstractID,CachedTreeNode}()
GraphUtils.cached_element(n::GeomNode) = n.cached_geom
set_cached_geom!(n::GeomNode,geom) = update_element!(n,geom)

get_base_geom(n::GeomNode) = n.base_geom
get_transform_node(n::GeomNode) = n.parent
for op in transform_node_accessor_interface
    @eval $op(g::GeomNode,args...) = $op(get_transform_node(g))
end
for op in transform_node_mutator_interface
    @eval $op(g::GeomNode,args...) = $op(get_transform_node(g),args...)
end
function GraphUtils.set_parent!(a::GeomNode,b::GeomNode)
    set_parent!(get_transform_node(a),get_transform_node(b))
end
function GraphUtils.propagate_forward!(t::TransformNode,n::GeomNode)
    transformed_geom = transform(get_base_geom(n),global_transform(n))
    set_cached_geom!(n,transformed_geom)
end


"""
    get_cached_geom(n::GeomNode)

If `n.cached_geom` is out of date, transform it according to
`global_transform(n)`. Updates both the global transform and geometry if 
necessary.
"""
get_cached_geom(n::GeomNode) = GraphUtils.get_cached_value!(n)
get_cached_geom(tree::AbstractCustomTree,v) = get_cached_geom(get_node(tree,v))
distance_lower_bound(a::GeomNode{G},b::GeomNode{G}) where {G<:Union{BallType,RectType}} = distance_lower_bound(get_cached_geom(a),get_cached_geom(b))
const geom_node_accessor_interface = [
    transform_node_accessor_interface...,
    :get_base_geom, :get_cached_geom,
    ]
const geom_node_mutator_interface = [
    transform_node_mutator_interface...
]

"""
    Base.copy(n::GeomNode)

Shares `n.base_geom`, deepcopies `n.transform_node`, and copies `n.cached_geom`
(see documenation for `Base.copy(::CachedElement)`).
"""
Base.copy(n::GeomNode) = GeomNode(
    get_unique_id(GeomID),
    n.base_geom,
    copy(get_transform_node(n)),
    copy(n.cached_geom)
)


export GeometryHierarchy

abstract type GeometryKey end
""" 
    BaseGeomKey <: GeometryKey

Points to any kind of geometry.
"""
struct BaseGeomKey <: GeometryKey end
"""
    PolyhedronKey <: GeometryKey

Points to a polyhedron.
"""
struct PolyhedronKey <: GeometryKey end
"""
    ZonotopeKey<: GeometryKey

Points to a zonotope. Might not be useful for approximating shapes that are very
asymmetrical.
"""
struct ZonotopeKey <: GeometryKey end
struct HyperrectangleKey <: GeometryKey end
struct HypersphereKey <: GeometryKey end
struct CylinderKey <: GeometryKey end

construct_child_approximation(::PolyhedronKey,geom)     = LazySets.overapproximate(geom,equatorial_overapprox_model())
construct_child_approximation(::HypersphereKey,geom)    = LazySets.overapproximate(geom,Ball2)
construct_child_approximation(::HyperrectangleKey,geom) = LazySets.overapproximate(geom,Hyperrectangle)

"""
    GeometryHierarchy

A hierarchical representation of geometry
Fields:
* graph - encodes the hierarchy of geometries
* nodes - geometry nodes
"""
@with_kw struct GeometryHierarchy <: AbstractCustomNTree{GeomNode,GeometryKey}
    graph               ::DiGraph               = DiGraph()
    nodes               ::Vector{GeomNode}      = Vector{GeomNode}()
    vtx_map             ::Dict{GeometryKey,Int} = Dict{GeometryKey,Int}()
    vtx_ids             ::Vector{GeometryKey}   = Vector{GeometryKey}() # maps vertex uid to actual graph node
end
# TODO How to safely copy? All nodes have the same parent.
function geom_hierarchy(geom::GeomNode)
    h = GeometryHierarchy()
    add_node!(h,geom,BaseGeomKey())
    return h
end

distance_lower_bound(a::GeometryHierarchy,b::GeometryHierarchy) = distance_lower_bound(
    get_node(a,HypersphereKey()),
    get_node(b,HypersphereKey())
    )
function has_overlap(a::GeometryHierarchy,b::GeometryHierarchy,leaf_id=HypersphereKey())
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
function add_child_approximation!(g::GeometryHierarchy,parent_id,child_id)
    @assert has_vertex(g,parent_id)
    @assert !has_vertex(g,child_id)
    node = get_node(g,parent_id)
    # geom = overapproximate(get_base_geom(node),model)
    geom = construct_child_approximation(child_id,get_base_geom(node))
    add_node!(g,
        GeomNode(geom,node.parent), # Share parent
        child_id
        ) # todo needs parent
    add_edge!(g,parent_id,child_id)
    return g
end


export construct_geometry_tree!
"""

TODO: Ensure that the parent member of the base geometry key is the parent 
of the SceneNode to which the Geometry Hierarchy is bound.
"""
function construct_geometry_tree!(g::GeometryHierarchy,geom::GeomNode)
    add_node!(g,geom,BaseGeomKey())
    add_child_approximation!(g,BaseGeomKey(),PolyhedronKey())
    add_child_approximation!(g,PolyhedronKey(),HyperrectangleKey())
    add_child_approximation!(g,PolyhedronKey(),HypersphereKey())
end
construct_geometry_tree!(g::GeometryHierarchy,geom) = construct_geometry_tree!(g,GeomNode(geom))

get_cached_geom(n::GeometryHierarchy) = get_cached_geom(n,BaseGeomKey())


export
    AssemblyID,
    TransportUnitID,
    SceneNode,
    RobotNode,
    ObjectNode,
    AssemblyNode,
        components,
        add_component!,
        child_transform,
    TransportUnitNode,
        robot_team,
        add_robot!,
    SceneTree,
        capture_child!,
        disband!

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

"""
    abstract type SceneNode

A node of the `SceneTree`. Each concrete subtype of `SceneNode` 
contains all of the following information:
- All required children (whether for a temporary or permanent edge)
- Required Parent (For a permanent edge only)
- Geometry of the entity represented by the node
- Unique ID accessible via `GraphUtils.node_id(node)`
- Current transform--both global (relative to base frame) and local (relative to 
current parent)
- Required transform relative to parent (if applicable)

"""
abstract type SceneNode end
abstract type SceneAssemblyNode <: SceneNode end
# SceneNode interface
GraphUtils.node_id(n::SceneNode) = n.id

for op in geom_node_accessor_interface
    @eval $op(n::SceneNode) = $op(n.geom)
    @eval $op(n::CustomNode) = $op(node_val(n))
end
for op in geom_node_mutator_interface
    @eval $op(n::SceneNode,args...) = $op(n.geom,args...)
    @eval $op(n::CustomNode,args...) = $op(node_val(n),args...)
end
GraphUtils.set_parent!(a::SceneNode,b::SceneNode) = set_parent!(a.geom,b.geom)
GraphUtils.set_parent!(a::CustomNode,b::CustomNode) = set_parent!(node_val(a),node_val(b))

"""
    required_transforms_to_children(n::SceneNode)

Returns a dictionary mapping child id to required relative transform.
"""
function required_transforms_to_children(n::SceneNode)
    return TransformDict{AbstractID}()
end

"""
    required_parent(n::SceneNode)

Returns the id of the required parent node. This is only applicable to 
`ObjectNode` and `AssemblyNode`.
"""
function required_parent end

"""
    required_transform_to_parent(n::SceneNode,parent_id)

Returns the required transform relative to the parent.
"""
function required_transform_to_parent end

"""
    Base.copy(n::N) where {N<:SceneNode}

Shares `n.base_geom`, deepcopies `n.transform_node`, and copies `n.cached_geom`
(see documenation for `Base.copy(::CachedElement)`).
"""
Base.copy(n::N) where {N<:SceneNode} = N(n,copy(n.geom))

struct RobotNode{R} <: SceneNode
    id::BotID{R}
    geom::GeomNode
    geom_hierarchy::GeometryHierarchy
end
RobotNode(id::BotID,geom) = RobotNode(id,geom,geom_hierarchy(geom))
RobotNode(n::RobotNode,geom) = RobotNode(n.id,geom)
struct ObjectNode <: SceneNode
    id::ObjectID
    geom::GeomNode
    geom_hierarchy::GeometryHierarchy
end
ObjectNode(id::ObjectID,geom) = ObjectNode(id,geom,geom_hierarchy(geom))
ObjectNode(n::ObjectNode,geom) = ObjectNode(n.id,geom)
has_component(n::SceneNode,id) = false
struct AssemblyNode <: SceneNode
    id::AssemblyID
    geom::GeomNode
    components::TransformDict{Union{ObjectID,AssemblyID}}
    geom_hierarchy::GeometryHierarchy
end
AssemblyNode(n::AssemblyNode,geom) = AssemblyNode(n.id,geom,n.components,geom_hierarchy(geom))
AssemblyNode(id,geom) = AssemblyNode(id,geom,TransformDict{Union{ObjectID,AssemblyID}}(),geom_hierarchy(geom))
components(n::AssemblyNode)         = n.components
add_component!(n::AssemblyNode,p)   = push!(n.components,p)
has_component(n::AssemblyNode,id)   = haskey(components(n),id)
child_transform(n::AssemblyNode,id) = components(n)[id]

struct TransportUnitNode <: SceneNode
    id::TransportUnitID
    geom::GeomNode
    assembly::Pair{AssemblyID,CoordinateTransformations.Transformation}
    robots::TransformDict{BotID} # must be filled with unique invalid ids
    geom_hierarchy::GeometryHierarchy
end
# TransportUnitNode(id::TransportUnitID,a::Pair,args...) = TransportUnitNode(id,
#     convert(Pair{AssemblyID,CoordinateTransformations.Transformation},a),
#     args...)
TransportUnitNode(n::TransportUnitNode,geom) = TransportUnitNode(n.id,geom,n.assembly,n.robots,geom_hierarchy(geom))
TransportUnitNode(id,geom,assembly) = TransportUnitNode(id,geom,assembly,TransformDict{BotID}(),geom_hierarchy(geom))
TransportUnitNode(id,geom,assembly_id::AssemblyID) = TransportUnitNode(id,geom,assembly_id=>identity_linear_map())
robot_team(n::TransportUnitNode) = n.robots
assembly_id(n::TransportUnitNode) = n.assembly.first
Base.convert(::Pair{A,B},pair) where {A,B} = convert(A,pair.first)=>convert(B,p.second) 

has_component(n::TransportUnitNode,id::AssemblyID)  = id == assembly_id(n)
has_component(n::TransportUnitNode,id::BotID)       = haskey(robot_team(n),id)
child_transform(n::TransportUnitNode,id::AssemblyID) = n.assembly.second
child_transform(n::TransportUnitNode,id::BotID)     = robot_team(n)[id]
add_robot!(n::TransportUnitNode,p)                  = push!(robot_team(n),p)
add_robot!(n::TransportUnitNode,r,t)                = add_robot!(n,r=>t)

# Necessary for copying
RobotNode{R}(n::RobotNode,args...) where {R} = RobotNode(n.id,args...)
# for T in (:ObjectNode,:AssemblyNode,:TransportUnitNode)
#     @eval $T{R}(n::$T,args...) where {R} = $T(n,args...)
# end

"""
    abstract type SceneTreeEdge

Concrete subtypes are `PermanentEdge` or `TemporaryEdge`. The idea is that 
the former is not allowed to be removed, whereas the latter stays in place.
"""
abstract type SceneTreeEdge end
struct PermanentEdge <: SceneTreeEdge end
struct TemporaryEdge <: SceneTreeEdge end
for U in (:TransportUnitID,:TransportUnitNode)
    for V in (:RobotID,:RobotNode,:AssemblyID,:AssemblyNode)
        @eval GraphUtils.make_edge(g,u::$U,v::$V) = TemporaryEdge()
    end
end
for U in (:AssemblyID,:AssemblyNode)
    for V in (:ObjectID,:ObjectNode,:AssemblyID,:AssemblyNode)
        @eval GraphUtils.make_edge(g,u::$U,v::$V) = PermanentEdge()
    end
end

"""
    SceneTree <: AbstractCustomNETree{SceneNode,SceneTreeEdge,AbstractID}

A tree data structure for describing the current configuration of a multi-entity 
system, as well as the desired final configuration. The configuration of the 
system is defined by the topology of the tree (i.e., which nodes are connected
by edges) and the configuration of each node (its transform relative to its 
parent node).
There are four different concrete subtypes of `SceneNode`
- `ObjectNode` refers to a single inanimate rigid body
- `RobotNode` refers to a single rigid body robot
- `AssemblyNode` refers to a rigid collection of multiple objects and/or 
subassemblies
- `TransportUnitNode` refers to a team of one or more robot that fit together as 
a transport team to move an assembly or object.
The edges of the SceneTree are either `PermanentEdge` or `TemporaryEdge`.
Temporary edges (removed after transport is complete)
- TransportUnitNode → RobotNode
- TransportUnitNode → ObjectNode
- TransportUnitNode → AssemblyNode
Permanent edges (cannot be broken after placement)
- AssemblyNode → AssemblyNode
- AssemblyNode → ObjectNode
"""
@with_kw struct SceneTree <: AbstractCustomNETree{SceneNode,SceneTreeEdge,AbstractID}
    graph               ::DiGraph               = DiGraph()
    nodes               ::Vector{SceneNode}     = Vector{SceneNode}()
    vtx_map             ::Dict{AbstractID,Int}  = Dict{AbstractID,Int}()
    vtx_ids             ::Vector{AbstractID}    = Vector{AbstractID}()
    inedges             ::Vector{Dict{Int,SceneTreeEdge}}   = Vector{Dict{Int,SceneTreeEdge}}()
    outedges            ::Vector{Dict{Int,SceneTreeEdge}}   = Vector{Dict{Int,SceneTreeEdge}}()
end
GraphUtils.add_node!(tree::SceneTree,node::SceneNode) = add_node!(tree,node,node.id)
GraphUtils.get_vtx(tree::SceneTree,n::SceneNode) = get_vtx(tree,n.id)
function Base.copy(tree::SceneTree)
    SceneTree(
        graph = deepcopy(tree.graph),
        nodes = map(copy, tree.nodes), # TODO This may cause problems with TransformNode
        vtx_map = deepcopy(tree.vtx_map),
        vtx_ids = deepcopy(tree.vtx_ids)
    )
end

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
    set_child!(tree,parent,get_vtx(tree,child),t,make_edge(tree,parent,child))
end
set_child!(tree::SceneTree,parent::SceneNode,args...) = set_child!(tree,node_id(parent),args...)
set_child!(tree::SceneTree,parent::AbstractID,child::SceneNode,args...) = set_child!(tree,parent,node_id(child),args...)

"""
    LightGraphs.rem_edge!(tree::SceneTree,u,v)

Edge may only be removed if it is not permanent.
"""
function LightGraphs.rem_edge!(tree::SceneTree,u,v)
    if !has_edge(tree,u,v)
        return true
    end
    @assert !isa(get_edge(tree,u,v),PermanentEdge) "Edge $u → $v is permanent!"
    GraphUtils.delete_edge!(tree,u,v)
end

"""
    disband!(tree::SceneTree,n::TransportUnitNode)

Disband a transport unit by removing the edge to its children.
"""
function disband!(tree::SceneTree,n::TransportUnitNode)
    for (r,tform) in robot_team(n)
        if !has_edge(tree,n,r)
            @warn "Disbanding $n -- $r should be attached, but is not"
        end
        rem_edge!(tree,n,r)
    end
    rem_edge!(tree,n,assembly_id(n))
    return true
end

"""
    capture_child!(tree::SceneTree,u,v,ttol,rtol)

If the transform error of `v` relative to the transform prescribed for it by `u`
is small, allow `u` to capture `v` (`v` becomes a child of `u`).
"""
function capture_child!(tree::SceneTree,u,v,ttol=1e-2,rtol=1e-2)
    nu = get_node(tree,u)
    nv = get_node(tree,v)
    @assert has_component(nu,node_id(nv)) "$nu cannot capture $nv"
    t = relative_transform(tree,u,v)
    t_des = child_transform(nu,node_id(nv))
    et = norm(t.translation - t_des.translation) # translation error
    er = norm(Rotations.rotation_error(t,t_des)) # rotation_error
    if et < ttol && er < rtol
        if !is_root_node(tree,v)
            p = get_node(tree,get_parent(tree,v)) # current parent
            @assert(isa(p, TransportUnitNode),
                "Trying to capture child $v from non-TransportUnit parent $p")
            # NOTE that there may be a time gap if the assembly has to be lifted
            # into place by a separate "robot"
            disband!(tree,p)
        end
        return set_child!(tree,u,v)
    end
    return false
end


### Collision Table
const CollisionStack = Dict{Int,Set{ID}} where {ID}
get_active_collision_ids(c::CollisionStack,t::Int) = get(c,t,valtype(c)())

"""
    CollisionTable <: AbstractCustomNGraph{Graph,CollisionStack,I}

A Data structure for efficient discrete-time collision checking between large
    numbers of objects.
An edge between two nodes indicates that those nodes are "eligible" to collide.
If there is no edge, we don't need to check collisions.
Each node of the graph is a `CollisionStack`--a kind of timer bank that keeps
track of when collision checking needs to be performed between two objects.
The timing depends on the distance between the objects at a given time step, as
well as the maximum feasible speed of each object.
"""
@with_kw struct CollisionTable{ID} <: AbstractCustomNGraph{Graph,CollisionStack{ID},ID}
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
