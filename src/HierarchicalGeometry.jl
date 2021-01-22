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
const cached_element_accessor_interface = [:is_up_to_date]
const cached_element_mutator_interface = [:update_element!,:set_element!,:set_up_to_date!]
"""
    Base.copy(e::CachedElement)

Shares the e.element, since it doesn't need to be replaced until `set_element!`
is called. Copies `e.is_up_to_date` to preserve the cache state.
"""
Base.copy(e::CachedElement) = CachedElement(e.element,copy(is_up_to_date(e)))
Base.convert(::Type{CachedElement{T}},e::CachedElement) where {T} = CachedElement{T}(get_element(e),is_up_to_date(e))


export
    TransformNode,
    local_transform,
    global_transform,
    set_local_transform!,
    set_global_transform!

mutable struct TransformNode
    local_transform::CoordinateTransformations.Transformation # transform from the parent frame to the child frame
    global_transform::CachedElement{CoordinateTransformations.Transformation}
end
tf_up_to_date(n::TransformNode) = is_up_to_date(n.global_transform)
set_tf_up_to_date!(n::TransformNode,val) = set_up_to_date!(n.global_transform,val)
function TransformNode()
    TransformNode(identity_linear_map(),CachedElement(identity_linear_map()))
end
function TransformNode(
    a::CoordinateTransformations.Transformation,
    b::CoordinateTransformations.Transformation
    )
    TransformNode(a,CachedElement(b))
end
local_transform(n::TransformNode) = n.local_transform
function global_transform(n::TransformNode)
    if !is_up_to_date(n.global_transform)
        @warn string("global_transform out of date!",
        " Update with set_global_transform!(n,tf)")
    end
    return get_element(n.global_transform)
end
function set_local_transform!(n::TransformNode,t)
    n.local_transform = t
end
function set_global_transform!(n::TransformNode,t)
    update_element!(n.global_transform,t)
    return t
end
const transform_node_accessor_interface = [:tf_up_to_date,:local_transform,:global_transform]
const transform_node_mutator_interface = [:set_tf_up_to_date!,:set_local_transform!,:set_global_transform!]

export
    # TransformTree,
    set_child!,
    get_parent_transform,
    update_transform_tree!

# @with_kw struct TransformTree{ID} <: AbstractCustomNTree{TransformNode,ID}
#     graph               ::DiGraph               = DiGraph()
#     nodes               ::Vector{TransformNode} = Vector{TransformNode}()
#     vtx_map             ::Dict{ID,Int}          = Dict{ID,Int}()
#     vtx_ids             ::Vector{ID}            = Vector{ID}()
# end
tf_up_to_date(tree,v) = tf_up_to_date(get_node(tree,n))
local_transform(tree,v) = local_transform(get_node(tree,n))
function get_parent_transform(g::AbstractCustomTree,v)
    vp = get_parent(g,v)
    if has_vertex(g,vp)
        return global_transform(g,vp)
    end
    # @warn string("Node $v has no parent.",
    # " Returning identity_linear_map() as parent transform")
    return identity_linear_map()
end
"""
    global_transform(tree,v)

Return global transform of node `v` in `tree`. If the current global transform
is out of date, backtrack until `v` and all of its ancestors are up to date.
"""
function global_transform(tree,v)
    n = get_node(tree,v)
    if !tf_up_to_date(n)
        set_global_transform!(n,get_parent_transform(tree,v) ∘ local_transform(n))
    end
    return global_transform(n)
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
"""
    propagate_tf_flag!(g::AbstractCustomTree,v,val=false)

Calls `set_tf_up_to_date!(get_node(g,vp),val)` on v and all outneighbors thereof.
Needs to be called by `set_local_transform!(g, ...)` unless
`update_transform_tree!(g, ...)` is called instead.
"""
function propagate_tf_flag!(g::AbstractCustomTree,v,val=false)
    n = get_node(g,v)
    if tf_up_to_date(n)
        return g
    end
    set_tf_up_to_date!(n,val) # Mark as
    for vp in outneighbors(g,v)
        propagate_tf_flag!(g,vp,val)
    end
    return g
end
function set_local_transform!(g::AbstractCustomTree,v,t,update=true)
    set_local_transform!(get_node(g,v),t)
    if update
        update_transform_tree!(g,v)
    else
        propagate_tf_flag!(g,v,false)
    end
    return g
end
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
    add_edge!(tree,parent,child,edge)
    @assert !is_cyclic(tree) "adding edge $(parent) → $(child) made tree cyclic"
    set_local_transform!(tree,child,t)
end
const transform_tree_mutator_interface = [:set_local_transform!,:set_global_transform!]

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
set_tf_up_to_date!(n::GeomNode,val) = set_tf_up_to_date!(n.transform_node,val)
function set_local_transform!(n::GeomNode,t)
    set_up_to_date!(n,false)
    set_local_transform!(n.transform_node,t)
end
function set_global_transform!(n::GeomNode,t)
    set_up_to_date!(n,false)
    set_global_transform!(n.transform_node,t)
end
"""
    get_cached_geom(n::GeomNode)

If `n.cached_geom` is out of date, transform it according to
`global_transform(n)`. Note that nothing can be done if `global_transform(n)` is
out of date, since this function call doesn't include the tree.
"""
function get_cached_geom(n::GeomNode)
    if !is_up_to_date(n)
        transformed_geom = transform(get_base_geom(n),global_transform(n))
        update_element!(n,transformed_geom)
    end
    return get_element(n.cached_geom)
end
"""
    get_cached_geom(tree,v)

returns the cached geometry associated with vertex `v`, updating both the
global transform and geometry if necessary.
"""
function get_cached_geom(tree,v)
    n = get_node(tree,v)
    if !is_up_to_date(n)
        transformed_geom = transform(get_base_geom(n),global_transform(tree,v))
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
@with_kw struct GeometryHierarchy <: AbstractCustomNTree{GeomNode,Symbol}
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
    SceneTree,
    capture_child!

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

for op in geom_node_accessor_interface
    @eval $op(n::SceneNode) = $op(n.geom)
    @eval $op(n::CustomNode) = $op(node_val(n))
end
for op in geom_node_mutator_interface
    @eval $op(n::SceneNode,val) = $op(n.geom,val)
    @eval $op(n::CustomNode,val) = $op(node_val(n),val)
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
components(n::AssemblyNode)         = n.components
add_component!(n::AssemblyNode,p)   = push!(n.components,p)
has_component(n::AssemblyNode,id)   = haskey(components(n),id)
child_transform(n::AssemblyNode,id) = components(n)[id]

struct TransportUnitNode <: SceneNode
    id::TransportUnitID
    geom::GeomNode
    assembly::Pair{AssemblyID,CoordinateTransformations.Transformation}
    robots::TransformDict{BotID} # must be filled with unique invalid ids
end
TransportUnitNode(n::TransportUnitNode,geom::GeomNode) = TransportUnitNode(n.id,geom,n.assembly,n.robots)
TransportUnitNode(id,geom,assembly) = TransportUnitNode(id,geom,assembly,TransformDict{BotID}())
TransportUnitNode(id,geom,assembly_id::AssemblyID) = TransportUnitNode(id,geom,assembly_id=>identity_linear_map())
robot_team(n::TransportUnitNode) = n.robots
assembly_id(n::TransportUnitNode) = n.assembly.first

has_component(n::TransportUnitNode,id::AssemblyID)  = id == assembly_id(n)
has_component(n::TransportUnitNode,id::BotID)       = haskey(robot_team(n),id)
child_transform(n::TransportUnitNode,id::AssemblyID) = n.assembly.second
child_transform(n::TransportUnitNode,id::BotID)     = robot_team(n)[id]
add_robot!(n::TransportUnitNode,p)                  = push!(robot_team(n),p)
add_robot!(n::TransportUnitNode,r,t)                = add_robot!(n,r=>t)

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
        nodes = map(copy, tree.nodes),
        vtx_map = deepcopy(tree.vtx_map), # TODO would a shallow copy be fine?
        vtx_ids = deepcopy(tree.vtx_ids) # TODO would a shallow copy be fine?
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
        set_child!(tree,u,v)
    end
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
