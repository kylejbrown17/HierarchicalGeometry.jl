# Test Transformations
let
    ball = LazySets.Ball2([1.0,0.0,0.0],1.0)
    bbox = overapproximate(ball,Hyperrectangle)
    hpoly = convert(LazySets.HPolytope,bbox)
    vpoly = convert(LazySets.VPolytope,bbox)
    poly = Polyhedra.polyhedron(hpoly)

    # setobject!(vis[:ball], GeometryBasics.Sphere(ball))
    # setobject!(vis[:bbox], GeometryBasics.HyperRectangle(bbox))
    # setobject!(vis[:polytope], Polyhedra.Mesh(poly), POLYHEDRON_MATERIAL)

    # Translation
    t = CoordinateTransformations.Translation(1.0,2.0,3.0)
    @test array_isapprox(ball.center.+t.translation.data,transform(ball,t).center)
    @test isapprox(ball.radius,transform(ball,t).radius)
    @test array_isapprox(bbox.center.+t.translation.data,transform(bbox,t).center)
    @test isapprox(bbox.radius,transform(bbox,t).radius)
    for (v1,v2) in zip(vertices_list(vpoly),vertices_list(transform(vpoly,t)))
        @test array_isapprox(v1+t.translation,v2)
    end
    for (v1,v2) in zip(vertices_list(vpoly),vertices_list(tovrep(transform(hpoly,t))))
        @test array_isapprox(v1+t.translation,v2)
    end

    # settransform!(vis[:ball], a)
    # settransform!(vis[:polytope], a)

    # delete!(vis[:transformedball])
    # setobject!(vis[:transformed_ball],GeometryBasics.Sphere(transform(ball,a)))
    # setobject!(vis[:transformed_bbox],GeometryBasics.HyperRectangle(transform(bbox,a)))
    # setobject!(vis[:transformed_polytope],Polyhedra.Mesh(Polyhedra.polyhedron(transform(vpoly,a))), POLYHEDRON_MATERIAL)
    # setobject!(vis[:transformed_polytope],Polyhedra.Mesh(Polyhedra.polyhedron(transform(hpoly,a))), POLYHEDRON_MATERIAL)


    # Rotation
    r = CoordinateTransformations.LinearMap(RotZ(π/4))
    @test array_isapprox(transform(ball.center,r),transform(ball,r).center)
    @test isapprox(ball.radius,transform(ball,r).radius)
    @test array_isapprox(transform(bbox.center,r),transform(bbox,r).center)
    @test isapprox(bbox.radius[1]*sqrt(2),transform(bbox,r).radius[1])
    for (v1,v2) in zip(vertices_list(vpoly),vertices_list(transform(vpoly,r)))
        @test array_isapprox(Vector(r(v1)),v2)
    end
    for (v1,v2) in zip(vertices_list(vpoly),vertices_list(tovrep(transform(hpoly,r))))
        @test array_isapprox(Vector(r(v1)),v2)
    end

    # settransform!(vis[:ball], t)
    # settransform!(vis[:polytope], t)

    # Combined
    a = compose(t,r) # Rotate first, then translate
    @test array_isapprox(transform(ball.center,a),transform(ball,a).center)
    @test isapprox(ball.radius,transform(ball,a).radius)
    @test array_isapprox(transform(bbox.center,a),transform(bbox,a).center)
    @test isapprox(bbox.radius[1]*sqrt(2),transform(bbox,a).radius[1])
    for (v1,v2) in zip(vertices_list(vpoly),vertices_list(transform(vpoly,a)))
        @test array_isapprox(Vector(a(v1)),v2)
    end
    for (v1,v2) in zip(vertices_list(vpoly),vertices_list(tovrep(transform(hpoly,a))))
        @test array_isapprox(Vector(a(v1)),v2)
    end
end
# Test Transform Tree
let
    # tree = TransformTree{Symbol}()
    tree = GraphUtils.CustomTree{TransformNode,Symbol}()
    root = TransformNode()
    add_node!(tree,root,:ONE)
    GraphUtils.add_child!(tree,deepcopy(root),:TWO,:ONE)
    GraphUtils.add_child!(tree,deepcopy(root),:THREE,:TWO)
    t = compose(CoordinateTransformations.Translation(1,0,0),CoordinateTransformations.LinearMap(RotZ(0)))
    for v in LightGraphs.vertices(tree)
        set_local_transform!(get_node(tree,v),t) # Not connected to tree, so can't globally update
    end
    for v in LightGraphs.vertices(tree)
        data = [global_transform(tree,v).translation.data...]
        @test array_isapprox(data,[0.0, 0.0, 0.0])
    end
    update_transform_tree!(tree,:ONE)
    for v in LightGraphs.vertices(tree)
        data = [global_transform(tree,v).translation.data...]
        @test array_isapprox(data,[v,0.0,0.0])
    end
end
let
    # tree = TransformTree{Symbol}()
    tree = GraphUtils.CustomTree{TransformNode,Symbol}()
    root = TransformNode()
    add_node!(tree,root,:ONE)
    GraphUtils.add_child!(tree,deepcopy(root),:TWO,:ONE)
    GraphUtils.add_child!(tree,deepcopy(root),:THREE,:ONE)
    t = compose(CoordinateTransformations.Translation(1,0,0),CoordinateTransformations.LinearMap(RotZ(0)))
    for v in LightGraphs.vertices(tree)
        set_local_transform!(tree,v,t)
    end
    @test array_isapprox([global_transform(tree,:ONE).translation.data...],[1.0,0,0])
    @test array_isapprox([global_transform(tree,:TWO).translation.data...],[2.0,0,0])
    @test array_isapprox([global_transform(tree,:THREE).translation.data...],[2.0,0,0])
    change_parent!(tree,:THREE,:TWO)
    @test has_edge(tree,:TWO,:THREE)
    @test !has_edge(tree,:ONE,:THREE)
    for v in LightGraphs.vertices(tree)
        set_local_transform!(tree,v,t)
    end
    @test array_isapprox([global_transform(tree,:ONE).translation.data...],[1.0,0,0])
    @test array_isapprox([global_transform(tree,:TWO).translation.data...],[2.0,0,0])
    @test array_isapprox([global_transform(tree,:THREE).translation.data...],[3.0,0,0])
end
# Test CachedElement
let
    c = HierarchicalGeometry.CachedElement(1)
    @test get_element(c) == 1
    set_up_to_date!(c,false)
    @test !is_up_to_date(c)
    update_element!(c,2)
    @test is_up_to_date(c)
    @test get_element(c) == 2
end
# Test GeomNode
let
    geom = Ball2(ones(3),1.0)
    n = GeomNode(geom)
    @test array_isapprox(get_base_geom(n).center,geom.center)
    set_up_to_date!(n,true)
    @test is_up_to_date(n)
    t = compose(CoordinateTransformations.Translation(1.0,0,0),CoordinateTransformations.LinearMap(RotZ(0)))
    set_global_transform!(n,t)
    @test !is_up_to_date(n)
    @test array_isapprox(get_cached_geom(n).center, geom.center .+ [t.translation.data...])
    @test is_up_to_date(n)
    set_global_transform!(n,t)
    @test !is_up_to_date(n)
end
# relative transformations
let
    a = compose(CoordinateTransformations.Translation(1,0,0),CoordinateTransformations.LinearMap(RotZ(0.0)))
    b = compose(CoordinateTransformations.Translation(1,4,0),CoordinateTransformations.LinearMap(RotZ(0.75*π)))
    rot_err = Rotations.rotation_error(a,b)
    t = HierarchicalGeometry.relative_transform(a,b)
    @test array_isapprox(inv(t.linear),b.linear;rtol=1e-12,atol=1e-12)

    a = compose(CoordinateTransformations.Translation(1,0,0),CoordinateTransformations.LinearMap(RotXYZ(0.0,0.0,0.0)))
    b = compose(CoordinateTransformations.Translation(1,4,0),CoordinateTransformations.LinearMap(RotXYZ(0.1,0.0,1.0*π)))
    t = HierarchicalGeometry.relative_transform(a,b)
end
