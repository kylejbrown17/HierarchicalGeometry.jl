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
