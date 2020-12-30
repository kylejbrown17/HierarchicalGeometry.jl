let
    m = GridDiscretization(SVector(zeros(3)...),SVector(ones(3)...))

    GridOccupancy(m,trues(2,2,2))
    o1 = GridOccupancy(m,trues(2,2,2),SVector(0,0,0))
    @test (o1 + 1).offset == SVector(1,1,1)

    o2 = GridOccupancy(m,trues(2,2,2),SVector(1,1,1))
    @test has_overlap(o1,o2)
    @test has_overlap(o1+[0,1,1],o2)
    @test !has_overlap(o1,o2 + 1)
    @test !has_overlap(o1 - 1,o2)
    rect = overapproximate(o1,Hyperrectangle)
    @test m.origin in rect
    @test !isempty(intersection(Singleton(m.origin),rect))
end
let
    m = GridDiscretization(SVector(zeros(3)...),SVector(ones(3)...))
    @test HierarchicalGeometry.cell_indices(m,[0.0,0.0,0.0]) == [0,0,0]
    HierarchicalGeometry.get_hyperrectangle(m,[0,0,0])

    @test HierarchicalGeometry.cell_indices(m,[1.0,0.0,0.0]) == [1,0,0]
    @test HierarchicalGeometry.cell_indices(m,[-1.0,0.0,0.0]) == [-1,0,0]
end
let
    geom = LazySets.Ball2(zeros(3),4.0)
    hpoly = overapproximate(geom,equatorial_overapprox_model())
    @test is_intersection_empty(hpoly,LazySets.translate(hpoly,[9.0,0.0,0.0]))
    @test !is_intersection_empty(hpoly,LazySets.translate(hpoly,[2.0,0.0,0.0]))
    vpoly = tovrep(hpoly)
    overapproximate(hpoly,Ball2(LazySets.center(hpoly),1.0))
    overapproximate(hpoly,Ball2)

    m = GridDiscretization(SVector(zeros(3)...),SVector(ones(3)...))
    occupancy = overapproximate(hpoly,m)
    @test !is_intersection_empty(occupancy,occupancy+2)
    @test is_intersection_empty(occupancy,occupancy+8)
end
let
    a = Hyperrectangle([0.0,0.0],[1.0,1.0])
    b = Hyperrectangle([3.0,0.0],[1.0,1.0])
    c = Hyperrectangle([3.0,3.0],[1.0,1.0])
    d = Hyperrectangle([2.0,0.0],[1.5,0.25])
    @test isapprox(LazySets.distance(a,b), 1.0)
    @test isapprox(LazySets.distance(a,c), sqrt(2))
    @test isapprox(LazySets.distance(a,d), -0.5)
end
let
    a = Ball2([0.0,0.0],1.0)
    b = Ball2([2.0,0.0],1.0)
    c = Ball2([1.0,0.0],0.5)
    d = Ball2([0.0,0.0],0.5)
    @test isapprox(LazySets.distance(a,b), 0.0)
    @test isapprox(LazySets.distance(a,c), -0.5)
    @test isapprox(LazySets.distance(a,d), -1.5)

    @test isapprox(LazySets.distance(a,b), HierarchicalGeometry.distance_lower_bound(a,b))

    # @test isapprox(LazySets.distance(GeometryCollection([a]),GeometryCollection([b,c])), -0.5)
end
let
    a = GeomNode(Hyperrectangle([0.0,0.0],[1.0,1.0]))
    b = GeomNode(Hyperrectangle([3.0,0.0],[1.0,1.0]))
    @test isapprox(HierarchicalGeometry.distance_lower_bound(a,b), 1.0)
end
# let
#     g = GeometryHierarchy()
#     geom = LazySets.Ball2(zeros(3),1.0)
#     construct_geometry_tree!(g,geom)
#
#     for v in LightGraphs.vertices(g)
#         n = get_node(g,v)
#         LazySets.translate(n.geom,[1.0,0.0,0.0])
#     end
#
#     g = GeometryHierarchy()
#     geom2 = LazySets.translate(geom,[3.0,0.0,0.0])
#     add_node!(g,GeomNode(geom2),:BaseGeom)
#     construct_geometry_tree!(g2,geom2)
#
#     distance_lower_bound(geom,geom2)
#
#     # add_node!(g,GeomNode(geom),:BaseGeom)
#     # add_child_approximation!(g,equatorial_overapprox_model(),:BaseGeom,:Polyhedron)
# end
# let
#     table = HierarchicalGeometry.CollisionTable{Int}()
#     geoms = map(i->construct_geometry_tree(GeometryHierarchy(),Ball2(zeros(3),1.0))), 1:3)
#     transforms = [[0.0,0.0,0.0],[4.0,0.0,0.0],[0.0,4.5,0.0],]
#
# end
