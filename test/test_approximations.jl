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

    # ball = LazySets.Ball2(zeros(3),4.0)
    # ball2 = LazySets.Ball2([1.0,0.0,0.0],1.0)
    # lazy_set = UnionSet(ball,ball2)
    geom = LazySets.Ball2(zeros(3),4.0)
    hpoly = overapproximate(geom,equatorial_overapprox_model())
    occupancy = overapproximate(hpoly,m)

    @test !is_intersection_empty(occupancy,occupancy+2)
    @test is_intersection_empty(occupancy,occupancy+8)

end
