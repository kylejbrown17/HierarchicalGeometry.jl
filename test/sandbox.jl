using HierarchicalGeometry
using LazySets
using Polyhedra
using MeshCat
# using GeometryTypes
using GeometryBasics

vis = Visualizer()
render(vis)

ball = LazySets.Ball2(zeros(3),1.0)
ball2 = LazySets.Ball2([1.0,0.0,0.0],1.0)
lazy_set = UnionSet(ball,ball2)

bbox = overapproximate(ball)
hpoly = convert(LazySets.HPolytope,bbox)
overapproximate(hpoly,Hyperrectangle)
overapproximate(hpoly,LazySets.BallInf)
vpoly = tovrep(hpoly)
poly = Polyhedra.polyhedron(hpoly)

translated_hpoly = LazySets.Translation(hpoly,[1.0,1.0,1.0])


# model = PolyhedronOverapprox(3,8)
model = HierarchicalGeometry.equatorial_overapprox_model()
# model = HierarchicalGeometry.equatorial_overapprox_model([0.0],0.0:π/2:2π)
hpoly = overapproximate(lazy_set,model)
poly = Polyhedra.polyhedron(hpoly)

setobject!(vis, GeometryBasics.Sphere(Point(ball.center...),ball.radius))
setobject!(vis, Polyhedra.Mesh(poly))

MeshCat.delete!(vis)
