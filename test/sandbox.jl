using HierarchicalGeometry
using LazySets
using Polyhedra
using MeshCat
# using GeometryTypes
using GeometryBasics

vis = Visualizer()
render(vis)

ball = LazySets.Ball2(zeros(3),2.0)

bbox = overapproximate(ball)
hpoly = convert(LazySets.HPolytope,bbox)
overapproximate(hpoly,Hyperrectangle)
overapproximate(hpoly,LazySets.BallInf)
vpoly = tovrep(hpoly)
poly = Polyhedra.polyhedron(hpoly)

translated_hpoly = LazySets.Translation(hpoly,[1.0,1.0,1.0])


ball2 = LazySets.Ball2([3.0,0.0,0.0],1.0)
lazy_set = UnionSet(ball,ball2)
model = HierarchicalGeometry.equatorial_overapprox_model()
hpoly = overapproximate(lazy_set,model)
poly = Polyhedra.polyhedron(hpoly)

setobject!(vis, GeometryBasics.Sphere(Point(ball.center...),ball.radius))
setobject!(vis, Polyhedra.Mesh(poly))

MeshCat.delete!(vis)
