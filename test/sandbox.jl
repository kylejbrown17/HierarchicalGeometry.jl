using HierarchicalGeometry
using LazySets
using Polyhedra
# using GeometryTypes
using GeometryBasics
using MeshCat
using Meshing
using Colors


POLYHEDRON_MATERIAL = MeshPhongMaterial(color=RGBA{Float32}(1, 0, 0, 0.5))

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

setobject!(vis[:ball], GeometryBasics.Sphere(Point(ball.center...),ball.radius))
setobject!(vis[:polytope], Polyhedra.Mesh(poly), POLYHEDRON_MATERIAL)

MeshCat.delete!(vis)
