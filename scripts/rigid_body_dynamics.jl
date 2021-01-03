using RigidBodyDynamics
using LinearAlgebra
using StaticArrays

g = -9.81
world = RigidBody{Float64}("world")
robot1 = RigidBody{Float64}("robot1")

