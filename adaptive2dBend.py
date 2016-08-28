#Example fenics file to show flow in a 2d bend with automatic mesh refinement
#Alex Martinez-Marchese

# Copyright (C) 2010 Marie E. Rognes
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Anders Logg 2011
#
# First added:  2010
# Last changed: 2011-10-04

from dolfin import *
from mshr   import *

#Dimenions of duct
width    = 1.0
lenghIn  = 4.0
lenghOut = 2.0
innerRad = 0.5

#Max flow speed at center of duct (for 2D laminar flow solution used in inlet)
Umax = 1.0

#Mesh resolution
res = 30

noSlipTol = 1e-5

class Noslip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and not ((x[0] < noSlipTol and x[1] > 0.0 and x[1] < width) or (x[1] > (width+innerRad+lenghOut) - noSlipTol and x[0] > lenghIn+innerRad and x[0] < lenghIn+innerRad+width))

class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < noSlipTol

class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > (width+innerRad+lenghOut) - noSlipTol


# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

# Allow approximating values for points that may be generated outside
# of domain (because of numerical inaccuracies)
parameters["allow_extrapolation"] = True

# Choose refinement algorithm (needed for solver to choose right refinement algorithm)
parameters["refinement_algorithm"] = "plaza_with_parent_facets"



# Material parameters
nu = Constant(0.002)

# Mesh
circleIn     = Circle(Point(lenghIn, innerRad+width), innerRad)
circleOut    = Circle(Point(lenghIn, innerRad+width), innerRad+width)
andRectagle  = Rectangle(Point(lenghIn,0.0), Point(lenghIn+innerRad+width,innerRad+width))
rectagleIn   = Rectangle(Point(0.0,0.0), Point(lenghIn,width))
rectagleOut  = Rectangle(Point(lenghIn+innerRad,width+innerRad), Point(lenghIn+innerRad+width,width+innerRad+lenghOut))
domain       = (circleOut-circleIn) * andRectagle + rectagleIn + rectagleOut

mesh     = generate_mesh(domain, res)
plot(mesh, title="2D mesh")
interactive()

# Create boundary subdomains
bc_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
bc_markers.set_all(3)
Inflow().mark(bc_markers, 1)
Outflow().mark(bc_markers, 2)
plot(bc_markers)
interactive()

# Define new measure with associated subdomains
ds = Measure("ds")[bc_markers]

# Define function spaces (Taylor-Hood)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define unknown and test function(s)
(v, q) = TestFunctions(W)
w = Function(W)
(u, p) = (as_vector((w[0], w[1])), w[2])

# Define variational forms for Navier-Stokes including boundary velocities perpendicular to boundaries (L term)
n = FacetNormal(mesh)
a = (nu*inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx()
a = a + inner(grad(u)*u, v)*dx()
L = - p*dot(v, n)*ds()
F = a - L

# Define boundary conditions
inletflow = Expression(("Umax*(1.0 - (x[1]-radius)*(x[1]-radius)/radius/radius)","0.0"), Umax=Umax, radius=width/2)
bcnp  = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Noslip())
bcin  = DirichletBC(W.sub(0), inletflow, Inflow())
bcout = DirichletBC(W.sub(1), Constant(0.0), Outflow())
bc = [bcnp, bcin, bcout]

# Define goal
M = p*ds(1)

# Define error tolerance (with respect to goal)
tol = 1e-04

# Compute Jacobian form
J = derivative(F, w)

# Define variational problem
pde = NonlinearVariationalProblem(F, w, bc, J)

# Define solver
solver = AdaptiveNonlinearVariationalSolver(pde, M)

# Solve to given tolerance
solver.solve(tol)

# Show all timings
list_timings()

# Solutions on coarsest and finest mesh:
mesh0 = mesh.root_node()
mesh1 = mesh.leaf_node()
plot(mesh0, title="Initial mesh")
plot(mesh1, title="Final mesh")

(u0, p0) = w.root_node().split()
(u1, p1) = w.leaf_node().split()
plot(u0, title="Velocity on initial mesh")
plot(u1, title="Velocity on final mesh")
plot(p0, title="Pressure on initial mesh")
plot(p1, title="Pressure on final mesh")
interactive()


