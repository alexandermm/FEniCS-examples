#Example fenics file to show flow in a 2d duct
#Alex Martinez-Marchese 

#Copyright (C) 2010 Marie E. Rognes
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

#Dimenions of duct
width = 1.0
lengh = 4.0

#Max flow speed at center of duct (for 2D laminar flow solution used in inlet)
Umax = 1.0

#Mesh resolution
res = 10

class Noslip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[1] < DOLFIN_EPS or x[1] > width - DOLFIN_EPS)

class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > lengh - DOLFIN_EPS

# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

# Allow approximating values for points that may be generated outside
# of domain (because of numerical inaccuracies)
parameters["allow_extrapolation"] = True

# Material parameters
nu = Constant(0.02)

# Mesh
mesh = RectangleMesh(Point(0.0, 0.0), Point(lengh, width), 10, 10)

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

# Compute Jacobian form
J = derivative(F, w)

# Define variational problem
pde = NonlinearVariationalProblem(F, w, bc, J)

# Define solver
solver = NonlinearVariationalSolver(pde)

# Solve to given tolerance
solver.solve()

# Show all timings
list_timings()

# Extract solutions on coarsest and finest mesh:
(u, p) = w.split()
plot(mesh, title="Mesh")
plot(u, title="Velocity")
plot(p, title="Pressure")
interactive()
