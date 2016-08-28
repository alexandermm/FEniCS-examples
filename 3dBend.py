#Example fenics file to show flow in a 3d bend
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
radius1 = 0.5
lengh1  = 1.0
radius2 = 0.6
lengh2  = 3.0
gap     = 0.7

#Max flow speed at center of duct (for 3D laminar flow solution used in inlet)
Umax = 1.0

#Tolerance to mark no slip boundary condition
nosliptol = 1.0e-2

class Noslip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt(x[1]*x[1] + x[2]*x[2]) > radius1 - nosliptol) and x[2] < radius1+nosliptol \
                           or (on_boundary and x[0] > lengh1 and x[2] < radius1+nosliptol) \
                           or (x[2] >= radius1+nosliptol and (sqrt((x[0]-lengh1)*(x[0]-lengh1) + x[1]*x[1]) > radius2 - nosliptol))

class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < nosliptol

class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] > lengh2-gap - nosliptol

# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

# Allow approximating values for points that may be generated outside
# of domain (because of numerical inaccuracies)
parameters["allow_extrapolation"] = True

# Material parameters
nu = Constant(0.2)

# Mesh
cylinder1 = Cylinder(Point(0.0, 0.0, 0.0), Point(lengh1, 0.0, 0.0), radius1, radius1)
cylinder2 = Cylinder(Point(lengh1, 0.0, -gap), Point(lengh1, 0.0, lengh2-gap), radius2, radius2)
domain    = cylinder1+cylinder2
res      = 20
mesh     = generate_mesh(domain, res)
plot(mesh, title="3D mesh")
interactive()

# Create boundary subdomains
bc_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
bc_markers.set_all(3)
Inflow().mark(bc_markers, 1)
Outflow().mark(bc_markers, 2)
plot(bc_markers)
#Can save marked mesh and view it in slices in Paraview
#File("boundarymarkers.pvd") << bc_markers
interactive()

# Define function spaces (Taylor-Hood)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define unknown and test function(s)
(v, q) = TestFunctions(W)
w = Function(W)
(u, p) = (as_vector((w[0], w[1], w[2])), w[3])

# Define variational forms for Navier-Stokes including boundary velocities perpendicular to boundaries (L term)
# Define variational forms
n = FacetNormal(mesh)
a = (nu*inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx()
a = a + inner(grad(u)*u, v)*dx()
L = - p*dot(v, n)*ds()
F = a - L

# Define boundary conditions
inletflow = Expression(("Umax*(1.0 - (x[1]*x[1] + x[2]*x[2])/radius/radius)","0.0","0.0"), Umax=Umax, radius=radius1)
bcnp  = DirichletBC(W.sub(0), Constant((0.0, 0.0, 0.0)), Noslip())
bcin  = DirichletBC(W.sub(0), inletflow, Inflow())
bcout = DirichletBC(W.sub(1), Constant(0.0), Outflow())
bc = [bcnp, bcin, bcout]

# Compute Jacobian form
J = derivative(F, w)

# Define variational problem
pde = NonlinearVariationalProblem(F, w, bc, J)

# Define solver
solver = NonlinearVariationalSolver(pde)
#Stop Newton solver early to see if boundary condition is not just being copied to other parts of domain
#solver.parameters['newton_solver']['relative_tolerance'] = 1E-1

# Solve to given tolerance
solver.solve()

# Show all timings
list_timings()

# Extract solutions on coarsest and finest mesh:
(u, p) = w.split()
plot(mesh, title="Mesh")
plot(u, title="Velocity")
File("u.pvd") << u
plot(p, title="Pressure")
interactive()
