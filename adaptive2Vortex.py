#Example fenics file to show flow due to two counter-rotating vortices with automatic mesh refinement
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


#Cone measurements
rad = 0.01
vfr = 0.004

#Peak core speed
Upeak = 1.0


meshtol = 1e-5

class InR1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt(x[0]*x[0] + x[1]*x[1]) < (vfr + meshtol))

class InR2(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt((x[0]-2.0*rad)*(x[0]-2.0*rad) + x[1]*x[1]) < (vfr + meshtol))

class NoSlip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (not (sqrt(x[0]*x[0] + x[1]*x[1]) < (vfr + meshtol))) and (not (sqrt((x[0]-2.0*rad)*(x[0]-2.0*rad) + x[1]*x[1]) < (vfr + meshtol)))

class PPoint(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], rad) and x[1] < -rad + meshtol


# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

# Allow approximating values for points that may be generated outside
# of domain (because of numerical inaccuracies)
parameters["allow_extrapolation"] = True

# Choose refinement algorithm (needed for solver to choose right refinement algorithm)
parameters["refinement_algorithm"] = "plaza_with_parent_facets"


# Material parameters
nu = Constant(0.00002)

# Mesh
mesh = Mesh('2vortexmesh.xml')

# Create boundary subdomains
bc_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
bc_markers.set_all(0)
InR1().mark(bc_markers, 1)
InR2().mark(bc_markers, 1)
NoSlip().mark(bc_markers, 2)
#plot(bc_markers)

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
F = (nu*inner(grad(u), grad(v)) - div(v)*p + q*div(u) + inner(grad(u)*u, v))*dx()

# Define boundary conditions
coreu1 = Expression(("Upeak*x[1]/vfr","-Upeak*x[0]/vfr"), Upeak=Upeak, vfr=vfr)
bcu1   = DirichletBC(W.sub(0), coreu1, InR1())

coreu2 = Expression(("-Upeak*x[1]/vfr","Upeak*(x[0]-2.0*rad)/vfr"), Upeak=Upeak, vfr=vfr, rad=rad)
bcu2   = DirichletBC(W.sub(0), coreu2, InR2())

bcp    = DirichletBC(W.sub(1), Constant(0.0), PPoint(), "pointwise")
noslip = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoSlip())

bc = [bcu1, bcu2, bcp, noslip]

# Define goal
M = p*ds(1)

# Define error tolerance (with respect to goal)
tol = 1e-06

# Compute Jacobian form
J = derivative(F, w)

# Define variational problem
pde = NonlinearVariationalProblem(F, w, bc, J)

# Define solver
solver = AdaptiveNonlinearVariationalSolver(pde, M)
#solver = NonlinearVariationalSolver(pde)

# Solve to given tolerance
solver.solve(tol)
#solver.solve()

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
#(u, p) = w.split()
#plot(u, title="Velocity field")
#plot(p, title="Pressure")
interactive()


