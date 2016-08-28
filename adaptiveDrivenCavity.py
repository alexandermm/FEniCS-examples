#Example fenics file to show flow in a driven cavity with automatic mesh refinement
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



# Mesh
n = 20
mesh = UnitSquareMesh(n,n,"crossed")

#Lid speed
U = 0.2



class NoSlip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS)

class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (x[1] > 1.0 - DOLFIN_EPS)


# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

# Allow approximating values for points that may be generated outside
# of domain (because of numerical inaccuracies)
parameters["allow_extrapolation"] = True

# Choose refinement algorithm (needed for solver to choose right refinement algorithm)
parameters["refinement_algorithm"] = "plaza_with_parent_facets"


# Material parameters
nu = Constant(0.05)


# Create boundary subdomains
bc_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
bc_markers.set_all(0)
NoSlip().mark(bc_markers, 1)
Lid().mark(bc_markers, 2)

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
F = (nu*inner(grad(u), grad(v)) - div(v)*p + q*div(u) + inner(grad(u)*u, v))*dx() + p*dot(v, n)*ds()


# Define boundary conditions
lidu   = Expression(("Upeak","0.0"), Upeak=U)
bcu    = DirichletBC(W.sub(0), lidu, Lid())
bcp    = DirichletBC(W.sub(1), Constant(0.0), "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")
noslip = DirichletBC(W.sub(0), Constant((0.0, 0.0)), NoSlip())
bc = [bcu, bcp, noslip]

# Define goal
M = p*ds(2)

# Define error tolerance (with respect to goal)
tol = 1.e-05

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
plot(p0, title="Pressure on initial mesh")
plot(p1, title="Pressure on final mesh")
plot(u0, title="Velocity on initial mesh")
plot(u1, title="Velocity on final mesh")
interactive()


