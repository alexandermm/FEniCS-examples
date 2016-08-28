# FEniCS-examples
Solutions to the Navier-Stokes equations using [FEniCS] (https://fenicsproject.org/) using the adaptive and plain nonlinear variational solver. The solver is able to solve for any steady laminar flow using finite elements. The adaptive version can subdivide elements based on an error estimate of an averaged parameter. For example, one can set the flow profile in one boundary such as an inlet and subdivide the mesh until the averaged pressure on that boundary is 10% of the actual estimated value. The pressure can be set to zero at the outlet.

File descriptions:
* 2dDuct.py: pipe flow on a straight 2d duct. Solves the nonlinear Navier-Stokes equations using the non adaptive nonlinear solver. The solution is the expected analytic parabolic velocity profile along the duct, with a linear pressure drop.
* 3dDuct.py: same as above, but in 3D, with the mesh generated using the meshing library msch, included with FEniCS. This problem also has a simple analytic solution for checking the program.
