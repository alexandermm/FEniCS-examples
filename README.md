# FEniCS-examples
Example code for solving the Navier-Stokes (NS) equations using [FEniCS] (https://fenicsproject.org/), with the adaptive and plain nonlinear variational solver. The solver is able to solve for any steady laminar flow using finite elements. The adaptive version can subdivide elements based on an error estimate of an averaged parameter. For example, one can set the flow profile in one boundary such as an inlet and subdivide the mesh until the averaged pressure on that boundary is 10% of the actual estimated value. The pressure can be set to zero at the outlet. More details are available [here] (https://fenicsproject.org/featured/2011/automated_error_control.html). The paper describing the adaptive algorithm is [here] (http://arxiv.org/abs/1205.3096).

## Non-adaptive solutions:
* 2dDuct.py: pipe flow on a straight 2d duct. Solves the full nonlinear incompressible NS equations using the non adaptive nonlinear solver. Note that to get the right solution, the variational form of the NS equations has an extra term, in order to enforce fluid moving normal to inlets and outlets. More on this [here] (https://harishnarayanan.org/research/navier-stokes/). The solution is the expected analytic parabolic velocity profile along the duct, with a linear pressure drop.
* 3dDuct.py: same as above, but in 3D, with the mesh generated using the meshing library msch, included with FEniCS. This problem also has a simple analytic solution for checking the program.

##Adaptive solutions:
* Hello.
* Hello.
