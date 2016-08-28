# FEniCS-examples
Example code for solving the Navier-Stokes (NS) equations using [FEniCS] (https://fenicsproject.org/), with the adaptive and plain nonlinear variational solver. The solver is able to solve for any steady laminar flow using finite elements. The adaptive version can subdivide elements based on an error estimate of an averaged parameter. For example, one can set the flow profile in one boundary such as an inlet and subdivide the mesh until the averaged pressure on that boundary is 10% of the actual estimated value. The pressure can be set to zero at the outlet. More details are available [here] (https://fenicsproject.org/featured/2011/automated_error_control.html). The paper describing the adaptive algorithm is [here] (http://arxiv.org/abs/1205.3096).

## Non-adaptive solutions:
* 2dDuct.py: pipe flow on a straight 2d duct. Solves the full nonlinear incompressible NS equations using the non adaptive nonlinear solver. Note that to get the right solution, the variational form of the NS equations has an extra term, in order to enforce fluid moving normal to inlets and outlets. More on this [here] (https://harishnarayanan.org/research/navier-stokes/). The solution is the expected analytic parabolic velocity profile along the duct, with a linear pressure drop.
* 3dDuct.py: same as above, but in 3D, with the mesh generated using the meshing library [mshr] (https://bitbucket.org/fenics-project/mshr/wiki/Home), included with FEniCS. This problem also has a simple analytic solution for checking the program.
* 2dBend.py: pipe flow in an elbow. The code shows how to combine the geometric shape primitives of the mshr library.
* 3dBend.py: same as 2dBend.py but in 3D and using a sharp bend due to two intersecting cylinders.

##Adaptive solutions:
* adaptive2dBend.py: same code as 2dBend.py but including a few more lines for using the adaptive variational solver. 
* adaptive3dBend.py: same code as 3dBend.py. The code needs a while to run, and needs sufficient memory due to the final number of elements.
* adaptiveDrivenCavity.py: shows how to run the known driven cavity problem using the adaptive variational solver. Same code as the one found in the FeniCS library with minor modifications.
* adaptive2Vortex.py: shows one way to include 2D vortices in FEniCS, how to use the adaptive variational solver with vortices and how to fix the pressure at a single node (corresponding to a stagnation point in the flow).
* 2VortexMesh.py: shows how to make an initialy refined mesh using [PyDistMesh] (https://github.com/bfroehle/pydistmesh) and writing it to a FEniCS compatible mesh file format.
