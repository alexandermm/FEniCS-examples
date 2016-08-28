#Example file to show how to generate mesh with pydistmesh and write it to a file type compatible with fenics
#Alex Martinez-Marchese

# Python imports

import numpy as np
import matplotlib.pyplot as plt

# Local imports
import dolfin as df
import distmesh as dm

#Cone measurements
rad = 0.01
vfr = 0.004

#Functions
def fdout(p):
    return dm.dunion(dm.dunion(dm.dcircle(p,0,0,rad), dm.drectangle(p,0.0,2*rad, -rad,rad)), dm.dcircle(p,2*rad,0,rad))

def fd(p):
    return dm.ddiff(dm.ddiff(fdout(p), dm.dcircle(p,0,0,vfr)), dm.dcircle(p,2*rad,0,vfr))

def fh(p):
    return dm.dunion(dm.dunion(0.0004-0.3*fdout(p), 0.0004+0.3*dm.dcircle(p,0,0,vfr)), 0.0004+0.3*dm.dcircle(p,2*rad,0,vfr))





#Make mesh
np.random.seed(1) # Always the same results
plt.ion()
p, t = dm.distmesh2d(fd, fh, 0.0004, (-0.01,-0.03, 0.03,0.03), 0.001, [(rad,-rad),(2*rad,-rad), (rad,rad),(2*rad,rad)])


# Write mesh as xml file
numVertices = p.shape[0]
numCells    = t.shape[0]

editor = df.MeshEditor()
mesh   = df.Mesh()

dim = 2

editor.open(mesh, 2, 2)            # top. and geom. dimension are both 3
editor.init_vertices(numVertices)  # number of vertices
editor.init_cells(numCells)        # number of cells

for x in range(0, numVertices):
	editor.add_vertex(x, p[x][:])

for x in range(0, numCells):
	editor.add_cell(x, np.array(t[x][:], dtype=np.uintp))

editor.close()

#Plot mesh using dolfin
df.plot(mesh)
df.interactive()

#Write to file
df.File('twovortexmesh.xml') << mesh



