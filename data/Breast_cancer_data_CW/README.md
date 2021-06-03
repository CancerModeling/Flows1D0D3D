# Chengyue's comments

Instruction for the data:

1. "TestCase_VesselNetwork.mat": matlab file storing vascular network. Variable structure as following
  Network -- Links -- n1: index of the node at one end of the vessel
                   -- n2: index of the node at the other end of the vessel
                   -- dx: step width of discretization along the vessel path (Unit = [cm])
                   -- YXZ: 3D coordinations of discretized points along the vessel (Unit = [cm])
                   -- radii: local vascular radii (Unit = [cm])
                   -- Ktrans: spatial-resolved material property, volume transform constant, corresponding to the vascular permeability (Unit = ml blood / min / ml tissue)
                   -- junc: label of ends of vessels being terminal (NaN) or branching points (number)
                   -- pv0: initial condition of blood pressure, indicating blood flow direction along the network
                   
          -- Nodes -- conn_L: indices of vessels connected to this node
                   -- conn_N: indices of nodes connected to this nodes
                   -- YXZ: 3D coordinations of this node (Unit = [cm])
                   -- term: label indicating this node being a terminal (1) or branching point (0)
                  
2. Convert xml to vtu
```py
import meshio as ms

# load xml file
fname = 'TestCase_breastwholemesh_dless'
mesh = ms.read(fname + '.xml')

# voxel to meter
mesh.points = 0.0267 * 0.01 * mesh.points

# save
ms.write(fname + '.vtu', mesh)
```

> Modified Nov 04, 2020
     Adding matlab files storing originally generated mesh information
"TestCase_breastwholemesh_dless.mat": matlab file storing the entire breast mesh. Included variables as following
      node2: Node coordinates of the tetrahedral mesh. Each row is a node. Columns 1 - 3 refer to x-, y-, z-coordinates, respectively
      elem2: Element list of the tetrahedral mesh. Each row is a element. Columns 1 - 4 refer to the indices of vertices (nodes). Last column is the ID of isolated regions
      face2: Mesh surface element list of the tetrahedral mesh. Each row is a facet. Columns 1 - 3 refer to the indices of vertices. Last column denotes the boundary ID (all 1 in this dataset)
      
"TestCase_tumormesh_dless.mat": matlab file storing the tumor mesh. Variables are the same as "TestCase_breastwholemesh_dless.mat". 
   

# Script to create exodus mesh file from the mat mesh file
`exodus_mesh.py` takes in the mat input file and outputs mesh in exodus format for libmesh.
```sh
python3 exodus_mesh.py --in_file TestCase_breastwholemesh_dless_coarse_mesh.mat 
```

Script depends on `python3, numpy, meshio, scipy, vtk`. You can install these dependencies in a conda environment using 
```sh
conda env create -f pymesh.yml
conda activate pymesh
python3 exodus_mesh.py --in_file TestCase_breastwholemesh_dless_coarse_mesh.mat 
```
