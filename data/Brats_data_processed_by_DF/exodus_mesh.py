# dependencies
import os
import sys
import numpy as np
from scipy.io import loadmat
import meshio as ms
import argparse
import vtk

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Convert mesh from .mat to .e format')
  parser.add_argument('--in_file',
                      default='TestCase_breastwholemesh_dless_coarse_mesh.mat',
                      type=str,
                      help='input .mat file')
  parser.add_argument('--voxel_size',
                      default=0.0267,
                      type=float,
                      help='size of voxel (default is 0.0267 cm)')
  args = parser.parse_args()

  mat_file = args.in_file
  out_file = mat_file.replace('.mat', '') + '_voxel_size_%Lf' % (args.voxel_size)
  
  # read data
  data = loadmat(mat_file)

  # get nodes and elems
  nodes = args.voxel_size * data['node'] # scale the nodal coordinates
  elements = data['elem'][:, 0:4] - 1 # subtract 1 to start indexing from 0
  print('node and element data shapes: {}, {}'.format(nodes.shape, elements.shape))

  # create meshio mesh object
  elements = elements.astype('uint32')
  mesh = ms.Mesh(nodes, [('tetra', elements)])
  print(mesh)

  # save as exodus format
  # first create .xdmf file
  mesh.write(out_file + '.xdmf')

  vtkFemReader = vtk.vtkXdmfReader() 
  vtkFemReader.SetFileName(out_file + '.xdmf')
  vtkFemReader.Update() 

  vtkFEMImageWriter = vtk.vtkExodusIIWriter() 
  # vtkFEMImageWriter.SetFileTypeToBinary() 
  vtkFEMImageWriter.SetInputData( vtkFemReader.GetOutput() )
  vtkFEMImageWriter.SetFileName( out_file + '.e')
  vtkFEMImageWriter.Update() 

  # clear .xdmf file
  os.remove(out_file + '.xdmf' )
  os.remove(out_file + '.h5' )



