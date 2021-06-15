# dependencies
import os
import sys
import numpy as np
from scipy.io import loadmat
import meshio as ms
import argparse
import vtk

if __name__ == "__main__":
  
  vtkFemReader = vtk.vtkGenericDataObjectReader() 
  vtkFemReader.SetFileName('vesseldilate.vtk')
  vtkFemReader.Update() 

  vtkFEMImageWriter = vtk.vtkExodusIIWriter() 
  # vtkFEMImageWriter.SetFileTypeToBinary() 
  vtkFEMImageWriter.SetInputData( vtkFemReader.GetOutput() )
  vtkFEMImageWriter.SetFileName('vesseldilate.e')
  vtkFEMImageWriter.Update() 

  # clear .xdmf file
  # os.remove(out_file + '.xdmf' )
  # os.remove(out_file + '.h5' )



