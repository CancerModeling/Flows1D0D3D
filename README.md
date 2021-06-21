# Flows1D0D3D

Application which (currently) simulates bloodflows with coupled 1D, 0D and 3D models.
Its future goal is to simulate and study breast cancer models.

## Directory structure

- apps: contains the demo applications copied to bin/macrocirculation
- data/
  - Brats_data_processed_by_DF: ?
  - Breast_cancer_data_CW : ?
  - meshes: Contains the 1D meshes and boundary conditions
- tools: contains utility scripts
  - mesh_creation: scripts for creating the 1D mesh 
  - script
    - docker: dockerfiles
    - install: automatic install scripts to facilitate setup
- src: the source code of the library 
- tests: A few catch2 unittests for critical components

## Run unittests

The unittests can be run right after building the project with ctest
```
cmake <cmake-options-here> ..
make 
ctest
```

## Petsc solver types
We can select the following solvers with option `-ksp_type`
	
	- cg - Conjugate Gradient
	- gmres - Generalized Minimal Residual method
	- bcgs - Stabilized version of BiConjugate Gradient
	- tfqmr - A transpose free QMR (quasi minimal residual)
	- cgs - Conjugate Gradient Squared

more can be found in the petsc documentation.
