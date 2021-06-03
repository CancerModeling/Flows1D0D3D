# About data
- `Brats18_2013_13_1_t1ce.nii.gz` - T1CE brats data
- `vesseldilate.nii.gz` - Vascular domain segmented by David
- `mri_mesh.m` - Matlab script for creating mesh from T1CE data using iso2mesh
- `mesh4.mat` - Mesh obtained by above matlab script
- `exodus_mesh.py` - Script to convert mesh from .mat to .e format. For dependencies and other information, see directory `Breast_cancer_data_CW`.

We use mesh from T1 CE image and model Darcy's flow. We consider heterogeneous parameter with one constant in vascular domain and other in extra-vascular domain. 
