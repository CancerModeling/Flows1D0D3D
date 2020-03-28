# Naming convention

- Net : Network-Tumor coupled model

- Ava : Avascular tumor model

- FV : Finite volume discretization

- FE : Finite element discretization

- FVFE : Finite volume finite element mixed discretization

- FC : Pressure and nutrient in 3d and 1d are solved simultaneously using block matrix

## Models implemented

### NetFVFE

Network-Tumor coupled modeled solved using mixed finite volume (FV) and finite element (FE)

Network-tumor coupled model. In this we solve nutrient and pressure in both 1d and 3d using finite volume method. Rest of the species are solved using bi-linear elements.

### NetFV

Network-tumor coupled model. All equations are discretized using finite volume method.

### NetFC

Network-tumor coupled model. Mixed FV-FE discretization. 1D-3D coupled field are solved simultaneously using block matrix.

### AvaFV

Avascular model which is trimmed down version of NetTum-FV. This allows one to test/debug the 3D tumor equations faster without having to solve 1D equations. 

### tests

Collection of tests for LibMesh features, geometrical functions, parallel communication

## Running the model

To run model `<model name>`:
  
  1. Either have `model_name = <model name>` in input.in, or do not provide the field `model_name` in input.in. 
  
  2. Run the model by calling `build/bin/<model name>` executible as follows
  ```sh
  <directory to build>/bin/<model name> input.in <libmesh options>
  ```
  where `<libmesh options>` are solver type, preconditioner type, etc. For example, `<libmesh options> = " -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 10"`

Currently, we have choices `AvaFV`, `NetFV`, `NetFVFE`, `NetFC`, `Test` for `<model name>`.
