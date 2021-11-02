cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "Macrocirculation")
set(ProjectLib "Lib${ProjectName}")
project(${ProjectName} CXX)

# Create executible in main directory of build directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/macrocirculation/)

# ****************************************************************************
# Create executable
# ****************************************************************************

# nonlinear flow solver line:
add_executable(${ProjectName}NonlinearFlowLine nonlinear_flow_line.cpp)
target_link_libraries(${ProjectName}NonlinearFlowLine ${ProjectLib})

# nonlinear flow solver bifurcation:
add_executable(${ProjectName}NonlinearFlowBifurcation nonlinear_flow_bifurcation.cpp)
target_link_libraries(${ProjectName}NonlinearFlowBifurcation ${ProjectLib})

# test the convergence:
add_executable(${ProjectName}ConvergenceStudy convergence_study.cpp)
target_link_libraries(${ProjectName}ConvergenceStudy ${ProjectLib})

# test the 1d network:
add_executable(${ProjectName}BreastGeometry breast_geometry.cpp)
target_link_libraries(${ProjectName}BreastGeometry ${ProjectLib})

# test the 1d network:
add_executable(${ProjectName}Nonlinear1DSolver nonlinear_1d_solver.cpp)
target_link_libraries(${ProjectName}Nonlinear1DSolver ${ProjectLib})
target_link_libraries(${ProjectName}Nonlinear1DSolver cxxopts)

add_executable(${ProjectName}CombinedGeometry combined_geometry.cpp)
target_link_libraries(${ProjectName}CombinedGeometry ${ProjectLib})
target_link_libraries(${ProjectName}CombinedGeometry cxxopts)

add_executable(${ProjectName}Full1D0D3DModel full_1D0D3D_model.cpp)
target_link_libraries(${ProjectName}Full1D0D3DModel ${ProjectLib})
target_link_libraries(${ProjectName}Full1D0D3DModel cxxopts fmt)

add_executable(${ProjectName}Linearized1D0D3DModel linearized_1D0D3D_model.cpp)
target_link_libraries(${ProjectName}Linearized1D0D3DModel ${ProjectLib})
target_link_libraries(${ProjectName}Linearized1D0D3DModel cxxopts fmt)

add_executable(${ProjectName}Full1D0D3DModelPKJ full_1D0D3D_model_pkj.cpp)
target_link_libraries(${ProjectName}Full1D0D3DModelPKJ ${ProjectLib})
target_link_libraries(${ProjectName}Full1D0D3DModelPKJ cxxopts fmt)

add_executable(${ProjectName}Calibration calibration.cpp)
target_link_libraries(${ProjectName}Calibration ${ProjectLib})
target_link_libraries(${ProjectName}Calibration cxxopts)

add_executable(${ProjectName}LinearFlow linear_flow.cpp)
target_link_libraries(${ProjectName}LinearFlow ${ProjectLib})
target_link_libraries(${ProjectName}LinearFlow ${LIBMESH_LIBRARIES}) # libmesh for eigen

add_executable(${ProjectName}LinearFlowBreastGeometry linear_flow_breast_geometry.cpp)
target_link_libraries(${ProjectName}LinearFlowBreastGeometry ${ProjectLib})
target_link_libraries(${ProjectName}LinearFlowBreastGeometry ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}LinearFlowBreastGeometry cxxopts)

add_executable(${ProjectName}ImplicitTransport implicit_transport.cpp)
target_link_libraries(${ProjectName}ImplicitTransport ${ProjectLib})
target_link_libraries(${ProjectName}ImplicitTransport ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}ImplicitTransport cxxopts)

add_executable(${ProjectName}ImplicitTransportWith0D implicit_transport_with_0d.cpp)
target_link_libraries(${ProjectName}ImplicitTransportWith0D ${ProjectLib})
target_link_libraries(${ProjectName}ImplicitTransportWith0D ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}ImplicitTransportWith0D cxxopts)

add_executable(${ProjectName}WindkesselTransportInOutTest windkessel_transport_in_out_test.cpp)
target_link_libraries(${ProjectName}WindkesselTransportInOutTest ${ProjectLib})
target_link_libraries(${ProjectName}WindkesselTransportInOutTest ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}WindkesselTransportInOutTest cxxopts)

add_executable(${ProjectName}WindkesselTransportThroughCompartments windkessel_transport_through_compartments.cpp)
target_link_libraries(${ProjectName}WindkesselTransportThroughCompartments ${ProjectLib})
target_link_libraries(${ProjectName}WindkesselTransportThroughCompartments ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}WindkesselTransportThroughCompartments cxxopts)

add_executable(${ProjectName}ExplicitImplicitScheme explicit_implicit_scheme.cpp)
target_link_libraries(${ProjectName}ExplicitImplicitScheme ${ProjectLib})
target_link_libraries(${ProjectName}ExplicitImplicitScheme ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}ExplicitImplicitScheme cxxopts)

add_executable(${ProjectName}ExplicitImplicitWithTransport explicit_implicit_with_transport.cpp)
target_link_libraries(${ProjectName}ExplicitImplicitWithTransport ${ProjectLib})
target_link_libraries(${ProjectName}ExplicitImplicitWithTransport ${LIBMESH_LIBRARIES}) # libmesh for eigen
target_link_libraries(${ProjectName}ExplicitImplicitWithTransport cxxopts)

add_executable(${ProjectName}Transport transport_test.cpp)
target_link_libraries(${ProjectName}Transport ${ProjectLib})

# Darcy's flow in 3D tissue domain
add_executable(${ProjectName}Darcy3D tissue_domain_darcy.cpp)
target_link_libraries(${ProjectName}Darcy3D ${ProjectLib})
target_link_libraries(${ProjectName}Darcy3D cxxopts fmt)

# Test nifti reader
add_executable(${ProjectName}NiftiRW test_nifti_reader.cpp)
target_link_libraries(${ProjectName}NiftiRW ${ProjectLib})
target_link_libraries(${ProjectName}NiftiRW cxxopts)

# Model Darcy's flow with heterogeneous parameter field
add_executable(${ProjectName}DarcyHetero darcy_test_heterogeneous.cpp)
target_link_libraries(${ProjectName}DarcyHetero ${ProjectLib})
target_link_libraries(${ProjectName}DarcyHetero cxxopts)

add_executable(${ProjectName}LegacyVTK legacy_vtk_reader.cpp)
target_link_libraries(${ProjectName}LegacyVTK ${ProjectLib})
target_link_libraries(${ProjectName}LegacyVTK cxxopts)

# test perfusion partitioning
add_executable(${ProjectName}PartPerfusionTest partioning_perfusion.cpp)
target_link_libraries(${ProjectName}PartPerfusionTest ${ProjectLib})
target_link_libraries(${ProjectName}PartPerfusionTest cxxopts nanoflann fmt)

# project voxel data on mesh
add_executable(${ProjectName}VoxelsToMesh voxels_to_mesh.cpp)
target_link_libraries(${ProjectName}VoxelsToMesh ${ProjectLib})
target_link_libraries(${ProjectName}VoxelsToMesh cxxopts nanoflann fmt)

# darcy's flow with regularized outlets
add_executable(${ProjectName}DarcyRegOutlets darcy_flow_regularized_outlets.cpp)
target_link_libraries(${ProjectName}DarcyRegOutlets ${ProjectLib})
target_link_libraries(${ProjectName}DarcyRegOutlets cxxopts nanoflann fmt)

# multiscale in time full simulation
add_executable(${ProjectName}MultiscaleTime3D1D0D multiscale_time_3d1d0d.cpp)
target_link_libraries(${ProjectName}MultiscaleTime3D1D0D ${ProjectLib})
target_link_libraries(${ProjectName}MultiscaleTime3D1D0D cxxopts fmt)

foreach(TargetName ${ProjectName}NonlinearFlowLine
                   ${ProjectName}NonlinearFlowBifurcation
                   ${ProjectName}ConvergenceStudy
                   ${ProjectName}BreastGeometry
                   ${ProjectName}Nonlinear1DSolver
                   ${ProjectName}CombinedGeometry
                   ${ProjectName}Transport
                   ${ProjectName}Darcy3D
                   ${ProjectName}NiftiRW
                   ${ProjectName}DarcyHetero
                   ${ProjectName}LegacyVTK
                   ${ProjectName}PartPerfusionTest
                   ${ProjectName}LinearFlow
                   ${ProjectName}VoxelsToMesh
                   ${ProjectName}DarcyRegOutlets
                   ${ProjectName}Full1D0D3DModel
                   ${ProjectName}Full1D0D3DModelPKJ
                   ${ProjectName}MultiscaleTime3D1D0D)

    foreach(FolderName output
            output_darcy3d
            output_darcy_hetero
            output_part_perfusion_test
            output_voxels_to_mesh
            output_darcy_flow_reg_outlets
            output_full_1d0d3d_pkj
            output_multiscale_time_3d1d0d)
        add_custom_command(
                TARGET ${TargetName} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory
                $<TARGET_FILE_DIR:${TargetName}>/${FolderName})
    endforeach()

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization $<TARGET_FILE_DIR:${TargetName}>/visualization
            COMMENT "Symlinking folder: ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization")

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../data $<TARGET_FILE_DIR:${TargetName}>/data
            COMMENT "Copying symlinking folder: data")

    file(GLOB visualization_scripts RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/*)
    foreach(f ${visualization_scripts})
        MESSAGE(${f})
        add_custom_command(
                TARGET ${TargetName} PRE_BUILD
                COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/${f} $<TARGET_FILE_DIR:${TargetName}>/${f}
                COMMENT "Symlinking : ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/${f} to $<TARGET_FILE_DIR:${TargetName}>/${f}")
    endforeach()
endforeach()
