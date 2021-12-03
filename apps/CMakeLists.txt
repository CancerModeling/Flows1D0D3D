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
target_link_libraries(${ProjectName}NonlinearFlowLine ${ProjectLib} PkgConfig::PETSC) 

# nonlinear flow solver bifurcation:
add_executable(${ProjectName}NonlinearFlowBifurcation nonlinear_flow_bifurcation.cpp)
target_link_libraries(${ProjectName}NonlinearFlowBifurcation ${ProjectLib} PkgConfig::PETSC)

# test the convergence:
add_executable(${ProjectName}ConvergenceStudy convergence_study.cpp)
target_link_libraries(${ProjectName}ConvergenceStudy ${ProjectLib} PkgConfig::PETSC)

# test the 1d network:
add_executable(${ProjectName}BreastGeometry breast_geometry.cpp)
target_link_libraries(${ProjectName}BreastGeometry ${ProjectLib} PkgConfig::PETSC)

# test the 1d network:
add_executable(${ProjectName}Nonlinear1DSolver nonlinear_1d_solver.cpp)
target_link_libraries(${ProjectName}Nonlinear1DSolver ${ProjectLib} PkgConfig::PETSC)
target_link_libraries(${ProjectName}Nonlinear1DSolver cxxopts)

add_executable(${ProjectName}CombinedGeometry combined_geometry.cpp)
target_link_libraries(${ProjectName}CombinedGeometry ${ProjectLib} PkgConfig::PETSC)
target_link_libraries(${ProjectName}CombinedGeometry cxxopts)

add_executable(${ProjectName}Calibration calibration.cpp)
target_link_libraries(${ProjectName}Calibration ${ProjectLib} PkgConfig::PETSC)
target_link_libraries(${ProjectName}Calibration cxxopts)

add_executable(${ProjectName}LinearFlow linear_flow.cpp)
target_link_libraries(${ProjectName}LinearFlow ${ProjectLib} PkgConfig::PETSC)

add_executable(${ProjectName}LinearFlowBreastGeometry linear_flow_breast_geometry.cpp)
target_link_libraries(${ProjectName}LinearFlowBreastGeometry ${ProjectLib})
target_link_libraries(${ProjectName}LinearFlowBreastGeometry PkgConfig::PETSC)
target_link_libraries(${ProjectName}LinearFlowBreastGeometry cxxopts)

add_executable(${ProjectName}ImplicitTransport implicit_transport.cpp)
target_link_libraries(${ProjectName}ImplicitTransport ${ProjectLib})
target_link_libraries(${ProjectName}ImplicitTransport PkgConfig::PETSC)
target_link_libraries(${ProjectName}ImplicitTransport cxxopts)

add_executable(${ProjectName}ImplicitTransportWith0D implicit_transport_with_0d.cpp)
target_link_libraries(${ProjectName}ImplicitTransportWith0D ${ProjectLib})
target_link_libraries(${ProjectName}ImplicitTransportWith0D PkgConfig::PETSC)
target_link_libraries(${ProjectName}ImplicitTransportWith0D cxxopts)

add_executable(${ProjectName}WindkesselTransportInOutTest windkessel_transport_in_out_test.cpp)
target_link_libraries(${ProjectName}WindkesselTransportInOutTest ${ProjectLib})
target_link_libraries(${ProjectName}WindkesselTransportInOutTest PkgConfig::PETSC)
target_link_libraries(${ProjectName}WindkesselTransportInOutTest cxxopts)

add_executable(${ProjectName}WindkesselTransportThroughCompartments windkessel_transport_through_compartments.cpp)
target_link_libraries(${ProjectName}WindkesselTransportThroughCompartments ${ProjectLib})
target_link_libraries(${ProjectName}WindkesselTransportThroughCompartments PkgConfig::PETSC)
target_link_libraries(${ProjectName}WindkesselTransportThroughCompartments cxxopts)

add_executable(${ProjectName}ExplicitImplicitScheme explicit_implicit_scheme.cpp)
target_link_libraries(${ProjectName}ExplicitImplicitScheme ${ProjectLib})
target_link_libraries(${ProjectName}ExplicitImplicitScheme PkgConfig::PETSC)
target_link_libraries(${ProjectName}ExplicitImplicitScheme cxxopts)

add_executable(${ProjectName}ExplicitImplicitWithTransport explicit_implicit_with_transport.cpp)
target_link_libraries(${ProjectName}ExplicitImplicitWithTransport ${ProjectLib})
target_link_libraries(${ProjectName}ExplicitImplicitWithTransport PkgConfig::PETSC)
target_link_libraries(${ProjectName}ExplicitImplicitWithTransport cxxopts)

add_executable(${ProjectName}Transport transport_test.cpp)
target_link_libraries(${ProjectName}Transport ${ProjectLib})

foreach (TargetName ${ProjectName}NonlinearFlowLine
        ${ProjectName}NonlinearFlowBifurcation
        ${ProjectName}ConvergenceStudy
        ${ProjectName}BreastGeometry
        ${ProjectName}Nonlinear1DSolver
        ${ProjectName}CombinedGeometry
        ${ProjectName}Transport
        ${ProjectName}LinearFlow
        )

    foreach (FolderName output
            output_darcy3d
            output_darcy_hetero
            output_part_perfusion_test
            output_voxels_to_mesh
            output_darcy_flow_reg_outlets
            output_full_1d0d3d_pkj
            output_3d_model
            output_multiscale_time_3d1d0d)
        add_custom_command(
                TARGET ${TargetName} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory
                $<TARGET_FILE_DIR:${TargetName}>/${FolderName})
    endforeach ()

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization $<TARGET_FILE_DIR:${TargetName}>/visualization
            COMMENT "Symlinking folder: ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization")

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../data $<TARGET_FILE_DIR:${TargetName}>/data
            COMMENT "Copying symlinking folder: data")

    file(GLOB visualization_scripts RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/*)
    foreach (f ${visualization_scripts})
        MESSAGE(${f})
        add_custom_command(
                TARGET ${TargetName} PRE_BUILD
                COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/${f} $<TARGET_FILE_DIR:${TargetName}>/${f}
                COMMENT "Symlinking : ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/${f} to $<TARGET_FILE_DIR:${TargetName}>/${f}")
    endforeach ()
endforeach ()
