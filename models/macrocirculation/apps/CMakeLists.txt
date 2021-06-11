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

# advection solver:
add_executable(${ProjectName}AdvectionImplicit advection_implicit.cpp)
target_link_libraries(${ProjectName}AdvectionImplicit ${ProjectLib})

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
add_executable(${ProjectName}Abstract33VesselsGeometry abstract_33_vessels_geometry.cpp)
target_link_libraries(${ProjectName}Abstract33VesselsGeometry ${ProjectLib})
target_link_libraries(${ProjectName}Abstract33VesselsGeometry cxxopts)

add_executable(${ProjectName}CombinedGeometry combined_geometry.cpp)
target_link_libraries(${ProjectName}CombinedGeometry ${ProjectLib})

add_executable(${ProjectName}Calibration calibration.cpp)
target_link_libraries(${ProjectName}Calibration ${ProjectLib})
target_link_libraries(${ProjectName}Calibration cxxopts)

add_executable(${ProjectName}LinearFlow linear_flow.cpp)
target_link_libraries(${ProjectName}LinearFlow ${ProjectLib})
target_link_libraries(${ProjectName}LinearFlow Eigen3::Eigen)

add_executable(${ProjectName}Transport transport_test.cpp)
target_link_libraries(${ProjectName}Transport ${ProjectLib})

# Darcy's flow in 3D tissue domain
add_executable(${ProjectName}Darcy3D tissue_domain_darcy.cpp)
target_link_libraries(${ProjectName}Darcy3D ${ProjectLib})
target_link_libraries(${ProjectName}Darcy3D cxxopts)

foreach(TargetName ${ProjectName}NonlinearFlowLine
                   ${ProjectName}NonlinearFlowBifurcation
                   ${ProjectName}ConvergenceStudy
                   ${ProjectName}AdvectionImplicit
                   ${ProjectName}BreastGeometry
                   ${ProjectName}Abstract33VesselsGeometry
                   ${ProjectName}CombinedGeometry
                   ${ProjectName}Transport
                   ${ProjectName}LinearFlow
                   ${ProjectName}Darcy3D)
    add_custom_command(
            TARGET ${TargetName} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory $<TARGET_FILE_DIR:${TargetName}>/output)

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/animation.py $<TARGET_FILE_DIR:${TargetName}>/animation.py
            COMMENT "Copying file: ${CMAKE_CURRENT_SOURCE_DIR}/animation.py")

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/plot_vessels.py $<TARGET_FILE_DIR:${TargetName}>/plot_vessels.py
            COMMENT "Copying file: ${CMAKE_CURRENT_SOURCE_DIR}/plot_vessels.py")

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/plot_windkessel.py $<TARGET_FILE_DIR:${TargetName}>/plot_windkessel.py
            COMMENT "Copying file: ${CMAKE_CURRENT_SOURCE_DIR}/plot_windkessel.py")

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../data $<TARGET_FILE_DIR:${ProjectName}CombinedGeometry>/data
            COMMENT "Copying symlinking folder: data")
endforeach()
