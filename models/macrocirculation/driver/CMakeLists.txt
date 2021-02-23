cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "Macrocirculation")
set(ProjectLib "Lib${ProjectName}")
project(${ProjectName} CXX)

# Create executible in main directory of build directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories("${CMAKE_SOURCE_DIR}/models/utils")
include_directories("${CMAKE_SOURCE_DIR}/models/${ProjectName}")

# ****************************************************************************
# Create executable
# ****************************************************************************

# advection:
add_executable(${ProjectName}Advection main_advection.cpp)
target_link_libraries(${ProjectName}Advection ${ProjectLib})

# nonlinear flow solver line:
add_executable(${ProjectName}NonlinearFlowLine main_nonlinear_flow_line.cpp)
target_link_libraries(${ProjectName}NonlinearFlowLine ${ProjectLib})

# nonlinear flow solver bifurcation:
add_executable(${ProjectName}NonlinearFlowBifurcation main_nonlinear_flow_bifurcation.cpp)
target_link_libraries(${ProjectName}NonlinearFlowBifurcation ${ProjectLib})

# test bifurcation flow:
add_executable(${ProjectName}BifurcationFlow test_bifurcation_flow.cpp)
target_link_libraries(${ProjectName}BifurcationFlow ${ProjectLib})

# test the convergence:
add_executable(${ProjectName}Convergence test_convergence.cpp)
target_link_libraries(${ProjectName}Convergence ${ProjectLib})

# test the communicator:
add_executable(${ProjectName}Communicator test_communicator.cpp)
target_link_libraries(${ProjectName}Communicator ${ProjectLib})
