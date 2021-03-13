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

# test bifurcation flow:
add_executable(${ProjectName}BifurcationFlow test_bifurcation_flow.cpp)
target_link_libraries(${ProjectName}BifurcationFlow ${ProjectLib})

# test the convergence:
add_executable(${ProjectName}ConvergenceStudy convergence_study.cpp)
target_link_libraries(${ProjectName}ConvergenceStudy ${ProjectLib})

foreach(TargetName ${ProjectName}NonlinearFlowLine ${ProjectName}NonlinearFlowBifurcation ${ProjectName}ConvergenceStudy)
    add_custom_command(
            TARGET ${TargetName} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory $<TARGET_FILE_DIR:${TargetName}>/output)

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/animation.py $<TARGET_FILE_DIR:${TargetName}>/animation.py
            COMMENT "Copying file: ${CMAKE_CURRENT_SOURCE_DIR}/animation.py")
endforeach()
