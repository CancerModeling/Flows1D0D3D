cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "NoisyCH")
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
# Create executible
# ****************************************************************************
add_executable(${ProjectName} main.cpp)
target_link_libraries(${ProjectName} ${ProjectLib})
#target_link_libraries(${ProjectName} ${ProjectLib}
#        ${LIBMESH_LIBRARIES}
#        ${PETSC_LIBRARIES} ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# ****************************************************************************
# Create a default input file to execute
# ****************************************************************************
add_custom_command(
        TARGET ${ProjectName} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/noisy_ch_input.in $<TARGET_FILE_DIR:${ProjectName}>/noisy_ch_input.in
        COMMENT "Copying config file: ${CMAKE_CURRENT_SOURCE_DIR}/noisy_ch_input.in")