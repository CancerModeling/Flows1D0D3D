cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(Test CXX)

# Create executible in main directory of build directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories("${CMAKE_SOURCE_DIR}/models/utils")
include_directories("${CMAKE_SOURCE_DIR}/models/tests")

# ****************************************************************************
# Add subprojects
# ****************************************************************************
add_subdirectory(TestComm)
add_subdirectory(TestMesh)
add_subdirectory(TestGeom)

AUX_SOURCE_DIRECTORY(./ SOURCES)

# ****************************************************************************
# Create executible
# ****************************************************************************
add_executable(Test ${SOURCES})
target_link_libraries(Test TestComm TestMesh TestGeom Utils)
#        ${LIBMESH_LIBRARIES} ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
