cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "test1")

# Create executible in main directory of build directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/apps/sensitivity/)


# ****************************************************************************
# Create executible
# ****************************************************************************
add_executable(${ProjectName} main.cpp)
target_link_libraries(${ProjectName} LibNetFVFE
        ${LIBMESH_LIBRARIES}
        ${PETSC_LIB} ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
