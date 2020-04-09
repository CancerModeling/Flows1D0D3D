cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(LibTG CXX)

# Create executible in main directory of build directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories("${PROJECT_SOURCE_DIR}/utils")
include_directories(${CMAKE_SOURCE_DIR}/external/gmm)

# ****************************************************************************
# Add subprojects
# ****************************************************************************
add_subdirectory(utils)
add_subdirectory(AvaFV)
add_subdirectory(NetFV)
add_subdirectory(NetFVFE)
add_subdirectory(NetFC)

# # ****************************************************************************
# # Create executible
# # ****************************************************************************
# add_executable(LibTG main.cpp ../docs/temp.txt)
# target_link_libraries(LibTG Utils LibAvaLima LibAngLima LibAvaMech LibNetTum
#         ${LIBMESH_LIBRARIES}
#         ${PETSC_LIB} ${MPI_CXX_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
