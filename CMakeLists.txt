cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(TumorModels CXX)

# ****************************************************************************
# Append directory containing important modules
# ****************************************************************************
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# ****************************************************************************
# Package search
# ****************************************************************************
# MPI
find_package(MPI REQUIRED)      # Built-in in CMake
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)

# pthreads
find_package(Threads REQUIRED)

# libmesh
find_package(LibMesh REQUIRED)

# petsc
# Note: -DPETSC_DIR should point to the same directory that libmesh was compiled with
find_package(Petsc REQUIRED)

# debug output flags
message("LIBMESH_DEFINITIONS = ${PC_LIBMESH_CFLAGS_DEFS}")
message("LIBMESH_FLAGS = ${PC_LIBMESH_CFLAGS_FLAGS}")
message("LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")
message("LIBMESH_LIBRARIES = ${LIBMESH_LIBRARIES}")
message("PETSC_LIB = ${PETSC_LIB}")


# VTK
find_package(VTK REQUIRED)
#find_package(VTK COMPONENTS vtkCommonCore QUIET)

if (VTK_VERSION VERSION_LESS "6")
    message(FATAL_ERROR "ERROR: Only VTK library 6+ supported")
endif ()
message("VTK Lib: ${VTK_LIBRARIES}")
message("VTK Dir: ${VTK_DIR}")
message("VTK Include: ${VTK_INCLUDE_DIRS}")

# issue error and exit
if ((NOT MPI_FOUND) AND (NOT LIBMESH_FOUND))
    message(FATAL_ERROR "ERROR: MPI and LIBMESH are required")
endif ()

#if ("${PETSC_DIR}" STREQUAL "")
#    message(FATAL_ERROR "ERROR: Please specify petsc directory")
#endif ()

# ****************************************************************************
# Build, compiler, linker flags
# ****************************************************************************
# enable_language(CXX)
# set(CMAKE_BUILD_TYPE RelWithDebInfo)

#message("Build type: ${CMAKE_BUILD_TYPE}")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_FLAGS} -DNDEBUG")
    elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_FLAGS} -DNDEBUG")
    endif()
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_FLAGS} -xHost -std=c++14")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_FLAGS} -DNDEBUG -std=c++11 -stdlib=libc++")
endif ()
#
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")


# ****************************************************************************
# Include directories
# ****************************************************************************
# include_directories(${PETSC_DIR} ${EIGEN_DIR})
# Note that LIBMESH_INCLUDE_DIRS contains the path to petsc and eigen directory
include_directories(${LIBMESH_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

# ****************************************************************************
# Add subprojects
# ****************************************************************************
set(LIBTG_BUILD_FLAG "Build_All" CACHE STRING "Specify which if all or subset
of models to build")

message("LIBTG_BUILD_FLAG: ${LIBTG_BUILD_FLAG}")

message("Gnerating build files for LibTG ")
add_subdirectory(models)

set(TEST_BUILD TRUE CACHE BOOL "Build Test")
if(${TEST_BUILD})
    message("Gnerating build files for Test ")
    add_subdirectory(models/tests)
endif()

message("Adding applications ")
add_subdirectory(applications)


INSTALL (
        DIRECTORY ${CMAKE_SOURCE_DIR}/models
        DESTINATION include
        FILES_MATCHING PATTERN "*.h*")

INSTALL (
        DIRECTORY ${CMAKE_SOURCE_DIR}/external
        DESTINATION include
        FILES_MATCHING PATTERN "*.h*")