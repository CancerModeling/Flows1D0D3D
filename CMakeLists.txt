#////////////////////////////////////////////////////////////////////////////////
#//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
#//
#//  Distributed under the Boost Software License, Version 1.0. (See accompanying
#//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#////////////////////////////////////////////////////////////////////////////////

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(TumorModels LANGUAGES CXX)

set(CMAKE_PROJECT_DESCRIPTION "Collection of tumor growth models")

# version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_UPDATE 0)

# ****************************************************************************
# Project setting
# ****************************************************************************
# cmake scripts to help search dependencies
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# Create executible in main directory of build directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# set default build type "Release"
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# set flag for information dump during cmake build
set(Enable_CMAKE_Debug_Build TRUE CACHE BOOL "Output important information during cmake for isolating build issues")

# set flag for documentation
set(Enable_Documentation FALSE CACHE BOOL "Generate documentation")

# set flag that enables or disables tests
set(Enable_Tests FALSE CACHE BOOL "Build test executibles and perform tests")

# add flag which allows user to restrict building only specific models
set(LIBTG_BUILD_FLAG "Build_NetFVFE" CACHE STRING "Specify which if all or subset
of models to be built")

# ****************************************************************************
# Testing
# ****************************************************************************
enable_testing()

if (${Enable_CMAKE_Debug_Build})
    message("Enable_Tests = ${Enable_Tests}")
endif()

# ****************************************************************************
# Package search
# ****************************************************************************

# **** bash ****#
find_program (BASH_PROGRAM bash REQUIRED)

# **** MPI ****#
find_package(MPI REQUIRED)
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)

# **** pthreads ****#
find_package(Threads REQUIRED)

# **** libmesh ****#
find_package(LibMesh REQUIRED)

# **** petsc ****#
# Note: -DPETSC_DIR should point to the same directory that libmesh was compiled with
find_package(Petsc REQUIRED)

# json
if (USE_EXTERNAL_JSON)
    find_package(nlohmann_json 3.2.0 REQUIRED)
else ()
    include(FetchContent)
    set(JSON_BuildTests OFF)
    FetchContent_Declare(json
            GIT_REPOSITORY https://github.com/nlohmann/json.git
            GIT_TAG v3.7.3)
    FetchContent_GetProperties(json)
    if(NOT json_POPULATED)
        FetchContent_Populate(json)
        add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
endif()

# cxxopts
include(FetchContent)
FetchContent_Declare(cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts
        GIT_TAG a20bda6179b41c0cf49d557b0b21c3fb4f6c09b7)
FetchContent_GetProperties(cxxopts)
if(NOT cxxopts_POPULATED)
    FetchContent_Populate(cxxopts)
    add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# fmt
FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 7.1.3)
FetchContent_GetProperties(fmt)
if(NOT fmt_POPULATED)
    FetchContent_Populate(fmt)
    add_subdirectory(${fmt_SOURCE_DIR} ${fmt_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# debug output flags
if (${Enable_CMAKE_Debug_Build})
message("LIBMESH_DEFINITIONS = ${PC_LIBMESH_CFLAGS_DEFS}")
message("LIBMESH_FLAGS = ${PC_LIBMESH_CFLAGS_FLAGS}")
message("LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")
message("LIBMESH_LIBRARIES = ${LIBMESH_LIBRARIES}")
message("PETSC_LIBRARIES = ${PETSC_LIBRARIES}")
endif()

# **** vtk ****#
find_package(VTK REQUIRED)

if (VTK_VERSION VERSION_LESS "6")
    message(FATAL_ERROR "ERROR: Only VTK library 6+ supported")
endif ()

if (${Enable_CMAKE_Debug_Build})
    message("VTK_DIR = ${VTK_DIR}")
    message("VTK_LIBRARIES = ${VTK_LIBRARIES}")
    message("VTK_INCLUDE_DIRS = ${VTK_INCLUDE_DIRS}")
endif()

# issue error and exit
if ((NOT MPI_FOUND) AND (NOT LIBMESH_FOUND))
    message(FATAL_ERROR "ERROR: MPI and LIBMESH are required")
endif ()

#if ("${PETSC_DIR}" STREQUAL "")
#    message(FATAL_ERROR "ERROR: Please specify petsc directory")
#endif ()

# nanoflann
FetchContent_Declare(nanoflann
        GIT_REPOSITORY https://github.com/jlblancoc/nanoflann
        GIT_TAG v1.3.2)
FetchContent_GetProperties(nanoflann)
if(NOT nanoflann_POPULATED)
    FetchContent_Populate(nanoflann)
    add_subdirectory(${nanoflann_SOURCE_DIR} ${nanoflann_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
# ****************************************************************************
# Build, compiler, linker flags
# ****************************************************************************
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG ")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -fPIC -L /usr/lib/libmesh-1.6.1/lib/libtimpi_opt.so.2")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_FLAGS} -xHost -std=c++14")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -std=c++11 -stdlib=libc++")
endif ()

if (${Enable_CMAKE_Debug_Build})
    message("Build type: ${CMAKE_BUILD_TYPE}")
    message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
    message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_CXX_STANDARD = ${CMAKE_CXX_STANDARD}")
endif()


# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories(${LIBMESH_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

# ****************************************************************************
# Add subprojects
# ****************************************************************************
if (${Enable_CMAKE_Debug_Build})
    message("LIBTG_BUILD_FLAG: ${LIBTG_BUILD_FLAG}")
endif()

if (${Enable_CMAKE_Debug_Build})
    message("Gnerating build files for LibTG ")
endif()
add_subdirectory(models)

if(${Enable_Tests})
    message("Gnerating build files for Test ")
    add_subdirectory(models/tests)
    add_subdirectory(tests)
endif()

if (${Enable_CMAKE_Debug_Build})
    message("Adding applications ")
endif()
add_subdirectory(applications)

# documentation
if(${Enable_Documentation})
    if (${Enable_CMAKE_Debug_Build})
        message("Adding Documentation ")
    endif()
    add_subdirectory(docs)
endif()

# ****************************************************************************
# Post build operations
# ****************************************************************************
# test directory
if(${Enable_Tests})
    if (${Enable_CMAKE_Debug_Build})
        message("Copying test data to build directory ")
    endif()

    add_custom_target("copydata" ALL)
    add_custom_command(
            TARGET "copydata" POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/tests/test_data/
            ${PROJECT_BINARY_DIR}/tests/test_data/)
endif()
