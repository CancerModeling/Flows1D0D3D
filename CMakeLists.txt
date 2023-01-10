#////////////////////////////////////////////////////////////////////////////////
#//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
#//
#//  Distributed under the Boost Software License, Version 1.0. (See accompanying
#//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#////////////////////////////////////////////////////////////////////////////////

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(TumorModels LANGUAGES CXX)

set(CMAKE_PROJECT_DESCRIPTION "Collection of tumor growth models")

# version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_UPDATE 0)

# pybind11 needs position independent code:
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
endif ()

# set flag for information dump during cmake build
set(Enable_CMAKE_Debug_Build TRUE CACHE BOOL "Output important information during cmake for isolating build issues")

# set flag for documentation
set(Enable_Documentation TRUE CACHE BOOL "Generate documentation")

# set flag that enables or disables tests
set(Enable_Tests TRUE CACHE BOOL "Build test executibles and perform tests")

# ****************************************************************************
# Package search
# ****************************************************************************

# **** bash ****#
find_program(BASH_PROGRAM bash REQUIRED)

# **** MPI ****#
find_package(MPI REQUIRED)
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
# for backward compatibility it make sense to add
# see https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
if(NOT TARGET MPI::MPI_CXX)
    add_library(MPI::MPI_CXX IMPORTED INTERFACE)

    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_FLAGS})
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
endif()

# **** pthreads ****#
find_package(Threads REQUIRED)

# **** petsc ****#
# Note: -DPETSC_DIR should point to the same directory that libmesh was compiled with
find_package(Petsc REQUIRED)

# **** eigen ***#
find_package (Eigen3 3.3 REQUIRED)

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
    if (NOT json_POPULATED)
        FetchContent_Populate(json)
        add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif ()
endif ()

# cxxopts
include(FetchContent)
FetchContent_Declare(cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts
        GIT_TAG a20bda6179b41c0cf49d557b0b21c3fb4f6c09b7)
FetchContent_GetProperties(cxxopts)
if (NOT cxxopts_POPULATED)
    FetchContent_Populate(cxxopts)
    add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()

# ****************************************************************************
# Build, compiler, linker flags
# ****************************************************************************
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG ")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -fPIC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -std=c++14")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -std=c++11 -stdlib=libc++")
endif ()

if (${Enable_CMAKE_Debug_Build})
    message("Build type: ${CMAKE_BUILD_TYPE}")
    message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
    message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_CXX_STANDARD = ${CMAKE_CXX_STANDARD}")
endif ()


# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories(${MPI_INCLUDE_PATH})

add_library(lmgmm INTERFACE)
target_include_directories(lmgmm INTERFACE external/gmm)

# ****************************************************************************
# Add subfolders
# ****************************************************************************
if (${Enable_CMAKE_Debug_Build})
    message("LIBTG_BUILD_FLAG: ${LIBTG_BUILD_FLAG}")
endif ()

add_subdirectory(src)

set(Enable_Apps ON CACHE BOOL "Includes the apps in our build process")
if (${Enable_Apps})
    add_subdirectory(apps)
endif()


# ****************************************************************************
# Testing
# ****************************************************************************
if (${Enable_CMAKE_Debug_Build})
    message("Enable_Tests = ${Enable_Tests}")
endif ()
if (${Enable_Tests})
    enable_testing()

    add_subdirectory(tests)

    find_package(Catch2 QUIET)
    # if cmake cannot find the package, we download it
    if (NOT Catch2_FOUND)
        # get catch2 for testing
        include(FetchContent)
        FetchContent_Declare(
                catch2_content
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_TAG v2.13.9
        )
        FetchContent_MakeAvailable(catch2_content)
    endif()
endif ()

# documentation
if (${Enable_Documentation})
    if (${Enable_CMAKE_Debug_Build})
        message("Adding Documentation ")
    endif ()
    add_subdirectory(docs)
endif ()

# python bindings
set(Enable_Python_Bindings OFF CACHE BOOL "Adds the python bindings to our build process")
if (${Enable_Python_Bindings})
    add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11)
    add_subdirectory(${CMAKE_SOURCE_DIR}/python)
endif ()
