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
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/gmm")
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/aixlog")
include_directories(${VTK_INCLUDE_DIRS})

# our csv parser, included as a system dependencies to silence the warnings:
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/fast-cpp-csv-parser")

# ****************************************************************************
# Add subprojects
# ****************************************************************************
add_subdirectory(utils)

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}"
        STREQUAL "Build_Prashant" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_AvaFV")
    add_subdirectory(AvaFV)
endif ()

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}"
        STREQUAL "Build_Prashant" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_NetFV")
    add_subdirectory(NetFV)
endif ()

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}"
        STREQUAL "Build_Prashant" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_NetFVFE")
    add_subdirectory(NetFVFE)
endif ()

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}"
        STREQUAL "Build_Prashant" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_NetFCFVFE")
    add_subdirectory(NetFCFVFE)
endif ()

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}"
        STREQUAL "Build_Tobias" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_NetFC")
    add_subdirectory(NetFC)
endif ()

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_TwoSpecies")
    add_subdirectory(TwoSpecies)

endif ()


if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_NetPresNut")
    add_subdirectory(NetPresNut)

endif ()

if ("${LIBTG_BUILD_FLAG}" STREQUAL "Build_All" OR "${LIBTG_BUILD_FLAG}"
        STREQUAL "Build_Prashant" OR "${LIBTG_BUILD_FLAG}" STREQUAL
        "Build_NetFVFEExp")
    add_subdirectory(NetFVFEExp)
endif ()