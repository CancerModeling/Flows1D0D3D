cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(LibTG CXX)

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories("${PROJECT_SOURCE_DIR}/utils")
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/gmm")
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/aixlog")
include_directories(${VTK_INCLUDE_DIRS})

# our csv parser, included as a system dependencies to silence the warnings:
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/external/fast-cpp-csv-parser")

add_subdirectory(macrocirculation)
