cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
project(Applications CXX)

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories("${CMAKE_SOURCE_DIR}/models/utils")
include_directories("${CMAKE_SOURCE_DIR}/models/")
include_directories("${CMAKE_SOURCE_DIR}/external/gmm")
include_directories("${CMAKE_SOURCE_DIR}/external/aixlog")
include_directories(${VTK_INCLUDE_DIRS})

# ****************************************************************************
# Add subprojects
# ****************************************************************************
add_subdirectory(sensitivity)
