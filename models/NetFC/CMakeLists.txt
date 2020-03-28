cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "NetFC")
set(ProjectLib "Lib${ProjectName}")

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories(systems)
include_directories(util)
include_directories(inp)
include_directories(util)
include_directories(net)
include_directories(${CMAKE_SOURCE_DIR}/external/gmm)

# ****************************************************************************
# Add subprojects
# ****************************************************************************
add_subdirectory(driver)

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(util SOURCES)
AUX_SOURCE_DIRECTORY(inp SOURCES)
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(systems SOURCES)
AUX_SOURCE_DIRECTORY(net SOURCES)

add_library(${ProjectLib} ${SOURCES})
target_link_libraries(${ProjectLib} Utils ${VTK_LIBRARIES})
