cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "NetFVFE")
set(ProjectLib "Lib${ProjectName}")

# ****************************************************************************
# Include directories
# ****************************************************************************
include_directories(systems)

# ****************************************************************************
# Add subprojects
# ****************************************************************************
add_subdirectory(driver)

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(systems SOURCES)

add_library(${ProjectLib} ${SOURCES})
target_link_libraries(${ProjectLib} Utils)
