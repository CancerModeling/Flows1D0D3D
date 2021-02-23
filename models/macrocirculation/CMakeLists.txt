cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "Macrocirculation")
set(ProjectLib "Lib${ProjectName}")

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(src/macrocirculation SOURCES)
AUX_SOURCE_DIRECTORY(src/macrocirculation/communication SOURCES)

add_library(${ProjectLib} ${SOURCES})
target_link_libraries(${ProjectLib} Utils)
target_include_directories(${ProjectLib} INTERFACE ./src)

# ****************************************************************************
# Add apps
# ****************************************************************************
add_subdirectory(apps)
