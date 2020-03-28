cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)

add_library(TestGeom ${SOURCES})

