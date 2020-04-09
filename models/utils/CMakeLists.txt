cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

include_directories(${CMAKE_SOURCE_DIR}/external/gmm)

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(uinp/ SOURCES)
AUX_SOURCE_DIRECTORY(umodel/ SOURCES)
AUX_SOURCE_DIRECTORY(usystem/ SOURCES)
AUX_SOURCE_DIRECTORY(unet/ SOURCES)

add_library(Utils ${SOURCES})
