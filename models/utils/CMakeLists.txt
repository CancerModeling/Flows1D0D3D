cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(uinp/ SOURCES)
AUX_SOURCE_DIRECTORY(umodel/ SOURCES)
AUX_SOURCE_DIRECTORY(usystem/ SOURCES)
AUX_SOURCE_DIRECTORY(unet/ SOURCES)
AUX_SOURCE_DIRECTORY(rw/ SOURCES)

add_library(Utils ${SOURCES})

target_link_libraries(Utils
        ${LIBMESH_LIBRARIES} ${PETSC_LIBRARIES} ${VTK_LIBRARIES} ${MPI_CXX_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT})

INSTALL(TARGETS Utils
        EXPORT Utils
        RUNTIME DESTINATION "bin"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib")
INSTALL(EXPORT Utils DESTINATION "lib")