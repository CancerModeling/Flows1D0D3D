# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectLib "LibMacrocirculation")

# ****************************************************************************
# Create library
# ****************************************************************************
AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(./macrocirculation SOURCES)
AUX_SOURCE_DIRECTORY(./macrocirculation/communication SOURCES)

add_library(${ProjectLib} ${SOURCES})
target_link_libraries(${ProjectLib} PUBLIC ${LIBMESH_LIBRARIES}
                                           ${PETSC_LIBRARIES}
                                           ${VTK_LIBRARIES}
                                           ${MPI_CXX_LIBRARIES}
                                           ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(${ProjectLib} INTERFACE ./)
target_link_libraries(${ProjectLib} PRIVATE nlohmann_json::nlohmann_json)
