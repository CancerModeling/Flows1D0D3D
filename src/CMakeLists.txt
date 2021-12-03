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
AUX_SOURCE_DIRECTORY(./macrocirculation/petsc SOURCES)

add_library(${ProjectLib} ${SOURCES})
set_target_properties(${ProjectLib} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(${ProjectLib} PUBLIC ${MPI_CXX_LIBRARIES}
                                           ${CMAKE_THREAD_LIBS_INIT}
                                           PkgConfig::PETSC
                                           Eigen3::Eigen)
target_include_directories(${ProjectLib} INTERFACE ./)
target_link_libraries(${ProjectLib} PRIVATE nlohmann_json::nlohmann_json)
