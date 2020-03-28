#(C) https://github.com/capitalaslash/cmake-modules/blob/master/FindLibmesh.cmake
# (BSD - License)

# - Try to find Libmesh
# Once done this will define
#
#  LIBMESH_FOUND          - Libmesh has been successfully found
#  LIBMESH_INCLUDE_DIRS   - Libmesh include directories
#  LIBMESH_LIBRARIES      - Libmesh libraries
#  LIBMESH_DEFINITIONS    - Libmesh definitions
#  LIBMESH_FLAGS          - Libmesh flags
#  LIBMESH_VERSION_STRING - Libmesh version
#
#  Usage:
#  find_package(Libmesh)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE NONE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER)
message("BUILD_TYPE_UPPER = " ${BUILD_TYPE_UPPER})
if(${BUILD_TYPE_UPPER} MATCHES DEBUG)
  set(METHOD dbg)
else()
  set(METHOD opt)
endif()
#set(METHOD opt)
message(STATUS "linking against ${METHOD} libmesh library")

# required for LIBMESH_DEFINITIONS
find_package(PkgConfig REQUIRED)

set(LIBMESH_DIR LIBMESH_DIR-NOTFOUND CACHE PATH "Libmesh installation directory")

if(LIBMESH_DIR)
  set(ENV{PKG_CONFIG_PATH} "${LIBMESH_DIR}/lib/pkgconfig")
#  message($ENV{PKG_CONFIG_PATH})
endif()

pkg_check_modules(PC_LIBMESH libmesh-${METHOD})
message(STATUS "PC_LIBMESH_FOUND = ${PC_LIBMESH_FOUND}")
message(STATUS "PC_LIBMESH_LIBRARIES = ${PC_LIBMESH_LIBRARIES}")
message(STATUS "PC_LIBMESH_LIBRARY_DIRS = ${PC_LIBMESH_LIBRARY_DIRS}")

# distinguish flags and definitions (-D...)
#message(STATUS "PC_LIBMESH_CFLAGS_OTHER = ${PC_LIBMESH_CFLAGS_OTHER}")
foreach(FLAG ${PC_LIBMESH_CFLAGS_OTHER})
  if(${FLAG} MATCHES "^[-][D].+")
    set(PC_LIBMESH_CFLAGS_DEFS "${PC_LIBMESH_CFLAGS_DEFS} ${FLAG}")
  else()
    set(PC_LIBMESH_CFLAGS_FLAGS "${PC_LIBMESH_CFLAGS_FLAGS} ${FLAG}")
  endif()
endforeach()
set(LIBMESH_DEFINITIONS ${PC_LIBMESH_CFLAGS_DEFS})
set(LIBMESH_FLAGS ${PC_LIBMESH_CFLAGS_FLAGS})

find_path(LIBMESH_INCLUDE_DIR libmesh/libmesh.h
  HINTS ${PC_LIBMESH_INCLUDEDIR} ${PC_LIBMESH_INCLUDE_DIRS}
  PATH_SUFFIXES libmesh
)

find_library(LIBMESH_LIBRARY
  NAMES mesh_${METHOD} libmesh
  HINTS ${PC_LIBMESH_LIBDIR} ${PC_LIBMESH_LIBARY_DIRS}
)

set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARY})
set(LIBMESH_INCLUDE_DIRS ${LIBMESH_INCLUDE_DIR})

find_program( LIBMESH_CONFIG_EXECUTABLE
    NAMES libmesh-config
    HINTS ${LIBMESH_DIR} $ENV{LIBMESH_DIR}
    PATH_SUFFIXES bin
    DOC "libmesh-config executable"
)
mark_as_advanced( LIBMESH_CONFIG_EXECUTABLE )

exec_program( ${LIBMESH_CONFIG_EXECUTABLE}
  ARGS --include
  OUTPUT_VARIABLE LMC_INC_FLAG
  RETURN_VALUE LMC_INC_RET
)
string(REPLACE " " ";" LMC_INC_LIST ${LMC_INC_FLAG})
foreach( IPATH ${LMC_INC_LIST} )
  string(REGEX REPLACE "^-I" "" IPATH ${IPATH})
  string(REGEX REPLACE "//" "/" IPATH ${IPATH})
  list(APPEND LM_INC ${IPATH})
endforeach()
set(LIBMESH_INCLUDE_DIRS ${LM_INC})

if(PC_LIBMESH_VERSION)
  set(LIBMESH_VERSION_STRING ${PC_LIBMESH_VERSION})
endif()

# handle the QUIETLY and REQUIRED arguments and set LIBMESH_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libmesh
  FOUND_VAR LIBMESH_FOUND
  REQUIRED_VARS LIBMESH_LIBRARIES LIBMESH_INCLUDE_DIRS LIBMESH_DEFINITIONS
    LIBMESH_FLAGS
  VERSION_VAR LIBMESH_VERSION_STRING
)

mark_as_advanced(
  LIBMESH_INCLUDE_DIR
  LIBMESH_LIBRARIES
  LIBMESH_CONFIG_EXECUTABLE
)
