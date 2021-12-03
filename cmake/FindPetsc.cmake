# Copyright (c) 2019    Prashant K. Jha
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)
find_package(PkgConfig)

set(PETSC_DIR $ENV{PETSC_DIR} CACHE PATH "Petsc installation directory")
set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/lib/pkgconfig")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

if (NOT PETSC_LIBRARIES)
    message(FATAL_ERROR "PETSC_LIBRARIES Library not found: Specify the PETSC_DIR where petsc is located")
endif ()
