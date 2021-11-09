#find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE src/main.cpp)
target_link_libraries(_core PUBLIC LibMacrocirculation)

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

# see: https://github.com/pybind/cmake_example/issues/11
set_target_properties(_core PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
set_target_properties(_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS _core DESTINATION .)
