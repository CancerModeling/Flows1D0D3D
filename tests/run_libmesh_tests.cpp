#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

#include "libmesh/libmesh.h"

namespace test_utils {

// global initializer of libMesh
std::unique_ptr<libMesh::LibMeshInit> lmInit = nullptr;

} // namespace test_utils

int main(int argc, char **argv) {
  test_utils::lmInit = std::make_unique<libMesh::LibMeshInit>(argc, argv);
  int result = Catch::Session().run(argc, argv);
  test_utils::lmInit = nullptr;
  return result;
}