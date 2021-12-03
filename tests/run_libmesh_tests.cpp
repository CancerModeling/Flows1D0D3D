#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

#include "macrocirculation/petsc/petsc_init.hpp"

namespace test_utils {

// global initializer of libMesh
std::unique_ptr<macrocirculation::PetscInit> petscInit = nullptr;

} // namespace test_utils

int main(int argc, char **argv) {
  test_utils::petscInit = std::make_unique<macrocirculation::PetscInit>(argc, argv);
  int result = Catch::Session().run(argc, argv);
  test_utils::petscInit = nullptr;
  return result;
}