#include <mpi.h>

#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

int main( int argc, char* argv[] )
{
   MPI_Init( &argc, &argv );
   int result = Catch::Session().run( argc, argv );
   MPI_Finalize();
   return result;
}
