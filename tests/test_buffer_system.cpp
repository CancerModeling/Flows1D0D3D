#include <mpi.h>
#include <iostream>

#include "catch2/catch.hpp"
#include "macrocirculation/communication/buffer.hpp"

namespace mc = macrocirculation;

TEST_CASE( "TestBufferSystem", "[BufferSystem]" )
{
   int rank, size;
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &size );

  if ( size < 2 )
   {
      std::cerr << "Executing this test only works with 2 processes or more!" << std::endl;
      return;
   }

   mc::BufferSystem bs( MPI_COMM_WORLD, 42 );

   // the rank to which we send
   std::size_t dst_rank = ( rank + 1 ) % size;
   // the rank from which we receive information
   std::size_t src_rank = rank > 0 ? ( rank - 1 ) % size : size - 1;

   for ( uint k = 0; k < 8; k += 1 )
   {
      bs.clear();

      // send
      {
         auto& buf = bs.get_send_buffer( dst_rank );
         buf << rank << k << ( 42. + k );
      }

      bs.start_communication();
      bs.end_communication();

      // receive
      {
         auto& buf = bs.get_receive_buffer( src_rank );

         // unpack data
         size_t other_rank, a;
         double b;
         buf >> other_rank;
         buf >> a;
         buf >> b;

         // check correctness
         REQUIRE( other_rank == src_rank );
         REQUIRE( a == k );
         REQUIRE( b == 42. + k );
      }
   }
}
