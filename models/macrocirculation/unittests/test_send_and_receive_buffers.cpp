#include "catch2/catch.hpp"

#include "macrocirculation/communication/buffer.hpp"

namespace mc = macrocirculation;

TEST_CASE("TestSendAndRecevieBuffers", "[send_and_receive_buffer]") {
   mc::SendBuffer sb;
   sb << 4.1 << 5U << 42.3;
   mc::ReceiveBuffer rb(sb.get_buffer());

   double a;
   rb >> a;
   REQUIRE(a == 4.1);

   unsigned int b;
   rb >> b;
   REQUIRE(b == 5U);

   double c;
   rb >> c;
   REQUIRE(c == 42.3);
}
