////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_BUFFER_HPP
#define TUMORMODELS_BUFFER_HPP

#include <cassert>
#include <mpi.h>
#include <vector>

namespace macrocirculation {

class SendBuffer {
public:
  SendBuffer();

  template<typename AtomicType>
  SendBuffer &operator<<(const AtomicType &value);

  std::vector<char> &get_buffer();

  void clear();

  char *ptr();

  std::size_t size();

private:
  std::vector<char> d_buffer;
};

class ReceiveBuffer {
public:
  ReceiveBuffer() = default;

  explicit ReceiveBuffer(std::vector<char> buffer);

  template<typename AtomicType>
  ReceiveBuffer &operator>>(AtomicType &value);

  std::vector<char> &get_buffer();

  void clear();

  char *ptr();

  std::size_t size();

  bool empty();

private:
  std::vector<char> d_buffer;
  std::size_t d_current_position{};
};

class BufferSystem {
public:
  explicit BufferSystem(MPI_Comm comm, std::size_t tag);

  ReceiveBuffer &get_receive_buffer(std::size_t from);

  SendBuffer &get_send_buffer(std::size_t to);

  void clear();

  void start_communication();

  void end_communication();

private:
  std::size_t d_num_processes;
  std::size_t d_rank;
  std::size_t d_tag;

  std::vector<ReceiveBuffer> d_receive_buffers;
  std::vector<SendBuffer> d_send_buffers;

  std::vector<MPI_Request> d_send_requests;
  std::vector<MPI_Request> d_receive_requests;
};

// implementations of template functions:
template<typename AtomicType>
SendBuffer &SendBuffer::operator<<(const AtomicType &value) {
  const char *serialized = reinterpret_cast<const char *>(&value);
  for (std::size_t k = 0; k < sizeof(AtomicType); k += 1)
    d_buffer.push_back(serialized[k]);
  return *this;
}

template<typename AtomicType>
inline ReceiveBuffer &ReceiveBuffer::operator>>(AtomicType &value) {
  assert(d_current_position + sizeof(AtomicType) <= d_buffer.size());
  value = *(reinterpret_cast<AtomicType *>(&d_buffer[d_current_position]));
  d_current_position += sizeof(AtomicType);
  return *this;
}

} // namespace macrocirculation

#endif
