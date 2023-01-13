////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "./buffer.hpp"

#include "./mpi.hpp"

namespace macrocirculation {

SendBuffer::SendBuffer()
    : d_buffer(0) {}

std::vector<char> &SendBuffer::get_buffer() { return d_buffer; }

void SendBuffer::clear() { d_buffer.clear(); }

char *SendBuffer::ptr() { return &d_buffer[0]; }

std::size_t SendBuffer::size() { return d_buffer.size(); }

ReceiveBuffer::ReceiveBuffer(std::vector<char> buffer)
    : d_buffer(std::move(buffer)), d_current_position(0) {}

std::vector<char> &ReceiveBuffer::get_buffer() { return d_buffer; }

void ReceiveBuffer::clear() {
  d_buffer.clear();
  d_current_position = 0;
}

char *ReceiveBuffer::ptr() { return &d_buffer[0]; }

std::size_t ReceiveBuffer::size() { return d_buffer.size(); }

bool ReceiveBuffer::empty() { return d_buffer.size() == d_current_position; }

BufferSystem::BufferSystem(MPI_Comm comm, std::size_t tag)
    : d_num_processes(mpi::size(comm)),
      d_rank(mpi::rank(comm)),
      d_tag(tag),
      d_receive_buffers(d_num_processes),
      d_send_buffers(d_num_processes),
      d_send_requests(d_num_processes),
      d_receive_requests(d_num_processes) {}

ReceiveBuffer &BufferSystem::get_receive_buffer(std::size_t from) {
  return d_receive_buffers.at(from);
}

SendBuffer &BufferSystem::get_send_buffer(std::size_t to) {
  return d_send_buffers.at(to);
};

void BufferSystem::clear() {
  for (std::size_t k = 0; k < d_num_processes; k += 1) {
    // check that the receive buffer is empty now and no information is unused.
    assert(d_receive_buffers[k].empty());
    d_receive_buffers[k].clear();
    d_send_buffers[k].clear();
  }
}

void BufferSystem::start_communication() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // send the data in the send buffers:
  for (std::size_t recipient = 0; recipient < d_send_buffers.size(); recipient += 1) {
    if (static_cast<size_t>(rank) == recipient)
      continue;

    auto &buffer = get_send_buffer(recipient);
    CHECK_MPI_SUCCESS(
      MPI_Isend(buffer.ptr(), buffer.size(), MPI_BYTE, recipient, d_tag, MPI_COMM_WORLD, &d_send_requests[recipient]));
  }

  // resize the receive buffers, depending on how much data arrives:
  for (std::size_t sender = 0; sender < d_receive_buffers.size(); sender += 1) {
    if (static_cast<size_t>(rank) == sender)
      continue;

    MPI_Status status;
    CHECK_MPI_SUCCESS(MPI_Probe(sender, d_tag, MPI_COMM_WORLD, &status));

    int size;
    CHECK_MPI_SUCCESS(MPI_Get_count(&status, MPI_BYTE, &size));

    auto &buffer = get_receive_buffer(sender);
    buffer.get_buffer().resize(size);
  }

  // initiate receiving:
  for (std::size_t sender = 0; sender < d_receive_buffers.size(); sender += 1) {
    if (static_cast<size_t>(rank) == sender)
      continue;

    auto &buffer = get_receive_buffer(sender);
    CHECK_MPI_SUCCESS(
      MPI_Irecv(buffer.ptr(), buffer.size(), MPI_BYTE, sender, d_tag, MPI_COMM_WORLD, &d_receive_requests[sender]));
  }

  // the data we send to ourselves is just copied:
  {
    auto &sbuf = get_send_buffer(rank);
    auto &rbuf = get_receive_buffer(rank);
    rbuf.get_buffer().resize(sbuf.get_buffer().size());
    for (size_t k = 0; k < sbuf.size(); k += 1)
      rbuf.get_buffer()[k] = sbuf.get_buffer()[k];
  }
}

void BufferSystem::end_communication() {
  assert(d_send_requests.size() == d_receive_requests.size());

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (std::size_t recipient = 0; recipient < d_send_buffers.size(); recipient += 1) {
    if (static_cast<size_t>(rank) == recipient)
      continue;

    MPI_Status status;
    CHECK_MPI_SUCCESS(MPI_Wait(&d_send_requests[recipient], &status));
  }

  for (std::size_t sender = 0; sender < d_send_requests.size(); sender += 1) {
    if (static_cast<size_t>(rank) == sender)
      continue;

    MPI_Status status;
    CHECK_MPI_SUCCESS(MPI_Wait(&d_receive_requests[sender], &status));
  }
}
} // namespace macrocirculation
