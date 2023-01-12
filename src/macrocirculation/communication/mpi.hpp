////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_COMMUNICATION_MPI_HPP
#define TUMORMODELS_COMMUNICATION_MPI_HPP

#include <mpi.h>
#include <stdexcept>
#include <string>

namespace macrocirculation {

#define CHECK_MPI_SUCCESS(code) (check_mpi_success(code, __FILE__, __LINE__))

inline void check_mpi_success(int mpi_code, const std::string &file, int line) {
  if (mpi_code != MPI_SUCCESS)
    throw std::runtime_error("wrong mpi error code at " + file + " on line " + std::to_string(line));
}

namespace mpi {
inline int size(MPI_Comm comm) {
  int size;
  CHECK_MPI_SUCCESS(MPI_Comm_size(comm, &size));
  return size;
}

inline int rank(MPI_Comm comm) {
  int rank;
  CHECK_MPI_SUCCESS(MPI_Comm_rank(comm, &rank));
  return rank;
}

inline MPI_Comm create_local() {
  MPI_Comm new_comm;
  MPI_Comm_split(MPI_COMM_WORLD, rank(MPI_COMM_WORLD), 0, &new_comm);
  return new_comm;
} 

} // namespace mpi

} // namespace macrocirculation

#endif
