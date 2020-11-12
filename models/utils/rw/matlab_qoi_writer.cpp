#include "matlab_qoi_writer.hpp"

#include "qoi.hpp"

namespace util {

MatlabQoIWriter::MatlabQoIWriter(libMesh::Parallel::Communicator *comm, std::string filename)
    : d_comm(comm),
      d_filename(std::move(filename)) {}

void MatlabQoIWriter::write(const util::QoIVec &qoi) const {
  // only rank zero is allowed to write
  if (d_comm->rank() != 0)
    return;

  std::fstream filematlab;
  filematlab.open(d_filename, std::ios::out);

  const auto names = qoi.get_names();
  const auto data = qoi.get_all();

  const auto num_timesteps = data.size();
  const auto num_quantities = names.size();

  for (int qidx = 0; qidx < num_quantities; qidx += 1) {
    filematlab << names[qidx] << " = [ ";
    for (int tidx = 0; tidx < num_timesteps; tidx += 1) {
      filematlab << std::scientific << data[tidx][qidx] << " ";
    }
    filematlab << "];" << std::endl;
  }
}

} // namespace util
