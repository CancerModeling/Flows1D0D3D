#include "csv_qoi_writer.hpp"

#include "qoi.hpp"

namespace util {

CSVQoIWriter::CSVQoIWriter(libMesh::Parallel::Communicator *comm, std::string filename)
    : d_comm(comm),
      d_filename(std::move(filename)) {}

void CSVQoIWriter::write(const util::QoIVec &qoi) const {
  // only rank zero is allowed to write
  if (d_comm->rank() != 0)
    return;

  std::fstream filecsv;
  filecsv.open(d_filename, std::ios::out);

  const auto names = qoi.get_names();
  const auto data = qoi.get_all();

  const auto num_timesteps = data.size();
  const auto num_quantities = names.size();

  // write header of csv file
  for (int qidx = 0; qidx < num_quantities; qidx += 1) {
    filecsv << names[qidx];
    if (qidx < num_quantities - 1)
      filecsv << ",";
    else
      filecsv << "\n";
  }

  for (int tidx = 0; tidx < num_timesteps; tidx += 1) {
    for (int qidx = 0; qidx < num_quantities; qidx += 1) {
      filecsv << std::scientific << data[tidx][qidx];
      if (qidx < num_quantities - 1)
        filecsv << ",";
      else
        filecsv << "\n";
    }
  }
}

} // namespace util
