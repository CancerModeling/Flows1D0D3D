#ifndef TUMORMODELS_CSV_QOI_WRITER_HPP
#define TUMORMODELS_CSV_QOI_WRITER_HPP

#include "libmesh/libmesh.h"
#include "libmesh/parallel.h"
#include <fstream>

namespace util {

class QoIVec;

class CSVQoIWriter {
public:
  CSVQoIWriter(libMesh::Parallel::Communicator *comm, std::string filename);

  void write(const util::QoIVec &qoi) const;

private:
  std::string d_filename;
  libMesh::Parallel::Communicator *d_comm;
};

} // namespace util

#endif //TUMORMODELS_CSV_QOI_WRITER_HPP
