#ifndef TUMORMODELS_MATLAB_QOI_WRITER_HPP
#define TUMORMODELS_MATLAB_QOI_WRITER_HPP

#include "libmesh/libmesh.h"
#include "libmesh/parallel.h"
#include <fstream>

namespace util {

class QoIVec;

class MatlabQoIWriter {
public:
  MatlabQoIWriter(libMesh::Parallel::Communicator *comm, std::string filename);

  void write(const util::QoIVec &qoi) const;

private:
  std::string d_filename;
  libMesh::Parallel::Communicator *d_comm;
};

} // namespace util

#endif //TUMORMODELS_MATLAB_QOI_WRITER_HPP
