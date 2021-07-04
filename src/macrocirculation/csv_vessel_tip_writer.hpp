////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_CSV_VESSEL_TIP_WRITER_HPP
#define TUMORMODELS_CSV_VESSEL_TIP_WRITER_HPP

#include "mpi.h"

#include <string>
#include <memory>

namespace macrocirculation {

// forward declarations:
class GraphStorage;
class DofMap;
class PetscVec;

class CSVVesselTipWriter {
public:
  CSVVesselTipWriter(
    MPI_Comm comm,
    std::string output_directory,
    std::string filename,
    std::shared_ptr<GraphStorage> graph,
    std::shared_ptr<DofMap> dofmap);

  void write(double t, const PetscVec &u);

private:
  MPI_Comm d_comm;
  std::string d_output_directory;
  std::string d_filename;
  std::shared_ptr<GraphStorage> d_graph;
  std::shared_ptr<DofMap> d_dofmap;

  void reset_all_files();

  void write_meta_file();

  void update_time(double t);

  std::string get_file_path(size_t vertex_id) const;

  std::string get_meta_file_path() const;
};

}

#endif //TUMORMODELS_CSV_VESSEL_TIP_WRITER_HPP
