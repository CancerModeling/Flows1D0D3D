////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_CSV_VESSEL_TIP_WRITER_HPP
#define TUMORMODELS_CSV_VESSEL_TIP_WRITER_HPP

#include "mpi.h"

#include <memory>
#include <string>
#include <vector>

namespace macrocirculation {

// forward declarations:
class GraphStorage;
class DofMap;
class PetscVec;

/*! @brief Serializes the pressures at the vessel tips to csv files which are structured with a meta json-file containing the time steps and vessel information.
 *
 * The meta json file <filename>.json has the structure:
 * {
 *     "times": [ ... list of time steps ... ],
 *     "vertices": [ ...
 *        {
 *          "filepath": "<full filepath to the csv file >",
 *          "name": "<assigned name of the vertex>",
 *          "vertex_id": <the vertex id>,
 *          "neighbor_edge_id": "<the edge id of the neighboring edge, which can be observed in paraview>",
 *          "R1": <the first resistance in the model, corresponding to the 1D resistance>,
 *          ... some attributes depending on the choice of boundary condition ...
 *        }, ...
 *     ]
 * }
 *
 * The csv files form a number-of-time-steps x number-of-capacitors matrix.
 */
class CSVVesselTipWriter {
public:
  /*! @brief Constructs a vessel tip writer. Writes the json meta file and overwrites previous csv files.
   *
   * @param comm The communicator used for the parallel communication.
   * @param output_directory The directory to which we output both the csv files as well as the json meta file.
   * @param filename The filename without extension. Results in a <filename>.json file, containing meta information about the vessel tips
   *                 and <filename>_<vertex_id>_<type>.csv files for each vessel tip containing for
   *                    <type> == p: the pressures,
   *                    <type> == c: the concentrations and
   *                    <type> == v: the the volumes at each capacitor (column direction)
   *                 for each time step (row direction).
   * @param graph A graph to which the 0D models are attached.
   * @param dofmap A dof map which contains the 0D models.
   */
  CSVVesselTipWriter(
    MPI_Comm comm,
    std::string output_directory,
    std::string filename,
    std::shared_ptr<GraphStorage> graph,
    std::vector<std::shared_ptr<DofMap>> dofmaps,
    std::vector<std::string> types);

  CSVVesselTipWriter(
    MPI_Comm comm,
    std::string output_directory,
    std::string filename,
    std::shared_ptr<GraphStorage> graph,
    std::shared_ptr<DofMap> dofmaps);

  /*! @brief Writes the dofs of a 0D-Model for a Petsc vector.
   *
   * @param t The current time step to write.
   * @param u Dof-vector with pressures.
   */
  void write(double t, const PetscVec &u);

  /*! @brief Writes the dofs of a 0D-Model for a Petsc vector.
   *
   * @param t The current time step to write.
   * @param u Dof-vector with quantities.
   */
  void write(double t, const std::vector<std::reference_wrapper<const PetscVec>> &u);

  void write(double t,
             const std::vector<std::reference_wrapper<std::vector<double>>> &u_1,
             const std::vector<std::reference_wrapper<const PetscVec>> &u_2);

  /*! @brief Writes the dofs of a 0D-Model for a gmm vector.
   *
   * @param t The current time step to write.
   * @param u Dof-vector from which the degrees of freedom
   */
  void write(double t, const std::vector<double> &u);

private:
  MPI_Comm d_comm;
  std::string d_output_directory;
  std::string d_filename;
  std::shared_ptr<GraphStorage> d_graph;
  std::vector<std::shared_ptr<DofMap>> d_dofmaps;
  std::vector<std::string> d_types;

  /*! @brief Writes the empty string to all csv files and thereby deleting previously recorded data. * */
  void reset_all_files();

  /*! @brief Writes the meta json file. */
  void write_meta_file();

  /*! @brief Adds the given time to the meta json file. */
  void update_time(double t);

  /*! @brief Returns the csv filepath for a vertex with the given vertex_id for the given data type. */
  std::string get_file_path(size_t vertex_id, const std::string &type) const;

  /*! @brief Returns the csv filename for a vertex with the given vertex_id. */
  std::string get_file_name(size_t vertex_id, const std::string &type) const;

  /*! @brief Returns the json filepath for the meta file. */
  std::string get_meta_file_path() const;

  /*!
   * Generic function to write the data at the vessel tips for all kinds of vectors.
   *
   * @param t The current time step to write.
   * @param u Dof-vector from which the degrees of freedom
   */
  template<typename VectorType>
  void write_generic(const DofMap &dof_map, const VectorType &u, const std::string &type);

  void write_p_out();
};

} // namespace macrocirculation

#endif //TUMORMODELS_CSV_VESSEL_TIP_WRITER_HPP
