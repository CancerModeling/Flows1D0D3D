//
// Created by andreas on 30.06.21.
//

#include "csv_vessel_tip_writer.hpp"

#include <fstream>
#include <nlohmann/json.hpp>
#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "graph_storage.hpp"
#include "petsc/petsc_vec.hpp"
#include "vessel_formulas.hpp"

namespace macrocirculation {

CSVVesselTipWriter::CSVVesselTipWriter(
  MPI_Comm comm,
  std::string output_directory,
  std::string filename,
  std::shared_ptr<GraphStorage> graph,
  std::vector<std::shared_ptr<DofMap>> dofmaps,
  std::vector<std::string> types)
    : d_comm(comm),
      d_output_directory(std::move(output_directory)),
      d_filename(std::move(filename)),
      d_graph(std::move(graph)),
      d_dofmaps(std::move(dofmaps)),
      d_types(std::move(types)) {
  if (d_dofmaps.size() != d_types.size())
    throw std::runtime_error("number of types and dofmaps must coincide");

  reset_all_files();
}

CSVVesselTipWriter::CSVVesselTipWriter(MPI_Comm comm,
                                       std::string output_directory,
                                       std::string filename,
                                       std::shared_ptr<GraphStorage> graph,
                                       std::shared_ptr<DofMap> dofmaps)
    : CSVVesselTipWriter(comm,
                         std::move(output_directory),
                         std::move(filename),
                         std::move(graph),
                         std::vector<std::shared_ptr<DofMap>>({std::move(dofmaps)}),
                         std::vector<std::string>({"p"})) {}

void CSVVesselTipWriter::write(double t, const PetscVec &u) {
  if (d_dofmaps.size() > 1)
    throw std::runtime_error("CSVVesselTipWriter::write: method supported only for one substance");

  update_time(t);
  write_p_out();
  write_generic(*d_dofmaps.front(), u, d_types.front());
}

void CSVVesselTipWriter::write(double t, const std::vector<double> &u) {
  if (d_dofmaps.size() > 1)
    throw std::runtime_error("CSVVesselTipWriter::write: method supported only for one substance");

  update_time(t);
  write_p_out();
  write_generic(*d_dofmaps.front(), u, d_types.front());
}

void CSVVesselTipWriter::write(double t,
                               const std::vector<std::reference_wrapper<std::vector<double>>> &u_1,
                               const std::vector<std::reference_wrapper<const PetscVec>> &u_2) {
  if (d_dofmaps.size() != u_1.size() + u_2.size())
    throw std::runtime_error("CSVVesselTipWriter::write: not all expected quantities provided");

  update_time(t);
  write_p_out();
  for (size_t k = 0; k < u_1.size(); k += 1)
    write_generic(*d_dofmaps[k], u_1[k], d_types[k]);
  for (size_t k = 0; k < u_2.size(); k += 1)
    write_generic(*d_dofmaps[u_1.size() + k], u_2[k], d_types[u_1.size() + k]);
}

void CSVVesselTipWriter::write(double t, const std::vector<std::reference_wrapper<const PetscVec>> &quantities) {
  if (d_dofmaps.size() != quantities.size())
    throw std::runtime_error("CSVVesselTipWriter::write: not all expected quantities provided");

  update_time(t);
  write_p_out();
  for (size_t k = 0; k < quantities.size(); k += 1)
    write_generic(*d_dofmaps[k], quantities[k], d_types[k]);
}

template<typename VectorType>
void CSVVesselTipWriter::write_generic(const DofMap &dof_map, const VectorType &u, const std::string &type) {
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto v = *d_graph->get_vertex(v_id);

    if (v.is_vessel_tree_outflow() || v.is_windkessel_outflow() || v.is_rcl_outflow()) {
      std::ofstream f(get_file_path(v_id, type), std::ios::app);
      auto local_dof_map = dof_map.get_local_dof_map(v);
      const auto &dofs = local_dof_map.dof_indices();
      std::vector<double> values(dofs.size());
      extract_dof(dofs, u, values);
      for (auto value : values)
        f << value << " ";
      f << std::endl;
    }
  }
}

void CSVVesselTipWriter::write_p_out() {
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto v = *d_graph->get_vertex(v_id);

    double p_out = 0;
    if (v.is_windkessel_outflow())
      p_out = v.get_peripheral_vessel_data().p_out;
    else if (v.is_vessel_tree_outflow())
      p_out = v.get_vessel_tree_data().p_out;
    else
      continue;

    if (v.is_vessel_tree_outflow() || v.is_windkessel_outflow() || v.is_rcl_outflow()) {
      std::ofstream f(get_file_path(v_id, "p_out"), std::ios::app);
      f << p_out << " " << std::endl;
    }
  }
}

void CSVVesselTipWriter::reset_all_files() {
  // reset all the files
  write_meta_file();
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto &vertex = *d_graph->get_vertex(v_id);
    for (auto &type : d_types) {
      std::ofstream f(get_file_path(v_id, type), std::ios::out);
      f << "";
    }
    // write special file for p_out
    if (vertex.is_vessel_tree_outflow() || vertex.is_windkessel_outflow()) {
      std::ofstream f(get_file_path(v_id, "p_out"), std::ios::out);
      f << "";
    }
  }
}

void CSVVesselTipWriter::write_meta_file() {
  using json = nlohmann::json;

  json j;

  auto vertices_list = json::array();
  for (auto v_id : d_graph->get_vertex_ids()) {
    auto v = d_graph->get_vertex(v_id);
    auto e = d_graph->get_edge(v->get_edge_neighbors()[0]);

    if (v->is_vessel_tree_outflow() || v->is_windkessel_outflow() || v->is_rcl_outflow()) {
      std::string outflow_type;
      if (v->is_vessel_tree_outflow())
        outflow_type += "vessel_tree";
      else if (v->is_windkessel_outflow())
        outflow_type += "windkessel";
      else if (v->is_rcl_outflow())
        outflow_type += "rcl";
      else
        outflow_type += "unknown";

      // get number of dofs:
      int num_dofs = mpi::rank(d_comm) == e->rank() ? static_cast<int>(d_dofmaps.front()->get_local_dof_map(*v).num_local_dof()) : 0;
      // only the master process has to know the number of dof.
      if (e->rank() != 0) {
        if (mpi::rank(d_comm) == e->rank())
          MPI_Send(&num_dofs, 1, MPI_INT, 0, 101, d_comm);
        if (mpi::rank(d_comm) == 0)
          MPI_Recv(&num_dofs, 1, MPI_INT, static_cast<int>(e->rank()), 101, d_comm, nullptr);
      }

      std::map<std::string, std::string> filepaths;
      for (size_t k = 0; k < d_types.size(); k += 1) {
        filepaths["filepath_" + d_types[k]] = get_file_name(v_id, d_types[k]);
      }

      json vessel_obj = {
        {"vertex_id", v_id},
        {"name", v->get_name()},
        {"neighbor_edge_id", v->get_edge_neighbors()[0]},
        {"filepaths", filepaths},
        {"outflow_type", outflow_type},
        {"num_dofs", num_dofs},
      };

      if (v->is_vessel_tree_outflow()) {
        vessel_obj["resistance"] = v->get_vessel_tree_data().resistances;
        vessel_obj["capacitances"] = v->get_vessel_tree_data().capacitances;
        vessel_obj["radii"] = v->get_vessel_tree_data().radii;
        vessel_obj["R1"] = calculate_R1(e->get_physical_data());
        vessel_obj["furcation_number"] = v->get_vessel_tree_data().furcation_number;
        vessel_obj["filepath_p_out"] = get_file_name(v_id, "p_out");
      }

      if (v->is_windkessel_outflow()) {
        vessel_obj["R2"] = v->get_peripheral_vessel_data().resistance - calculate_R1(e->get_physical_data());
        vessel_obj["C"] = v->get_peripheral_vessel_data().compliance;
        vessel_obj["filepath_p_out"] = get_file_name(v_id, "p_out");
      }

      if (v->is_rcl_outflow()) {
        vessel_obj["R"] = v->get_rcl_data().resistances;
        vessel_obj["C"] = v->get_rcl_data().capacitances;
        vessel_obj["L"] = v->get_rcl_data().inductances;
      }

      vertices_list.push_back(vessel_obj);
    }
  }
  j["vertices"] = vertices_list;
  j["times"] = json::array();

  if (mpi::rank(d_comm) == 0) {
    std::ofstream f(get_meta_file_path(), std::ios::out);
    f << j.dump(1);
  }
}

void CSVVesselTipWriter::update_time(double t) {
  if (mpi::rank(d_comm) != 0)
    return;

  nlohmann::json j;
  {
    std::ifstream f(get_meta_file_path());
    f >> j;
  }

  j["times"].push_back(t);

  {
    std::ofstream f(get_meta_file_path(), std::ios::out);
    f << j.dump(1);
  }
}

std::string CSVVesselTipWriter::get_file_path(size_t vertex_id, const std::string &type) const {
  return d_output_directory + "/" + get_file_name(vertex_id, type);
}

std::string CSVVesselTipWriter::get_file_name(size_t vertex_id, const std::string &type) const {
  return d_filename + "_" + std::to_string(vertex_id) + "_" + type + ".csv";
}

std::string CSVVesselTipWriter::get_meta_file_path() const {
  return d_output_directory + "/" + d_filename + ".json";
}

} // namespace macrocirculation
