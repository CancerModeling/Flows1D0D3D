//
// Created by andreas on 30.06.21.
//

#include "csv_vessel_tip_writer.hpp"

#include <fstream>
#include <nlohmann/json.hpp>

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
  std::shared_ptr<DofMap> dofmap)
    : d_comm(comm),
      d_output_directory(std::move(output_directory)),
      d_filename(std::move(filename)),
      d_graph(std::move(graph)),
      d_dofmap(std::move(dofmap)) {
  reset_all_files();
}

void CSVVesselTipWriter::write(double t, const PetscVec &u)
{
  write_generic(t, u);
}

void CSVVesselTipWriter::write(double t, const std::vector< double > &u)
{
  write_generic(t, u);
}

template< typename VectorType >
void CSVVesselTipWriter::write_generic(double t, const VectorType &u) {
  update_time(t);
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    auto v = *d_graph->get_vertex(v_id);

    if (v.is_vessel_tree_outflow() || v.is_windkessel_outflow()) {
      std::ofstream f(get_file_path(v_id), std::ios::app);
      auto local_dof_map = d_dofmap->get_local_dof_map(v);
      const auto &dofs = local_dof_map.dof_indices();
      std::vector<double> values(dofs.size());
      extract_dof(dofs, u, values);
      for (auto value : values)
        f << value << " ";
      f << std::endl;
    }
  }
}

void CSVVesselTipWriter::reset_all_files() {
  // reset all the files
  write_meta_file();
  for (auto v_id : d_graph->get_active_vertex_ids(mpi::rank(d_comm))) {
    std::ofstream f(get_file_path(v_id), std::ios::out);
    f << "";
  }
}

void CSVVesselTipWriter::write_meta_file() {
  std::ofstream f(get_meta_file_path(), std::ios::out);

  using json = nlohmann::json;

  json j;

  auto vertices_list = json::array();
  for (auto v_id : d_graph->get_vertex_ids()) {
    auto v = d_graph->get_vertex(v_id);
    auto e = d_graph->get_edge(v->get_edge_neighbors()[0]);

    if (v->is_vessel_tree_outflow() || v->is_windkessel_outflow()) {
      std::string outflow_type;
      if (v->is_vessel_tree_outflow())
        outflow_type += "vessel_tree";
      else if (v->is_windkessel_outflow())
        outflow_type += "windkessel";
      else
        outflow_type += "unknown";
      auto num_dofs = d_dofmap->get_local_dof_map(*v).num_local_dof();

      json vessel_obj = {
        {"vertex_id", v_id},
        {"name", v->get_name()},
        {"neighbor_edge_id", v->get_edge_neighbors()[0]},
        {"filepath", get_file_path(v_id)},
        {"outflow_type", outflow_type},
        {"num_dofs", num_dofs},
      };

      if (v->is_vessel_tree_outflow()) {
        vessel_obj["resistance"] = v->get_vessel_tree_data().resistances;
        vessel_obj["capacitances"] = v->get_vessel_tree_data().capacitances;
        vessel_obj["R1"] = calculate_R1(e->get_physical_data());
        vessel_obj["furcation_number"] = v->get_vessel_tree_data().furcation_number;
        vessel_obj["p_out"] = v->get_vessel_tree_data().p_out;
      }

      if (v->is_windkessel_outflow()) {
        vessel_obj["R2"] = v->get_peripheral_vessel_data().resistance - calculate_R1(e->get_physical_data());
        vessel_obj["p_out"] = v->get_peripheral_vessel_data().p_out;
        vessel_obj["C"] = v->get_peripheral_vessel_data().compliance;
      }

      vertices_list.push_back(vessel_obj);
    }
  }
  j["vertices"] = vertices_list;
  j["times"] = json::array();

  f << j.dump(1);
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

std::string CSVVesselTipWriter::get_file_path(size_t vertex_id) const {
  return d_output_directory + "/" + d_filename + "_" + std::to_string(vertex_id) + ".csv";
}

std::string CSVVesselTipWriter::get_meta_file_path() const {
  return d_output_directory + "/" + d_filename + ".json";
}

} // namespace macrocirculation
