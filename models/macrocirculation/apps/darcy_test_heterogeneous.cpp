////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <cxxopts.hpp>
#include <memory>

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "macrocirculation/assembly_system.hpp"
#include "macrocirculation/base_model.hpp"
#include "macrocirculation/vtk_io.hpp"
#include "macrocirculation/nifti_reader.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedIntArray.h>

namespace mc = macrocirculation;

// define input, systems, model
namespace darcy3d {
// read and store input parameters
struct InputDeck {
  InputDeck(const std::string &filename = "") : d_K_vasc(1.), d_K_tiss(0.2), d_T(1.),
                                                d_dt(0.01), d_h(0.1), d_mesh_file("") {
    if (!filename.empty() and filename != "")
      read_parameters(filename);
  }

  void read_parameters(const std::string &filename) {
    GetPot input(filename);
    d_K_vasc = input("K_vasc", 1.);
    d_K_tiss = input("K_tiss", 0.2);
    d_T = input("T", 1.);
    d_dt = input("dt", 0.01);
    d_h = input("h", 0.1);
    d_mesh_file = input("mesh_file", "mesh.e");
  }

  std::string print_str() {
    std::ostringstream oss;
    oss << "K_vasc = " << d_K_vasc << "\n";
    oss << "K_tiss = " << d_K_tiss << "\n";
    oss << "T = " << d_T << "\n";
    oss << "dt = " << d_dt << "\n";
    oss << "h = " << d_h << "\n";
    oss << "mesh_file = " << d_mesh_file << "\n";
    return oss.str();
  }

  double d_K_vasc;
  double d_K_tiss;
  double d_T;
  double d_dt;
  double d_h;
  std::string d_mesh_file;
};

// forward declare model class (specific definition requires first defining systems)
class Model;

// define pressure system

// ic
lm::Number ic_p(const lm::Point &p, const lm::Parameters &es,
                const std::string &system_name,
                const std::string &var_name) { return 0.; }

void ic(lm::EquationSystems &es, const std::string &system_name) {
  auto &sys = es.get_system<lm::TransientLinearImplicitSystem>(
    system_name);
  if (system_name == "P")
    sys.project_solution(ic_p, nullptr,
                         es.parameters);
}

// assembly
class Pres : public mc::BaseAssembly {
public:
  Pres(Model *model, lm::MeshBase &mesh,
       lm::TransientLinearImplicitSystem &sys)
      : mc::BaseAssembly("P", mesh, sys, 1,
                         {sys.variable_number("p")}),
        d_model_p(model) {
    sys.attach_assemble_object(
      *this);                     // attach this element assembly object
    sys.attach_init_function(ic); // add ic
  }

  void assemble() override;
  Model *d_model_p;
};

// complete definition of model
class Model : public mc::BaseModel {
public:
  Model(lm::Parallel::Communicator *comm,
        InputDeck &input,
        lm::ReplicatedMesh &mesh,
        lm::EquationSystems &eq_sys,
        lm::TransientLinearImplicitSystem &pres,
        lm::ExplicitSystem &hyd_cond,
        mc::Logger &log)
      : mc::BaseModel(comm, mesh, eq_sys, log, "Darcy_3D"),
        d_input(input),
        d_pres(this, d_mesh, pres),
        d_hyd_cond(hyd_cond) {
    d_log("model created\n");
  };

  Pres &get_pres() { return d_pres; }

  void run() override{};

  void write_system(const unsigned int &t_step) override{};

  void solve_system() override{};

  void compute_qoi() override{};

  InputDeck &d_input;
  Pres d_pres;
  lm::ExplicitSystem &d_hyd_cond;
};

// read vascular domain from nifti file and create hydraulic parameter with one value in
// vascular domain and other outside the vascular domain
void create_heterogeneous_conductivity(std::string vasc_filename, double voxel_size, lm::MeshBase &mesh, lm::ExplicitSystem &hyd_cond, lm::EquationSystems &eq_sys) {

  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");

  // read file
  auto vasc_nifti = mc::NiftiReader(vasc_filename);
  auto img_dim = vasc_nifti.get_data_dimension();
  auto img_fields = vasc_nifti.get_point_fields();
  std::vector<double> img_data;
  vasc_nifti.read_point_data(img_fields[0], &img_data);
  std::vector<std::vector<std::vector<double>>> img_data_grid;
  vasc_nifti.read_point_data(img_fields[0], &img_data_grid);
  std::cout << "nifit data \n\n";
  std::cout << vasc_nifti.print_str() << "\n\n";
  size_t count = 0;
  for (auto a : img_data)
    if (a > 0.5)
      count++;
  std::cout << "\nvascular voxel count = " << count << "\n\n";

  auto vtu_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  vtu_writer->SetFileName("vascular_domain.vtu");
  auto points = vtkSmartPointer<vtkPoints>::New();

  auto vtu_writer2 = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  vtu_writer2->SetFileName("extravascular_domain.vtu");
  auto points2 = vtkSmartPointer<vtkPoints>::New();

  std::vector<unsigned int> dof_indices;
  for (const auto &elem : mesh.active_local_element_ptr_range()) {
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);
    auto x = elem->centroid() / voxel_size; // + lm::Point(0.5, 0.5, 0.5);

    // find voxel element (1d vector representation of the image data)
    int i = mc::locate_voxel_1d({x(0), x(1), x(2)}, img_dim);
    auto i_3d = mc::locate_voxel_3d({x(0), x(1), x(2)}, img_dim);

    // get image data
    auto a = img_data_grid[i_3d[0]][i_3d[1]][i_3d[2]];

    // debug
    // auto a = img_data[i + 352];
    //    std::cout << "center = " << elem->centroid()
    //              << ", elem id = " << elem->id()
    //              << ", dof = " << dof_indices[0]
    //              << ", x = " << x
    //              << ", i = " << i
    //              << ", i_3d = (" << i_3d[0] << ", " << i_3d[1] << ", " << i_3d[2]
    //              << "), i_check = " << mc::index_3d_1d(i_3d, img_dim)
    //              << ", data = " << a << "\n";
    auto xv = lm::Point(i_3d[0], i_3d[1], i_3d[2]) * voxel_size;
    if (a > 0.5) {
      //std::cout << "vascular domain: xv = " << xv << ", elem center = " << elem->centroid() << "\n";
      points->InsertNextPoint(xv(0), xv(1), xv(2));
    } else {
      points2->InsertNextPoint(xv(0), xv(1), xv(2));
    }

    // set parameter
    if (a > 0.5)
      hyd_cond.solution->set(dof_indices[0], input->d_K_vasc);
    else
      hyd_cond.solution->set(dof_indices[0], input->d_K_tiss);
  }
  hyd_cond.solution->close();
  hyd_cond.update();

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);
  vtu_writer->SetInputData(grid);
  vtu_writer->SetDataModeToAppended();
  vtu_writer->EncodeAppendedDataOn();
  vtu_writer->Write();

  auto grid2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid2->SetPoints(points2);
  vtu_writer2->SetInputData(grid2);
  vtu_writer2->SetDataModeToAppended();
  vtu_writer2->EncodeAppendedDataOn();
  vtu_writer2->Write();
}

} // namespace darcy3d

int main(int argc, char *argv[]) {
  auto sim_begin = std::chrono::steady_clock::now();

  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
  options.add_options()("input-file", "path to the input file",
                        cxxopts::value<std::string>()->default_value(""))                                       //
    ("final-time", "final simulation time", cxxopts::value<double>()->default_value("1."))                      //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                            //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.1"))                                  //
    ("hyd-cond-tiss", "hydraulic conductivity in tissue domain", cxxopts::value<double>()->default_value("0.1"))                       //
    ("hyd-cond-vasc", "hydraulic conductivity in vascular domain", cxxopts::value<double>()->default_value("1."))                       //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/darcy_test_tissue_mesh.e"))                            //
    ("vasc-nifti-file", "nifti file to inform vascular subdomain", cxxopts::value<std::string>()->default_value("data/darcy_test_vascular_domain.nii.gz"))                            //
    ("voxel-size", "voxel size used in creating the tissue mesh", cxxopts::value<double>()->default_value("1."))                            //
    ("output-dir", "directory for the output", cxxopts::value<std::string>()->default_value("./output_darcy_hetero/")) //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc
  auto args = options.parse(argc, argv);
  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }
  if (!args.unmatched().empty()) {
    std::cout << "The following arguments were unmatched: " << std::endl;
    for (auto &it : args.unmatched())
      std::cout << " " << it;
    std::cout
      << "\nAre they part of petsc or a different auxiliary library?"
      << std::endl;
  }

  auto filename = args["input-file"].as<std::string>();
  auto out_dir = args["output-dir"].as<std::string>();
  auto vasc_nifti_filename = args["vasc-nifti-file"].as<std::string>();
  auto voxel_size = args["voxel-size"].as<double>();

  // read input parameters
  auto input = darcy3d::InputDeck(filename);
  if (filename == "") {
    input.d_T = args["final-time"].as<double>();
    input.d_dt = args["time-step"].as<double>();
    input.d_h = args["mesh-size"].as<double>();
    input.d_K_vasc = args["hyd-cond-vasc"].as<double>();
    input.d_K_tiss = args["hyd-cond-tiss"].as<double>();
    input.d_mesh_file = args["mesh-file"].as<std::string>();
  }

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  log("input data \n" + input.print_str() + "\n");

  // create mesh
  log("creating mesh\n");
  lm::ReplicatedMesh mesh(*comm);
  long N = long(1. / input.d_h);
  if (input.d_mesh_file != "") {
    log("reading mesh from file\n");
    mesh.read(input.d_mesh_file);
  }
  else {
    log("unit cube mesh\n");
    lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                          1., 0., 1., lm::HEX8);
  }

  // create equation system
  log("creating equation system\n");
  lm::EquationSystems eq_sys(mesh);
  eq_sys.parameters.set<darcy3d::InputDeck *>("input_deck") = &input;
  eq_sys.parameters.set<lm::Real>("time_step") = input.d_dt;
  auto &pres = eq_sys.add_system<lm::TransientLinearImplicitSystem>("P");
  pres.add_variable("p", lm::FIRST);

  // create spatial field of hydraulic conductivity
  auto &hyd_cond = eq_sys.add_system<lm::ExplicitSystem>("K");
  hyd_cond.add_variable("k", lm::CONSTANT, lm::MONOMIAL);

  // create model that holds all essential variables
  log("creating model\n");
  auto model = darcy3d::Model(comm, input, mesh, eq_sys, pres, hyd_cond, log);
  model.d_dt = input.d_dt;

  eq_sys.init();

  // create heterogeneous property field
  log("setting up K field\n");
  auto bbox = lm::MeshTools::create_bounding_box(mesh);
  auto xc = 0.5 * bbox.min() + 0.5 * bbox.max();
  auto l = (bbox.min() - bbox.max()).norm();
  eq_sys.parameters.set<lm::Point>("center") = xc;
  eq_sys.parameters.set<double>("length") = l;

  darcy3d::create_heterogeneous_conductivity(vasc_nifti_filename, voxel_size, mesh, hyd_cond, eq_sys);

  // write
  log("writing to file\n");
  mc::VTKIO(mesh).write_equation_systems(out_dir + "output_0.pvtu", eq_sys);
  return 0;

  // time stepping
  do {
    // Prepare time step
    model.d_step++;
    model.d_time += model.d_dt;

    auto solve_clock = std::chrono::steady_clock::now();
    model.d_pres.solve();
    log("solve time = " + std::to_string(mc::time_diff(solve_clock, std::chrono::steady_clock::now())) + "\n");

    if (model.d_step % 20 == 0) {
      log("writing to file\n");
      mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(model.d_step/20) + ".pvtu", eq_sys);
    }
  } while (model.d_time < input.d_T);

  // write
  log("writing to file\n");
  mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(model.d_step/20) + ".pvtu", eq_sys);

  return 0;
}

// define assembly functions
void darcy3d::Pres::assemble() {
  auto &eq_sys = d_model_p->get_system();
  const auto &input = eq_sys.parameters.get<darcy3d::InputDeck *>("input_deck");
  const auto &xc = eq_sys.parameters.get<lm::Point>("center");
  const auto &l = eq_sys.parameters.get<double>("length");
  const double dt = d_model_p->d_dt;
  const double t = d_model_p->d_time;
  auto &hyd_cond = d_model_p->d_hyd_cond;
  std::vector<unsigned int> dof_indices;

  lm::Point source_xc = xc + lm::Point(0.1 * l, 0., 0.2 * l);

  // assemble
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    // init dof map
    init_dof(elem);
    hyd_cond.get_dof_map().dof_indices(elem, dof_indices);

    // init fe
    init_fe(elem);

    // get K at this element
    double elem_K = hyd_cond.current_solution(dof_indices[0]);

    double rhs = 0.;
    auto x = elem->centroid();
    if ((x - source_xc).norm() < 0.15 * l)
      rhs = std::sin(t / 0.1) * std::exp(-t / 10.);

    for (unsigned int qp = 0; qp < d_qrule.n_points(); qp++) {
      double lhs = d_JxW[qp] * elem_K;

      // Assembling matrix
      for (unsigned int i = 0; i < d_phi.size(); i++) {
        d_Fe(i) += d_JxW[qp] * rhs * d_phi[i][qp];

        for (unsigned int j = 0; j < d_phi.size(); j++)
          d_Ke(i, j) += lhs * d_dphi[j][qp] * d_dphi[i][qp];
      }
    } // loop over quad points

    // dirichlet bc
    {
      // The penalty value.
      const double penalty = 1.e10;

      for (auto s : elem->side_index_range())
        if (elem->neighbor_ptr(s) == nullptr) {
          init_face_fe(elem, s);

          for (unsigned int qp = 0; qp < d_qrule_face.n_points(); qp++) {
            // Matrix contribution
            for (unsigned int i = 0; i < d_phi_face.size(); i++) {
              d_Fe(i) += penalty * d_JxW_face[qp] * d_phi_face[i][qp];
              for (unsigned int j = 0; j < d_phi_face.size(); j++)
                d_Ke(i, j) += penalty * d_JxW_face[qp] * d_phi_face[i][qp] * d_phi_face[j][qp];
            }
          }
        }
    }

    d_dof_map_sys.heterogenously_constrain_element_matrix_and_vector(d_Ke, d_Fe,
                                                                     d_dof_indices_sys);
    d_sys.matrix->add_matrix(d_Ke, d_dof_indices_sys);
    d_sys.rhs->add_vector(d_Fe, d_dof_indices_sys);
  } // elem loop

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}