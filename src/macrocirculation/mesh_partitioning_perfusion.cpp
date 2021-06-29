////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "mesh_partitioning_perfusion.hpp"
#include "assembly_system.hpp"
#include "base_model.hpp"
#include "random_dist.hpp"
#include "tree_search.hpp"
#include "vtk_io.hpp"

#include <fmt/format.h>

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <memory>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace mc = macrocirculation;

void macrocirculation::estimate_flow_rate(const std::vector<double> &radii,
                                          double gamma,
                                          std::vector<double> &flow) {
  flow.clear();
  for (size_t i = 0; i < radii.size(); i++)
    flow.push_back(std::pow(radii[i], gamma));
}

void macrocirculation::estimate_vol_fraction(const std::vector<double> &radii,
                                             const std::vector<double> &flow,
                                             std::vector<double> &vol) {
  vol.clear();

  double flow_sum = 0.;
  for (auto q : flow) flow_sum += q;
  for (int i = 0; i < flow.size(); i++)
    vol.push_back(flow[i] / flow_sum);
}

void macrocirculation::create_partition_field(const lm::MeshBase &mesh,
                                              lm::ExplicitSystem &field,
                                              const std::vector<long> &partition) {

  std::vector<unsigned int> dof_indices;
  for (const auto &elem : mesh.active_local_element_ptr_range()) {
    field.get_dof_map().dof_indices(elem, dof_indices);
    int outlet_id = partition[elem->id()];
    field.solution->set(dof_indices[0], double(outlet_id+1));
  }
  field.solution->close();
  field.update();
}

void macrocirculation::create_perfusion_territory(const std::vector<lm::Point> &x,
                                                  const std::vector<double> &v,
                                                  lm::ReplicatedMesh &mesh,
                                                  std::vector<long> &partition,
                                                  int n_rand_loop,
                                                  lm::ExplicitSystem &partition_field,
                                                  lm::EquationSystems &eq_sys,
                                                  macrocirculation::Logger &log,
                                                  std::string out_dir,
                                                  int out_freq) {

  // initialize random number generator
  int seed = 0;
  srand(seed);

  // create list of element centers for tree search
  int nelems = mesh.n_elem();
  std::vector<lm::Point> elem_centers(nelems, lm::Point());
  double total_vol = 0.;
  for (const auto &elem : mesh.element_ptr_range()) {
    elem_centers[elem->id()] = elem->centroid();
    total_vol += elem->volume();
  }

  // create tree for search
  std::unique_ptr<mc::NFlannSearchKd> tree = std::make_unique<mc::NFlannSearchKd>(elem_centers);
  tree->set_input_cloud();

  // step 1 - sort the points based on increasing volume
  log("create_perfusion_territory() log: \n");
  log("  sorting points\n");
  auto imap = mc::sort_indexes(v);

  // print
  {
    log("    info before sorting\n");
    std::ostringstream oss;
    std::vector<size_t> temp;
    for (size_t i=0; i<x.size(); i++) temp.push_back(i);
    oss << "    pts ids = " << mc::print_str(temp);
    oss << "    vol fraction target = " << mc::print_str(v);
    log(oss.str());
  }

  // rearrange
  std::vector<lm::Point> pts = x;
  std::vector<double> vol = v; // convert from vol fraction to volume
  for (size_t i = 0; i < imap.size(); i++) {
    auto ii = imap[i];
    pts[i] = x[ii];
    vol[i] = v[ii] * total_vol;
  }

  // compute volume represented by points and cube root of vol as radius estimation
  double fos = 1.25;
  std::vector<double> vol_r(pts.size(), 0.);
  for (size_t i = 0; i < vol.size(); i++)
    vol_r[i] = fos * std::pow(3. * vol[i] / (4. * M_PI), 1. / 3.);

  // compute min vol and vol tolerance
  double min_vol = vol[0] * total_vol;
  double vol_tol = 0.01 * min_vol;

  // print
  {
    log("    info after sorting\n");
    std::ostringstream oss;
    oss << "    pts ids = " << mc::print_str(imap);
    oss << "    vol target = " << mc::print_str(vol);
    oss << "    radius target = " << mc::print_str(vol_r);
    log(oss.str());
  }

  // step 2 - brute force element assigning
  log("  brute-force assignment\n");
  std::vector<std::vector<long>> pts_elems(pts.size(), std::vector<long>());
  std::vector<double> actual_vol(pts.size(), 0.);

  // resize partition data and fill with invalid outlet id
  partition.resize(mesh.n_elem());
  for (size_t i = 0; i < mesh.n_elem(); i++)
    partition[i] = -1;

  // assign elements
  for (int i = 0; i < pts.size(); i++) {
    auto xi = pts[i];
    auto voli = vol[i];
    auto voli_r = vol_r[i];

    // find elements whose center is within radi distance of outlet point
    std::vector<size_t> neighs;
    std::vector<double> sqr_dist;
    auto search_status =
      tree->radius_search(xi, voli_r, neighs, sqr_dist);

    log("    search_status = " + std::to_string(search_status) + "\n");

    if (search_status == 0) {
      std::cerr << "Error: Not enough elements in the neighborhood of outlet point\n";
      exit(EXIT_FAILURE);
    }

    double sum_voli = 0.;
    for (int j = 0; j < neighs.size(); j++) {
      auto j_id = neighs[j];

      // check if this element is already taken
      if (partition[j_id] != -1)
        continue;

      auto elemj_vol = mesh.elem_ptr(j_id)->volume();
      if (sum_voli + elemj_vol < voli + vol_tol) {

        pts_elems[i].push_back(j_id);
        sum_voli += elemj_vol;

        // store the fact that j_id element now is owned by outlet point i
        partition[j_id] = i;
      }
    }
    actual_vol[i] = sum_voli;
  }

  // update conductivity parameter and write to file
  create_partition_field(mesh, partition_field, partition);
  log("  writing to file\n");
  mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(0) + ".pvtu", eq_sys);

  // step 3 - randomly readjust the element distribution
  log("  redistribution\n");
  std::vector<int> pts_id_vec;
  for (int i = 0; i < pts.size(); i++) pts_id_vec.push_back(i);

  std::vector<double> store_err;
  std::vector<double> count_no_outlet_elem;
  std::vector<double> vol_not_dist;
  std::vector<double> vol_err(pts.size(), 0.);
  int out_count = 1;
  for (int rloop = 0; rloop < n_rand_loop; rloop++) {
    log("    random loop = " + std::to_string(rloop) + "\n");

    // shuffle pts_id_vec
    std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
    std::vector<int> v(pts_id_vec.begin(), pts_id_vec.end());
    std::shuffle(v.begin(), v.end(), gen);
    pts_id_vec.assign(v.begin(), v.end());

    std::ostringstream oss;
    oss << "      shuffled list = " << mc::print_str(pts_id_vec);

    // compute errors
    double avg_err = 0.;
    int no_out_el = 0;
    double total_vol_not_dist = 0.;
    {
      // get volume percentage error
      for (int i = 0; i < pts.size(); i++) {
        vol_err[i] = 100. * std::abs(vol[i] - actual_vol[i]) / vol[i];
        avg_err += vol_err[i];
      }
      avg_err = avg_err / pts.size();
      store_err.push_back(avg_err);

      // get number of element without any outlet
      for (auto i : partition) {
        if (i == -1)
          no_out_el++;
      }
      count_no_outlet_elem.push_back(no_out_el);

      // get total volume which is not assigned
      for (auto i : actual_vol) total_vol_not_dist += i;
      total_vol_not_dist = 100. * (total_vol - total_vol_not_dist) / total_vol;
      vol_not_dist.push_back(total_vol_not_dist);

      // debug
      std::vector<size_t> temp;
      temp.clear();
      for (const auto &i : pts_elems) temp.push_back(i.size());
      oss << "      num elements = " << mc::print_str(temp);

      oss << "      vol actual = " << mc::print_str(actual_vol);
      oss << "      vol target = " << mc::print_str(vol);
      oss << "      radius target = " << mc::print_str(vol_r);

      oss << "      percent vol difference = " << mc::print_str(vol_err);
      oss << "      average percent error = " << avg_err << "\n";
      oss << "      number of no outlet elems = " << no_out_el << "\n";
      oss << "      percent volume not distributed = " << total_vol_not_dist << "\n";
      log(oss.str());
    }

    // loop over outlet points
    int num_pts_change = 0;
    int num_elem_moved = 0;
    for (int ii = 0; ii < pts_id_vec.size(); ii++) {
      int i = pts_id_vec[ii];

      auto xi = pts[i];
      auto voli = vol[i];
      auto radi = vol_r[i];
      auto voli_act = actual_vol[i];

      std::vector<long> elem_i_old = pts_elems[i];
      auto &elem_i = pts_elems[i];

      // check if we really need to process this outlet
      if (voli_act < voli + vol_tol and voli_act > voli - vol_tol)
        continue;

      log(fmt::format("      processing outlet end = {}\n", i));

      // loop over elements of this outlet and for each element find the neighboring element
      int num_neigh_elem = 10;
      double sum_voli_act = voli_act;
      for (auto e : elem_i_old) {
        auto xe = elem_centers[e];

        // if element e is very far from the point center, remove it from the list
        if ((xe - xi).norm() > 5 * radi) {
          auto vol_e = mesh.elem_ptr(e)->volume();
          elem_i.erase(std::find(elem_i.begin(), elem_i.end(), e));
          sum_voli_act -= vol_e;

          // free this element
          partition[e] = -1;

          continue;
        }

        std::vector<size_t> ei_neighs;
        std::vector<double> ei_sqr_dist;
        tree->nearest_search(xe, num_neigh_elem, ei_neighs, ei_sqr_dist);

        // loop over neighboring elements
        for (auto ee : ei_neighs) {
          auto vol_ee = mesh.elem_ptr(ee)->volume();
          auto ee_taken = partition[ee];

          // check if ee already exists in the list
          if (ee_taken == i)
            continue;

          // check if adding this element perturbs the volume a lot
          if (sum_voli_act + vol_ee > voli + vol_tol)
            continue;

          // move element ee from previous owner to outlet i
          pts_elems[i].push_back(ee);
          auto ee_outlet_id_old = partition[ee];
          partition[ee] = i;
          sum_voli_act += vol_ee;

          // remove ee from old owner
          if (ee_outlet_id_old != -1) {
            auto &elem_j = pts_elems[ee_outlet_id_old];
            elem_j.erase(std::find(elem_j.begin(), elem_j.end(), ee));
            actual_vol[ee_outlet_id_old] -= vol_ee;
          }
        }
      }

      actual_vol[i] = sum_voli_act;

      if (elem_i_old.size() < pts_elems[i].size()) {
        num_pts_change++;
        num_elem_moved += pts_elems[i].size() - elem_i_old.size();
      }
    }

    log(fmt::format("      num_pts_change = {}, num_elem_moved = {}\n", num_pts_change, num_elem_moved));

    if (rloop % out_freq == 0) {
      // update conductivity parameter and write to file
      create_partition_field(mesh, partition_field, partition);
      log("      writing to file\n");
      mc::VTKIO(mesh).write_equation_systems(out_dir + "output_" + std::to_string(out_count) + ".pvtu", eq_sys);
      out_count++;
    }
  }

  {
    std::ofstream oss;
    oss.open(out_dir + "average_error.txt");
    for (int i=0; i<store_err.size(); i++)
      oss << store_err[i] << ", "
          << count_no_outlet_elem[i] << ", "
          << vol_not_dist[i] << "\n";
    oss.close();
  }
}