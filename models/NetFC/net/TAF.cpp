////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "netUtil.hpp"
#include "network.hpp"
#include "nodes.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>

void netfc::Network::assemble3DSystemForTAF() {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  // 3D-1D TAF problem on a cube
  std::cout << " " << std::endl;
  std::cout << "3D TAF problem on a cube \\Omega = (0," << L_x << ")^3" << std::endl;

  double vol_elem = h_3D * h_3D * h_3D;

  double area_face = h_3D * h_3D;

  double K_3D = input.d_tissue_flow_K;

  for (int i = 0; i < A_TAF_3D.nrows(); i++) {

    A_TAF_3D[i].clear();
  }

  for (int i = 0; i < b_TAF_3D.size(); i++) {

    b_TAF_3D[i] = 0.0;
  }

  std::vector<std::vector<double>> directions;

  directions = defineDirections();

  double dt = d_model_p->d_dt;
  std::cout << "dt: " << dt << std::endl;

  double K_G = 1.0e-3;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        A_TAF_3D(index, index) += vol_elem;

        b_TAF_3D[index] = vol_elem * phi_TAF_old[index];

        // Get element center
        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        b_TAF_3D[index] = b_TAF_3D[index] + vol_elem * dt * sourceTermTAFTwoVessels(center);

        // Iterate over the interfaces
        for (int face = 0; face < 6; face++) {

          std::vector<double> center_neighbor = getCenterNeighbor(center, directions[face], h_3D);

          bool isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            int index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            double v = -K_3D * (P_3D1D[index_neighbor] - P_3D1D[index]) / h_3D;

            if (v > 0.0) {

              A_TAF_3D(index, index) += dt * area_face * v;

            } else {

              A_TAF_3D(index, index_neighbor) += dt * area_face * v;
            }

            A_TAF_3D(index, index) += dt * D_TAF * area_face / h_3D + vol_elem * dt * K_G;

            A_TAF_3D(index, index_neighbor) = -dt * D_TAF * area_face / h_3D;
          }
        }
      }
    }
  }
}

void netfc::Network::solve3DTAFProblem(int timeStep, double time) {

  // Solver
  gmm::iteration iter(5.0e-10);

  std::cout << " " << std::endl;
  std::cout << "Assemble 3D TAF matrix and right hand side" << std::endl;
  assemble3DSystemForTAF();

  phi_TAF = phi_TAF_old;

  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_TAF_3D, 50, 1e-6);

  std::cout << " " << std::endl;
  std::cout << "Solve linear system of equations (TAF)" << std::endl;
  gmm::bicgstab(A_TAF_3D, phi_TAF, b_TAF_3D, P, iter);

  phi_TAF_old = phi_TAF;

  if (timeStep % 2 == 0) {

    std::cout << " " << std::endl;
    std::cout << "Plot solutions" << std::endl;
    writeDataToVTK3D_TAF(phi_TAF, N_3D, h_3D, timeStep);
  }
}

double netfc::Network::sourceTermTAFTwoVessels(std::vector<double> coord) {

  double source_TAF = 0.0;

  std::vector<double> center_1;

  center_1.push_back(0.5);
  center_1.push_back(1.5);
  center_1.push_back(1.0);

  std::vector<double> dist_vec_1 = std::vector<double>(3, 0.0);

  for (int i = 0; i < 3; i++) {

    dist_vec_1[i] = center_1[i] - coord[i];
  }

  double dist_1 = gmm::vect_norm2(dist_vec_1);

  std::vector<double> center_2;

  center_2.push_back(1.5);
  center_2.push_back(0.5);
  center_2.push_back(1.0);

  std::vector<double> dist_vec_2 = std::vector<double>(3, 0.0);

  for (int i = 0; i < 3; i++) {

    dist_vec_2[i] = center_2[i] - coord[i];
  }

  double dist_2 = gmm::vect_norm2(dist_vec_2);

  int element_index = getElementIndex(coord, h_3D, N_3D);

  double nutrients = phi_sigma[element_index];

  if (dist_1 < 0.2) {

    if (nutrients < 0.00005) {

      source_TAF = 0.75;

    } else if (nutrients > 0.00005 && nutrients < 0.01) {

      source_TAF = 0.75 / 8.0 * (1.0 + 7.0 * std::pow((0.01 - nutrients) / (0.01 - 0.00005), 3.0));

    } else {

      source_TAF = 0.75 / 8.0;
    }
  }

  if (dist_2 < 0.2) {

    if (nutrients < 0.00005) {

      source_TAF = 0.75;

    } else if (nutrients > 0.00005 && nutrients < 0.01) {

      source_TAF = 0.75 / 8.0 * (1.0 + 7.0 * std::pow((0.01 - nutrients) / (0.01 - 0.00005), 3.0));

    } else {

      source_TAF = 0.75 / 8.0;
    }
  }

  return source_TAF;
}

void netfc::Network::writeDataToVTK3D_TAF(std::vector<double> phi_TAF_3D, int N_3D, double h_3D, int timeStep) {

  int numberOfCells = phi_TAF_3D.size();

  int numberOfPolygonData = 9 * numberOfCells;

  int numberOfPoints = (N_3D + 1) * (N_3D + 1) * (N_3D + 1);

  std::string path = scenario;
  path += "_TAF_3D_";
  path += std::to_string(timeStep);
  path.append(".vtk");

  std::fstream filevtk;
  filevtk.open(path, std::ios::out);
  filevtk << "# vtk DataFile Version 2.0" << std::endl;
  filevtk << "Pressure 3D1D coupled problem" << std::endl;
  filevtk << "ASCII" << std::endl;
  filevtk << "DATASET UNSTRUCTURED_GRID" << std::endl;

  filevtk << "POINTS " << numberOfPoints << " float" << std::endl;

  std::vector<std::vector<int>> indices;

  for (int i = 0; i < N_3D + 1; i++) { // x-loop

    for (int j = 0; j < N_3D + 1; j++) { // y-loop

      for (int k = 0; k < N_3D + 1; k++) { // z-loop

        filevtk << (double) i * h_3D << " " << (double) j * h_3D << " " << (double) k * h_3D << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        filevtk << "8"
                << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] + 0.5 * h_3D, h_3D)
                << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_TYPES " << numberOfCells << std::endl;

  for (int i = 0; i < numberOfCells; i++) {

    filevtk << 11 << std::endl;
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_DATA " << numberOfCells << std::endl;
  filevtk << "SCALARS Phi_TAF_(3D) float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        filevtk << phi_TAF_3D[index] << std::endl;
      }
    }
  }
}
