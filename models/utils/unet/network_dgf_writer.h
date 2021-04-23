////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////
#ifndef TUMORMODELS_NETWORK_DGF_WRITER_HPP
#define TUMORMODELS_NETWORK_DGF_WRITER_HPP

#include "network.hpp"
#include <cassert>

namespace util {
namespace unet {

/*!
 * @brief Writes the current 1D network into a dgf file, which we can read back into our simulation.
 *        Note: Only the at the d
 */
class network_dgf_writer {
public:
  network_dgf_writer(const libMesh::Parallel::Communicator *comm_p, std::string outfilename)
      : d_comm_p(comm_p),
        d_outfilename(std::move(outfilename)) {}

  /*!
   * @brief Writes the 1D network to disk, overwriting the previous results.
   *
   * @param VGM
   */
  void write(ListStructure<VGNode> &VGM) {
    if (d_comm_p->rank() != 0)
      return;

    std::string path = d_outfilename + "_last.dgf";

    write(VGM, path);
  }

private:
  void write(ListStructure<VGNode> &VGM, const std::string & path) {
    std::fstream filedgf;
    filedgf.open(path, std::ios::out);
    filedgf << "DGF" << std::endl;

    // write vertex data
    {
      filedgf << "Vertex" << std::endl;
      filedgf << "parameters 2" << std::endl;
      std::shared_ptr<VGNode> pointer = VGM.getHead();
      uint idx = 0;   // indices currently start at 1.
      while (pointer) {
        // make sure that our assumptions about the order of vertices are correct.
        assert(idx == pointer->index);

        const auto coord = pointer->coord;

        // if our node is dirichlet, then we export the real pressure and nutrient values,
        // otherwise we set them to zero.
        const auto pressure = pointer->typeOfVGNode == DirichletNode ? pointer->p_v : 0;
        const auto nutrient = pointer->typeOfVGNode == DirichletNode ? pointer->c_v : 0;

        filedgf << std::scientific << std::setprecision(16);
        filedgf << coord[0] << " " << coord[1] << " " << coord[2] << " " << pressure << " " << nutrient << std::endl;

        pointer = pointer->global_successor;
        idx += 1;
      }
    }

    filedgf << "#" << std::endl;

    // write simplex data
    {
      filedgf << "SIMPLEX" << std::endl;
      std::shared_ptr<VGNode> nodePtr = VGM.getHead();
      filedgf << "parameters 1" << std::endl;
      while (nodePtr) {

        for (int i = 0; i < nodePtr->neighbors.size(); i++) {
          std::shared_ptr<VGNode> neighborPtr = nodePtr->neighbors[i];

          // we just write the line segment if the index of nodePtr is smaller than the index of neighborPtr
          // otherwise we would get the line segments several times.
          if (nodePtr->index >= neighborPtr->index)
            continue;

          filedgf << nodePtr->index << " " << neighborPtr->index << " " << nodePtr->radii[i] << std::endl;
        }

        nodePtr = nodePtr->global_successor;
      }
    }

    filedgf << "#" << std::endl;
  }

private:
  const libMesh::Parallel::Communicator *d_comm_p;
  const std::string d_outfilename;
};
} // namespace unet
} // namespace util

#endif
