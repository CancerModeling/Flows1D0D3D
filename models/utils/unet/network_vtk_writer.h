////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////
#ifndef TUMORMODELS_NETWORK_VTK_WRITER_HPP
#define TUMORMODELS_NETWORK_VTK_WRITER_HPP

#include "network.hpp"

namespace util {
namespace unet {
class network_vtk_writer {
public:
  network_vtk_writer(const libMesh::Parallel::Communicator *comm_p, std::string outfilename)
      : d_comm_p(comm_p),
        d_outfilename(std::move(outfilename)) {}

  void write(ListStructure<VGNode> &VGM, const unsigned int timeStep, const unsigned int rank = 0) {
    if (d_comm_p->rank() != rank)
      return;

    const auto path = getPath(timeStep, rank);

    // we create stringstreams to write into different parts of our vtk file while traversing the graph concurrently
    // the part for our point data
    // each point is written several times, once for each segment
    std::stringstream pointPart;
    // the connectivity between our points
    std::stringstream linePart;
    // the radii are saved as CELLDATA wrt the lines:
    std::stringstream radiiPart;
    // the indices are saved as POINTDATA wrt the lines:
    std::stringstream indexPart;
    // the pressures are saved as POINTDATA wrt the lines:
    std::stringstream pressurePart;
    // the nutrients are saved as POINTDATA wrt the lines:
    std::stringstream nutrientPart;

    // we save the index of the last written of our network point
    std::size_t numberOfPoints = 0;
    std::size_t numberOfSegments = 0;
    std::shared_ptr<VGNode> nodePtr = VGM.getHead();
    while (nodePtr) {
      const int numberOfNeighbors = nodePtr->neighbors.size();
      const auto &coordCurNode = nodePtr->coord;

      for (int i = 0; i < numberOfNeighbors; i++) {
        std::shared_ptr<VGNode> neighborPtr = nodePtr->neighbors[i];

        // we just write the line segment if the index of nodePtr is smaller than the index of neighborPtr
        // otherwise we would get the line segments several times.
        if (nodePtr->index >= neighborPtr->index)
          continue;

        const auto &coordNeighbor = neighborPtr->coord;

        // write the two points into our file
        pointPart << coordCurNode[0] << " " << coordCurNode[1] << " " << coordCurNode[2] << std::endl;
        pointPart << coordNeighbor[0] << " " << coordNeighbor[1] << " " << coordNeighbor[2] << std::endl;

        // connect the two points
        linePart << "2 " << numberOfPoints << " " << numberOfPoints + 1 << std::endl;

        // write the radius
        radiiPart << nodePtr->radii[i] << std::endl;

        // write the pressure at the nodes
        pressurePart << nodePtr->p_v << std::endl;
        pressurePart << neighborPtr->p_v << std::endl;

        // write the nutrients at the nodes
        nutrientPart << nodePtr->c_v << std::endl;
        nutrientPart << neighborPtr->c_v << std::endl;

        // write the indices of the nodes
        indexPart << nodePtr->index << std::endl;
        indexPart << neighborPtr->index << std::endl;

        // we added to points for the segment
        numberOfPoints += 2;

        // one new segment
        numberOfSegments += 1;
      }

      nodePtr = nodePtr->global_successor;
    }

    // we write the header
    std::fstream filevtk;
    filevtk.open(path, std::ios::out);
    filevtk << "# vtk DataFile Version 2.0" << std::endl;
    filevtk << "Network Nutrient Transport" << std::endl;
    filevtk << "ASCII" << std::endl;
    filevtk << "DATASET POLYDATA" << std::endl;

    filevtk << "POINTS " << numberOfPoints << " float" << std::endl;
    filevtk << pointPart.rdbuf();
    filevtk << " " << std::endl;

    filevtk << "LINES " << numberOfSegments << " " << 3 * numberOfSegments << std::endl;
    filevtk << linePart.rdbuf();
    filevtk << " " << std::endl;

    filevtk << " " << std::endl;
    filevtk << "CELL_DATA " << numberOfSegments << std::endl;
    filevtk << "SCALARS radii float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    filevtk << radiiPart.rdbuf();
    filevtk << " " << std::endl;

    filevtk << " " << std::endl;
    filevtk << "POINT_DATA " << numberOfPoints << std::endl;
    filevtk << "SCALARS pressure_1d float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    filevtk << pressurePart.rdbuf();
    filevtk << "SCALARS nutrient_1d float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    filevtk << nutrientPart.rdbuf();
    filevtk << "SCALARS index float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;
    filevtk << indexPart.rdbuf();
  }

private:
  std::string getPath(const unsigned int timeStep, const unsigned int rank) {
    std::string path = d_outfilename + "_";
    if (rank != 0)
      path += "rank" + std::to_string(rank) + "_";
    path += std::to_string(timeStep);
    path += ".vtk";
    return path;
  }

private:
  const libMesh::Parallel::Communicator *d_comm_p;
  const std::string d_outfilename;
};
} // namespace unet
} // namespace util

#endif //TUMORMODELS_NETWORK_VTK_WRITER_HPP
