#ifndef TUMORMODELS_NETWORKVTKWRITEROLD_HPP
#define TUMORMODELS_NETWORKVTKWRITEROLD_HPP

#include "network.hpp"

namespace util {
namespace unet {
class NetworkVTKWriterOld {
public:
  NetworkVTKWriterOld(const libMesh::Parallel::Communicator *comm_p, std::string outfilename)
      : d_comm_p(comm_p),
        d_outfilename(std::move(outfilename)) {}

  void write(ListStructure<VGNode> &VGM, const unsigned int timeStep, const unsigned int rank = 0) {
    if (d_comm_p->rank() != rank)
      return;

    const auto path = getPath(timeStep, rank);

    std::fstream filevtk;
    filevtk.open(path, std::ios::out);
    filevtk << "# vtk DataFile Version 2.0" << std::endl;
    filevtk << "Network Nutrient Transport" << std::endl;
    filevtk << "ASCII" << std::endl;
    filevtk << "DATASET POLYDATA" << std::endl;

    filevtk << "POINTS " << VGM.getNumberOfNodes() << " float" << std::endl;

    std::shared_ptr<VGNode> pointer = VGM.getHead();

    while (pointer) {

      std::vector<double> coord;

      coord = pointer->coord;

      filevtk << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;

      pointer = pointer->global_successor;
    }

    int numberOfNodes = 0;

    pointer = VGM.getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      for (int i = 0; i < numberOfNeighbors; i++) {

        if (!pointer->edge_touched[i]) {

          numberOfNodes = numberOfNodes + 1;

          pointer->edge_touched[i] = true;

          pointer->neighbors[i]->markEdge(pointer->index);
        }
      }

      pointer = pointer->global_successor;
    }

    pointer = VGM.getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      for (int i = 0; i < numberOfNeighbors; i++) {

        pointer->edge_touched[i] = false;
      }

      pointer = pointer->global_successor;
    }

    int polygons = 3 * numberOfNodes;

    filevtk << " " << std::endl;
    filevtk << "LINES " << numberOfNodes << " " << polygons << std::endl;

    pointer = VGM.getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      for (int i = 0; i < numberOfNeighbors; i++) {

        if (!pointer->edge_touched[i]) {

          filevtk << "2 " << pointer->index << " " << pointer->neighbors[i]->index
                  << std::endl;

          auto length = util::dist_between_points(pointer->coord,
                                                  pointer->neighbors[i]->coord);
          //d_model_p->d_log("vessel length: " + std::to_string(length) + "\n",
          //                "debug");
          if (util::definitelyGreaterThan(length, 0.75))
            std::cout << "vessel is too long" << std::endl;

          pointer->edge_touched[i] = true;

          pointer->neighbors[i]->markEdge(pointer->index);
        }
      }

      pointer = pointer->global_successor;
    }

    pointer = VGM.getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      for (int i = 0; i < numberOfNeighbors; i++) {

        pointer->edge_touched[i] = false;
      }

      pointer = pointer->global_successor;
    }

    filevtk << " " << std::endl;
    filevtk << "CELL_DATA " << numberOfNodes << std::endl;
    filevtk << "SCALARS radii float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;

    pointer = VGM.getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      for (int i = 0; i < numberOfNeighbors; i++) {

        if (!pointer->edge_touched[i]) {

          filevtk << pointer->radii[i] << std::endl;

          pointer->edge_touched[i] = true;

          pointer->neighbors[i]->markEdge(pointer->index);
        }
      }

      pointer = pointer->global_successor;
    }

    pointer = VGM.getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      for (int i = 0; i < numberOfNeighbors; i++) {

        pointer->edge_touched[i] = false;
      }

      pointer = pointer->global_successor;
    }

    filevtk << " " << std::endl;
    filevtk << "POINT_DATA " << VGM.getNumberOfNodes() << std::endl;
    filevtk << "SCALARS pressure_1d float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;

    pointer = VGM.getHead();

    while (pointer) {

      filevtk << pointer->p_v << std::endl;

      pointer = pointer->global_successor;
    }

    filevtk << "SCALARS nutrient_1d float 1" << std::endl;
    filevtk << "LOOKUP_TABLE default" << std::endl;

    pointer = VGM.getHead();

    while (pointer) {

      filevtk << pointer->c_v << std::endl;

      pointer = pointer->global_successor;
    }
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

#endif //TUMORMODELS_NETWORKVTKWRITEROLD_HPP
