////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tests.hpp"
#include <libmesh/getpot.h>
#include <iostream>

int main(int argc, char* argv[]){

  // Open file with model setup
  std::string filename = "input.in";

  // Open file with model setup
  GetPot input_file(filename);

  // find which model to run
  std::string model_name = input_file("model_name", "AvaLima");

  // Initialize libMesh
  LibMeshInit init(argc, argv);

  if (model_name == "TestComm") {

    // run mpi comm test
    test::comm::run(argc, argv, &init.comm());
  } else if (model_name == "TestGeom") {

    // run geom test
    test::geom::rotation(argc, argv, &init.comm());

    test::geom::angle_test(argc, argv, &init.comm());

    test::geom::cylinder(argc, argv, &init.comm(), Point(5., 0., 0.),
                         Point(1., 1., -1.), 4, 10);

    test::geom::point_in_cylinder(argc, argv, &init.comm(), Point(5., 0., 0.),
                         Point(1., 1., -1.), 4, 10);

    test::geom::cube_to_ball(argc, argv, 2, 2., &init.comm());

    test::geom::cube_to_ball(argc, argv, 3, 2., &init.comm());

    test::geom::cylinder_intersection_mesh(argc, argv, &init.comm(), Point(5., 0., 0.),
                         Point(1., 1., -1.), 4, 10);
  } else if (model_name == "TestMesh") {

    // run mesh test
    if (false) {
    test::mesh::add_node_elem_test(argc, argv, &init.comm());

    test::mesh::add_node_elem_eq_sys_test_2(argc, argv, &init.comm());

    test::mesh::eq_sys_assemble(argc, argv, &init.comm());

    test::mesh::eq_sys_assemble_2(argc, argv, 0.5, &init.comm());

    test::mesh::eq_sys_assemble_3(argc, argv, 1.5, &init.comm());

    //    test::mesh::eq_sys_assemble_transient(argc, argv, "one_segment",
    //                                          10, 1.5, &init.comm());
    //
    //    test::mesh::eq_sys_assemble_transient(argc, argv, "two_segments",
    //                                          8, 0.9, &init.comm());

    test::mesh::eq_sys_assemble_transient(argc, argv, "three_segments_branch",
                                          4, 3.7, &init.comm());
  }

    test::mesh::elem_id_numbering(argc, argv, &init.comm());

    // test::mesh::memory_leak(argc, argv, &init.comm());
    
  } else if (model_name == "TestFV") {

    test::fv::run(argc, argv, &init.comm(), 1., true);

    test::fv::add_and_delete(argc, argv, &init.comm(), 1., true);
  }


  // End Code
  return 0;
}
