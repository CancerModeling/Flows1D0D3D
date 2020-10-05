////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "utilLibs.hpp"
#include <iostream>
#include <string>
#include <vector>

/*!
 * @brief Namespace to collect all test functions
 */
namespace test {

/*!
 * @brief Namespace to collect test associated to mpi communications
 */
namespace comm {

/*!
 * @brief Run comm test
 */
void run(int argc, char **argv, Parallel::Communicator *comm);

} // namespace comm

/*!
 * @brief Namespace to collect test associated to geometrical functions
 */
namespace geom {

/*!
 * @brief Run rotation test
 */
void rotation(int argc, char **argv, Parallel::Communicator *comm);

/*!
 * @brief Generate cylinder cross section at multiple points on cylinder
 *
 * @param x0 Center point in cross-section of cylinder at s = 0
 * @param a Axis (unit vector) of cylinder
 * @param R Radius of cylinder
 * @param L Length of cylinder
 */
void cylinder(int argc, char **argv, Parallel::Communicator *comm, Point x0,
              Point a, double R, double L);

/*!
 * @brief Generate cylinder cross section at multiple points on cylinder
 *
 * @param x0 Center point in cross-section of cylinder at s = 0
 * @param a Axis (unit vector) of cylinder
 * @param R Radius of cylinder
 * @param L Length of cylinder
 */
void point_in_cylinder(int argc, char **argv, Parallel::Communicator *comm,
                       Point x0, Point a, double R, double L);

/*!
 * @brief Generate cylinder cross section at multiple points on cylinder
 *
 * @param x0 Center point in cross-section of cylinder at s = 0
 * @param a Axis (unit vector) of cylinder
 * @param R Radius of cylinder
 * @param L Length of cylinder
 */
void angle_test(int argc, char **argv, Parallel::Communicator *comm);

/*!
 * @brief Transform cube into unit ball
 *
 * @param dim Dimension
 * @param a Size of the cube
 */
void cube_to_ball(int argc, char **argv, unsigned int dim, double a,
                  Parallel::Communicator *comm);

/*!
 * @brief Discretize cylinder surface and find the elements in mesh
 * containing these points
 *
 * @param x0 Center point in cross-section of cylinder at s = 0
 * @param a Axis (unit vector) of cylinder
 * @param R Radius of cylinder
 * @param L Length of cylinder
 */
void cylinder_intersection_mesh(int argc, char **argv,
                                Parallel::Communicator *comm, Point x0, Point a,
                                double R, double L);

} // namespace geom

/*!
 * @brief Namespace to test libmesh mesh features
 */
namespace mesh {

/*!
 * @brief Run test: start with empty mesh and successively add new nodes and
 * elements
 */
void add_node_elem_test(int argc, char **argv, Parallel::Communicator *comm);

/*!
 * @brief Run test: Same as add_node_elem_test() but now we add equation
 * system for pressure
 */
void add_node_elem_eq_sys_test(int argc, char **argv,
                               Parallel::Communicator *comm);

void add_node_elem_eq_sys_test_2(int argc, char **argv,
                                 Parallel::Communicator *comm);

/*!
 * @brief Run test: Same as add_node_elem_eq_sys_test() but now we add
 * aseembly function for simple 1-d flow equation on network and solve the
 * 1-d equation on changing mesh
 *
 * In this test, we start with empty mesh and do following
 *
 * 1. Add node at (0,0,0) and (1,1,1)
 * 2. Create edge element between above two nodes
 * 3. Assign nodes and element specific boundary ids
 * 4. Create a equation system for pressure on mesh with nodes and element as
 * above
 * 5. Set pressure dof at node 1 (i.e. at point (0,0,0)) to specific value
 * (Dirichlet boundary condition on this node) and set node 2 pressure dof
 * free.
 * 6. Assemble matrix and vector. Weak form is given by
 *
 * \int_0^L a(s) u v ds = \int_0^L f(s) v(s) ds
 *
 * where u is trial function, v is test function, and s is parametric
 * coordinate of point on line connecting (0,0,0) and (1,1,1).
 *
 * We let f(s) = 1 and a(s) = 1.
 *
 * Assuming linear interpolation, the exact numerical solution for two node
 * one element mesh can be obtained
 *
 * p(s2) = (3 - p(s1)) / 2.
 *
 * Where p(s1) is the Dirichlet boundary condition. We choose some value of p
 * (s1) and compare the numerical result p(s2) with above exact result.
 */
void eq_sys_assemble(int argc, char **argv, Parallel::Communicator *comm);

/*!
 * @brief Run test: Similar to eq_sys_assemble() except now the weak form on
 * each line element is
 *
 * \f[\int_0^L a(s) du/ds dv/ds ds = \int_0^L f(s) v(s) ds \f]
 *
 * where u is trial function, v is test function.
 *
 * We let f(s) = 1 and a(s) = 1.
 *
 * Let v1 = (0, 0, 0),
 *  v2 = (1, 1, 1),
 *  v3 = (1.5, 1 - \f$ \sqrt{5/2} \f$, 0.5),
 *  v4 = (2, 1 - 2\f$ \sqrt{5/2} \f$, 1).
 *
 *  Distance between v1-v2, v2-v3, v3-v4 are \f$ L = \sqrt{3} \f$.
 *
 * We consider three tests:
 *
 * - In first step, we consider single segments v1-v2, and compute the
 * numerical solution p2 and compare with the following numerical exact
 * solution:
 *
 * p2 = L^2/2 + p_in
 *
 * - In second step, we add one segment v2-v3 to the system, and compute the
 * numerical solution p2, p3 and compare with the following numerical exact
 * solution:
 *
 * p2 = 3L^2/2 + p_in, p3 = 3L^2 + p_in
 *
 * - In third step, we add another segment v3-v4 to the system and compute
 * numerical solution p2, p3, p4 and compare with the following numerical
 * exact solution:
 *
 *  p2 = 5L^2/2 + p_in, p3 = 4L^2 + p_in, p3 = 9L^2/2 + p_in.
 *
 * @param p_in Dirichlet bc at first vertex
 */
void eq_sys_assemble_2(int argc, char **argv, double p_in,
                       Parallel::Communicator *comm);

/*!
 * @brief Run test: Similar to eq_sys_assemble_2() with only third step
 * changed to following:
 *
 * Let v1, v2, v3 same as in eq_sys_assemble_2(). Let v4 is
 *
 * v4 = (0.5, 1 + \f$ \sqrt{5/2} \f$, 1.5)
 *
 * - In third step, we add segment v2-v4. Note unlike in eq_sys_assemble_2()
 * test, here the the result is that there is a branch at v2, i.e. it is
 * connected to three segments. We compute numerical solution p2, p3, p4 and
 * compare with the following numerical exact solution:
 *
 *  p2 = 5L^2/2 + p_in, p3 = 3L^2 + p_in, p3 = 3L^2 + p_in.
 *
 * @param p_in Dirichlet bc at first vertex
 */
void eq_sys_assemble_3(int argc, char **argv, double p_in,
                       Parallel::Communicator *comm);

/*!
 * @brief Run test: We now consider transient equation and compare the
 * solution. The weak form of the equation is
 *
 * \f[\int du/dt v ds + \int a(s) du/ds dv/ds ds = \int f(s) v(s) ds \f]
 *
 * We compare the results of explicit forward-Euler discretization in three
 * cases: in first case we have single segment, in second case we add one
 * segment to the system, and in third case we add another segment so that
 * second vertex becomes the branching vertex.
 *
 * Let v1, v2, v3 same as in eq_sys_assemble_2(). Let v4 is
 *
 * v4 = (0.5, 1 + \f$ \sqrt{5/2} \f$, 1.5)
 *
 * @param p_in Dirichlet bc at first vertex
 */
void eq_sys_assemble_transient(int argc, char **argv, std::string test_case,
                               unsigned int num_steps,
                               double p_in, Parallel::Communicator *comm);

void elem_id_numbering(int argc, char **argv, Parallel::Communicator *comm);

void memory_leak(int argc, char **argv, Parallel::Communicator *comm);

} // namespace mesh

/*!
 * @brief Namespace to test finite volume method
 */
namespace fv {

/*!
 * @brief Run test: start with empty mesh and successively add new nodes and
 * elements
 */
void run(int argc, char **argv, Parallel::Communicator *comm,
         double pressure_in, bool branch);

/*!
 * @brief Check deletion of element, adding new element at mid point etc
 */
void add_and_delete(int argc, char **argv, Parallel::Communicator *comm,
                    double pressure_in, bool branch);

} // namespace fv

} // namespace test
